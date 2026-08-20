import gzip
import json
from pathlib import Path

import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
import umap


PROJECT = "project-001-pulmonary-arterial-hypertension-transcriptomics"
ROOT = Path(r"G:\workdata\projects") / PROJECT
RAW = ROOT / "raw" / "current-data-snapshot" / "pah_scrna"
BASE = RAW / "gse210248" / "Human_PA_Integrated"
OUT = ROOT / "processed" / "gse210248_harmony"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42


def neighbor_metrics(z, sample, celltype, smc_label, k=30):
    ind = NearestNeighbors(n_neighbors=k + 1).fit(z).kneighbors(return_distance=False)[:, 1:]
    same_sample = np.mean(sample[ind] == sample[:, None])
    same_celltype = np.mean(celltype[ind] == celltype[:, None])
    smc = np.isin(smc_label, ["SMC 1", "SMC 2"])
    smc_rows = np.where(smc)[0]
    smc_neighbor_labels = smc_label[ind[smc_rows]]
    valid = np.isin(smc_neighbor_labels, ["SMC 1", "SMC 2"])
    cross = valid & (smc_neighbor_labels != smc_label[smc_rows, None])
    return {
        "same_sample_knn_fraction": float(same_sample),
        "same_celltype_knn_fraction": float(same_celltype),
        "smc_cross_subtype_neighbors_all_knn": float(cross.mean()),
        "smc_neighbors_among_all_knn": float(valid.mean()),
        "smc_cross_fraction_within_smc_neighbors": float(cross.sum() / max(valid.sum(), 1)),
    }


def safe_silhouette(z, labels, sample_size=5000):
    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return None
    return float(silhouette_score(z, labels, sample_size=min(sample_size, len(labels)), random_state=SEED))


def annotate_centroids(ax, xy, labels, fontsize=6):
    frame = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "label": labels})
    for label, g in frame.groupby("label"):
        ax.text(g.x.median(), g.y.median(), str(label), fontsize=fontsize, ha="center", va="center",
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.72})


def scatter_panel(ax, xy, labels, title, palette):
    labels = np.asarray(labels)
    for label in sorted(np.unique(labels)):
        keep = labels == label
        ax.scatter(xy[keep, 0], xy[keep, 1], s=1.2, alpha=0.48, linewidths=0,
                   color=palette[label], rasterized=True)
    annotate_centroids(ax, xy, labels)
    ax.set_title(title, loc="left", fontsize=9, weight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)


def main():
    meta = pd.read_csv(RAW / "GSE210248_Human_PA_Metadata.csv.gz", sep=";", decimal=",", index_col=0)
    with gzip.open(BASE / "barcodes.tsv.gz", "rt") as handle:
        barcodes = [x.strip() for x in handle]
    with gzip.open(BASE / "features.tsv.gz", "rt") as handle:
        genes = [x.rstrip().split("\t")[1] for x in handle]
    if barcodes != list(meta.index):
        raise ValueError("Barcode order differs between matrix and metadata")

    x = mmread(str(BASE / "matrix.mtx.gz")).tocsr().T.astype(np.float32)
    lib = np.asarray(x.sum(axis=1)).ravel()
    scale = np.divide(1e4, lib, out=np.zeros_like(lib), where=lib > 0)
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    mean = np.asarray(x.mean(axis=0)).ravel()
    var = np.asarray(x.power(2).mean(axis=0)).ravel() - mean ** 2
    valid = np.array([not (g.startswith("MT-") or g.startswith("RPL") or g.startswith("RPS")) for g in genes])
    hvg = np.where(valid)[0][np.argsort(var[valid])[-2000:]]
    pcs = TruncatedSVD(n_components=30, random_state=SEED).fit_transform(x[:, hvg])

    batch_meta = pd.DataFrame({"sample": meta["new.ident"].astype(str).values})
    harmony = hm.run_harmony(pcs, batch_meta, "sample", random_state=SEED, verbose=True)
    corrected = np.asarray(harmony.Z_corr)
    if corrected.shape[0] != pcs.shape[0] and corrected.shape[1] == pcs.shape[0]:
        corrected = corrected.T
    if corrected.shape != pcs.shape:
        raise ValueError(f"Unexpected Harmony matrix shape {corrected.shape}; expected {pcs.shape}")

    reducer = dict(n_neighbors=20, min_dist=0.30, n_components=2, metric="euclidean", random_state=SEED, n_jobs=1)
    raw_umap = umap.UMAP(**reducer).fit_transform(pcs)
    harmony_umap = umap.UMAP(**reducer).fit_transform(corrected)

    sample = meta["new.ident"].astype(str).to_numpy(dtype=object)
    celltype = meta["Cell_annotation"].astype(str).to_numpy(dtype=object)
    author_cluster = meta["integrated_snn_res.0.3"].astype(str).to_numpy(dtype=object)
    disease = meta["disease"].astype(str).to_numpy(dtype=object)
    metrics = []
    for name, z in [("uncorrected_pcs", pcs), ("harmony_pcs", corrected)]:
        row = {"representation": name}
        row.update(neighbor_metrics(z, sample, celltype, celltype))
        row["sample_silhouette"] = safe_silhouette(z, sample)
        row["celltype_silhouette"] = safe_silhouette(z, celltype)
        row["author_cluster_silhouette"] = safe_silhouette(z, author_cluster)
        smc = np.isin(celltype, ["SMC 1", "SMC 2"])
        row["smc1_smc2_silhouette"] = safe_silhouette(z[smc], celltype[smc], sample_size=int(smc.sum()))
        km = KMeans(n_clusters=len(np.unique(author_cluster)), random_state=SEED, n_init=20).fit_predict(z)
        row["kmeans_vs_author_ARI"] = float(adjusted_rand_score(author_cluster, km))
        row["kmeans_vs_author_NMI"] = float(normalized_mutual_info_score(author_cluster, km))
        metrics.append(row)
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(OUT / "harmony_batch_correction_metrics.csv", index=False)

    coords = pd.DataFrame({
        "barcode": barcodes, "raw_UMAP1": raw_umap[:, 0], "raw_UMAP2": raw_umap[:, 1],
        "harmony_UMAP1": harmony_umap[:, 0], "harmony_UMAP2": harmony_umap[:, 1],
        "cell_type": celltype, "author_cluster_res0.3": author_cluster,
        "sample": sample, "disease": disease,
    })
    coords.to_csv(OUT / "gse210248_raw_vs_harmony_umap.csv.gz", index=False, compression="gzip")

    labels = sorted(np.unique(celltype))
    cmap = plt.get_cmap("tab20")
    palette = {lab: cmap(i % 20) for i, lab in enumerate(labels)}
    smc_palette = {"SMC 1": "#2B6CB0", "SMC 2": "#D97706"}
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.0), constrained_layout=True)
    scatter_panel(axes[0, 0], raw_umap, celltype, "A  Uncorrected: all annotated cell types", palette)
    scatter_panel(axes[0, 1], harmony_umap, celltype, "B  Harmony by sample: all cell types", palette)
    smc = np.isin(celltype, ["SMC 1", "SMC 2"])
    scatter_panel(axes[1, 0], raw_umap[smc], celltype[smc], "C  Uncorrected: SMC1 and SMC2", smc_palette)
    scatter_panel(axes[1, 1], harmony_umap[smc], celltype[smc], "D  Harmony: SMC1 and SMC2", smc_palette)
    fig.suptitle("GSE210248 sample-level Harmony sensitivity analysis", fontsize=11, weight="bold")
    fig.text(0.5, 0.005, "UMAP is a visualization; quantitative preservation and mixing metrics are reported separately.",
             ha="center", fontsize=7, color="#444444")
    fig.savefig(OUT / "gse210248_raw_vs_harmony_umap_draft.png", dpi=120, bbox_inches="tight")
    fig.savefig(OUT / "gse210248_raw_vs_harmony_umap_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "gse210248_raw_vs_harmony_umap_draft.svg", bbox_inches="tight")
    plt.close(fig)

    params = {
        "project": PROJECT, "cells": int(x.shape[0]), "genes": int(x.shape[1]), "hvg": 2000,
        "normalization": "CP10K then log1p", "dimension_reduction": "TruncatedSVD 30 components",
        "batch_method": "Harmony", "batch_variable": "new.ident (six biological samples)",
        "protected_interpretation": "Disease was not supplied as a batch variable.",
        "umap": reducer, "metric_silhouette_sample_size": 5000, "random_seed": SEED,
        "limitations": "Sensitivity analysis, not exact Seurat integration/SNN reproduction; original integrated assay and full author parameters were unavailable."
    }
    (OUT / "harmony_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    print(metrics_df.to_string(index=False))
    print(json.dumps(params, indent=2))


if __name__ == "__main__":
    main()

