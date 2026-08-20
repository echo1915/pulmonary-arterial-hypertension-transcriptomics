"""Phase 1 reanalysis of GSE293580: QC, Harmony, clustering, and annotation evidence."""

from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import anndata as ad
import harmonypy as hm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import vstack


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import DATA_ROOT, INTERIM_ROOT, PROJECT_DATA_ROOT  # noqa: E402


RAW = DATA_ROOT / "pah_scrna" / "gse293580"
INTERIM = INTERIM_ROOT / "gse293580_scanpy"
PROCESSED = PROJECT_DATA_ROOT / "processed" / "gse293580_scanpy"
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis"
for folder in [INTERIM, PROCESSED, OUT]:
    folder.mkdir(parents=True, exist_ok=True)

SEED = 20260810
GROUPS = {
    "GSM8885793": ("BA-044", "Donor", "pool1"), "GSM8885794": ("BA-048", "Donor", "pool1"),
    "GSM8885795": ("BA-053", "Donor", "pool1"), "GSM8885796": ("BA-058", "Donor", "pool1"),
    "GSM8885797": ("BA-024", "SSc-PAH", "pool2"), "GSM8885798": ("BA-025", "SSc-PAH", "pool2"),
    "GSM8885799": ("ST-005", "SSc-PAH", "pool2"), "GSM8885800": ("ST-059", "SSc-PAH", "pool2"),
    "GSM8885801": ("BA-008", "IPAH", "pool3"), "GSM8885802": ("BA-014", "IPAH", "pool3"),
    "GSM8885803": ("CC-030", "IPAH", "pool3"), "GSM8885804": ("ST-042", "IPAH", "pool3"),
}

MARKERS = {
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG", "CLDN5", "RAMP2"],
    "SMC": ["ACTA2", "TAGLN", "MYH11", "CNN1", "MYL9", "TPM2"],
    "Pericyte": ["RGS5", "CSPG4", "MCAM", "PDGFRB", "ABCC9", "KCNJ8"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "AT1": ["AGER", "PDPN", "CAV1", "CAV2", "EMP2"],
    "AT2": ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "ABCA3"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS", "AIF1"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "LTB"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "B_cell": ["CD79A", "MS4A1", "CD37", "CD74", "CD22"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"],
    "Platelet": ["PPBP", "PF4", "GP9", "ITGA2B", "NRGN"],
    "Proliferative": ["MKI67", "TOP2A", "UBE2C", "CENPF", "STMN1"],
}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def read_lines(path: Path, field: int = 0) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t")[field] for line in handle]


def load_all() -> ad.AnnData:
    matrices, obs_rows, reference_genes = [], [], None
    for matrix_path in sorted(RAW.glob("GSM*_matrix.mtx.gz")):
        accession = matrix_path.name.split("_")[0]
        sample, condition, pool = GROUPS[accession]
        prefix = matrix_path.name.replace("matrix.mtx.gz", "")
        features_path = RAW / f"{prefix}features.tsv.gz"
        barcodes_path = RAW / f"{prefix}barcodes.tsv.gz"
        genes = read_lines(features_path, field=1)
        barcodes = read_lines(barcodes_path)
        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise ValueError(f"Feature order differs for {accession}")
        x = mmread(str(matrix_path)).tocsr().T.astype(np.float32)
        if x.shape != (len(barcodes), len(genes)):
            raise ValueError(f"Shape mismatch for {accession}: {x.shape}")
        matrices.append(x)
        obs_rows.extend({
            "cell_id": f"{sample}__{barcode}", "barcode": barcode, "accession": accession,
            "sample": sample, "condition": condition, "pool": pool,
        } for barcode in barcodes)
    obs = pd.DataFrame(obs_rows).set_index("cell_id")
    var = pd.DataFrame(index=pd.Index(reference_genes, name="gene"))
    var_names = pd.Index(reference_genes).astype(str)
    var.index = var_names
    adata = ad.AnnData(X=vstack(matrices, format="csr"), obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def cluster_marker_scores(adata: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    mean_rows = []
    for cluster in sorted(adata.obs[cluster_key].astype(str).unique(), key=lambda x: int(x)):
        mask = adata.obs[cluster_key].astype(str).eq(cluster).to_numpy()
        row = {"cluster": cluster, "cells": int(mask.sum())}
        for lineage, genes in MARKERS.items():
            present = [g for g in genes if g in adata.var_names]
            block = adata[mask, present].X
            row[lineage] = float(np.asarray(block.mean(axis=0)).mean()) if present else np.nan
            row[f"{lineage}_genes"] = len(present)
        mean_rows.append(row)
    scores = pd.DataFrame(mean_rows)
    lineage_cols = list(MARKERS)
    z = scores[lineage_cols].apply(lambda v: (v - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else 0)
    scores["top_lineage"] = z.idxmax(axis=1)
    scores["top_lineage_z"] = z.max(axis=1)
    for col in lineage_cols:
        scores[f"z_{col}"] = z[col]
    return scores


def main() -> None:
    np.random.seed(SEED)
    adata = load_all()
    if adata.n_obs != 13650:
        raise ValueError(f"Expected 13,650 GEO-filtered cells, found {adata.n_obs}")
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sample_qc = (adata.obs.groupby(["sample", "condition", "pool"], observed=True)
                 .agg(cells=("accession", "size"), median_genes=("n_genes_by_counts", "median"),
                      median_umis=("total_counts", "median"), median_pct_mt=("pct_counts_mt", "median"),
                      max_pct_mt=("pct_counts_mt", "max")).reset_index())
    sample_qc.to_csv(OUT / "gse293580_sample_qc.csv", index=False, encoding="utf-8-sig")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="sample")
    sc.pp.scale(adata, zero_center=False, max_value=10)
    sc.tl.pca(adata, n_comps=40, use_highly_variable=True, random_state=SEED)
    harmony = hm.run_harmony(adata.obsm["X_pca"], adata.obs[["sample"]], "sample", random_state=SEED, verbose=False)
    corrected = np.asarray(harmony.Z_corr)
    if corrected.shape[0] != adata.n_obs and corrected.shape[1] == adata.n_obs:
        corrected = corrected.T
    if corrected.shape[0] != adata.n_obs:
        raise ValueError(f"Unexpected Harmony shape {corrected.shape}")
    adata.obsm["X_pca_harmony"] = corrected
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(adata, min_dist=0.30, random_state=SEED)
    for resolution in [0.3, 0.5, 0.8]:
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}", random_state=SEED, flavor="igraph", n_iterations=2)
    primary = "leiden_0.5"
    sc.tl.rank_genes_groups(adata, groupby=primary, method="wilcoxon", use_raw=True, pts=True)
    marker_frame = sc.get.rank_genes_groups_df(adata, group=None)
    marker_frame.groupby("group", observed=True).head(25).to_csv(OUT / "gse293580_cluster_top25_markers.csv", index=False, encoding="utf-8-sig")
    scores = cluster_marker_scores(adata, primary)
    scores.to_csv(OUT / "gse293580_cluster_lineage_scores.csv", index=False, encoding="utf-8-sig")

    obs_cols = ["barcode", "accession", "sample", "condition", "pool", "n_genes_by_counts", "total_counts", "pct_counts_mt", "leiden_0.3", "leiden_0.5", "leiden_0.8"]
    coords = adata.obs[obs_cols].copy()
    coords["UMAP1"] = adata.obsm["X_umap"][:, 0]
    coords["UMAP2"] = adata.obsm["X_umap"][:, 1]
    coords.to_csv(PROCESSED / "gse293580_harmony_umap_clusters.csv.gz", compression="gzip")
    adata.write_h5ad(INTERIM / "gse293580_phase1_harmony_clustered.h5ad", compression="gzip")

    cluster_counts = adata.obs.groupby([primary, "sample", "condition"], observed=True).size().rename("cells").reset_index()
    cluster_counts.to_csv(OUT / "gse293580_cluster_sample_counts.csv", index=False, encoding="utf-8-sig")
    summary = {
        "cells": int(adata.n_obs), "genes": int(adata.n_vars), "samples": 12,
        "conditions": adata.obs.condition.value_counts().to_dict(),
        "clusters": {key: int(adata.obs[key].nunique()) for key in ["leiden_0.3", "leiden_0.5", "leiden_0.8"]},
        "primary_cluster_key": primary, "hvg": int(adata.var.highly_variable.sum()),
        "normalization": "per-cell total 10,000; log1p", "batch_correction": "Harmony 0.0.10 by patient/sample",
        "random_seed": SEED, "input_cells_note": "GEO matrices already contain the 13,650 cells reported after author QC.",
    }
    (OUT / "gse293580_phase1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    for ax, field, title in zip(axes, [primary, "condition", "sample"], ["Leiden clusters", "Disease group", "Patient/sample"]):
        labels = adata.obs[field].astype(str).to_numpy()
        categories = sorted(np.unique(labels))
        cmap = plt.get_cmap("tab20")
        for i, label in enumerate(categories):
            keep = labels == label
            ax.scatter(adata.obsm["X_umap"][keep, 0], adata.obsm["X_umap"][keep, 1], s=1.3, alpha=0.55, color=cmap(i % 20), rasterized=True)
        ax.set_title(title, loc="left", fontweight="bold"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("GSE293580 phase-1 reanalysis: Harmony-integrated 12-patient lung atlas", x=0.04, ha="left", fontsize=11.5, fontweight="bold")
    fig.text(0.04, 0.015, "UMAP is used for visualization. Final cell labels require marker-score and top-marker review; PDE8B was not used for clustering or annotation.", fontsize=7, color="#56616B")
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    fig.savefig(OUT / "gse293580_phase1_umap_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "gse293580_phase1_umap_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "gse293580_phase1_umap_draft.svg", bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    print(sample_qc.to_string(index=False))
    print(scores[["cluster", "cells", "top_lineage", "top_lineage_z"]].to_string(index=False))


if __name__ == "__main__":
    main()
