"""Phase 2 of GSE293580 reanalysis: unbiased mural-cell refinement.

PDE8B is deliberately excluded from clustering and annotation. This script
creates multi-resolution subclusters and marker/module evidence so that the
SMC boundary can be locked before patient-level hypothesis testing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import harmonypy as hm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


SEED = 20260810
INTERIM = PROJECT_DATA_ROOT / "interim" / "gse293580_scanpy"
PROCESSED = PROJECT_DATA_ROOT / "processed" / "gse293580_scanpy"
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis"
for folder in (INTERIM, PROCESSED, OUT):
    folder.mkdir(parents=True, exist_ok=True)

INPUT = INTERIM / "gse293580_phase1_harmony_clustered.h5ad"

MODULES = {
    "contractile_SMC": ["ACTA2", "TAGLN", "MYH11", "CNN1", "LMOD1", "SMTN", "DES", "MYL9", "TPM2"],
    "pericyte": ["RGS5", "CSPG4", "PDGFRB", "ABCC9", "KCNJ8", "COX4I2", "NOTCH3", "MCAM"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5"],
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def safe_score(adata, genes: list[str], name: str) -> None:
    present = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, present, score_name=name, random_state=SEED, use_raw=False)


def export_rank_table(adata, key: str, groupby: str, n: int = 30) -> pd.DataFrame:
    frame = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    frame["rank"] = frame.groupby("group").cumcount() + 1
    return frame.loc[frame["rank"] <= n, ["group", "rank", "names", "scores", "logfoldchanges", "pvals_adj"]]


def main() -> None:
    full = sc.read_h5ad(INPUT)
    # Phase-1 cluster 1 is the only cluster with concordant smooth-muscle and
    # pericyte marker enrichment. Selection is independent of PDE8B.
    mural = full[full.obs["leiden_0.5"].astype(str) == "1"].copy()
    # Phase 1 stores integer counts in the named layer and log-normalized
    # expression in .raw. Always rebuild subclustering from the count layer.
    mural.X = mural.layers["counts"].copy()
    mural.layers["counts"] = mural.X.copy()

    sc.pp.normalize_total(mural, target_sum=1e4)
    sc.pp.log1p(mural)
    mural.raw = mural
    sc.pp.highly_variable_genes(
        mural,
        n_top_genes=min(2000, mural.n_vars),
        batch_key="sample",
        flavor="seurat",
        subset=False,
    )
    sc.pp.scale(mural, max_value=10)
    sc.tl.pca(mural, n_comps=30, use_highly_variable=True, random_state=SEED)
    harmony = hm.run_harmony(mural.obsm["X_pca"], mural.obs, "sample", random_state=SEED, max_iter_harmony=20)
    mural.obsm["X_pca_harmony"] = harmony.Z_corr.T
    sc.pp.neighbors(mural, n_neighbors=15, n_pcs=25, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(mural, min_dist=0.30, random_state=SEED)

    resolutions = (0.2, 0.3, 0.4, 0.5, 0.6)
    resolution_summary = []
    for resolution in resolutions:
        key = f"mural_leiden_{resolution:.1f}"
        sc.tl.leiden(mural, resolution=resolution, key_added=key, random_state=SEED, flavor="igraph")
        counts = mural.obs[key].value_counts()
        resolution_summary.append({
            "resolution": resolution,
            "n_clusters": int(counts.size),
            "min_cluster_cells": int(counts.min()),
            "median_cluster_cells": float(counts.median()),
        })

    for module, genes in MODULES.items():
        safe_score(mural, genes, f"score_{module}")

    primary = "mural_leiden_0.4"
    sc.tl.rank_genes_groups(mural, groupby=primary, method="wilcoxon", pts=True, key_added="rank_mural")
    markers = export_rank_table(mural, "rank_mural", primary)
    markers.to_csv(OUT / "gse293580_mural_cluster_top30_markers.csv", index=False)

    score_cols = [f"score_{name}" for name in MODULES]
    profile = mural.obs.groupby(primary, observed=True)[score_cols].agg(["mean", "median"])
    profile.columns = [f"{a}_{b}" for a, b in profile.columns]
    sizes = mural.obs.groupby(primary, observed=True).size().rename("n_cells")
    profile = profile.join(sizes).reset_index()
    profile.to_csv(OUT / "gse293580_mural_cluster_module_profiles.csv", index=False)

    sample_counts = (
        mural.obs.groupby([primary, "condition", "sample"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    sample_counts.to_csv(OUT / "gse293580_mural_cluster_sample_counts.csv", index=False)
    pd.DataFrame(resolution_summary).to_csv(OUT / "gse293580_mural_resolution_summary.csv", index=False)

    coords = mural.obs[["sample", "condition", primary] + score_cols].copy()
    coords["UMAP1"] = mural.obsm["X_umap"][:, 0]
    coords["UMAP2"] = mural.obsm["X_umap"][:, 1]
    coords.to_csv(PROCESSED / "gse293580_mural_harmony_umap.csv.gz", compression="gzip")

    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.0))
    sc.pl.umap(mural, color=primary, ax=axes[0], show=False, title="Mural subclusters", frameon=False)
    sc.pl.umap(mural, color="condition", ax=axes[1], show=False, title="Condition", frameon=False)
    sc.pl.umap(mural, color="sample", ax=axes[2], show=False, title="Patient", frameon=False, legend_loc="none")
    fig.suptitle("GSE293580 mural-cell refinement (Harmony by patient)", fontsize=9, y=1.01)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"gse293580_phase2_mural_umap_draft.{suffix}", dpi=160 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)

    mural.write_h5ad(INTERIM / "gse293580_phase2_mural_clustered.h5ad", compression="gzip")
    summary = {
        "input": str(INPUT),
        "n_mural_cells": int(mural.n_obs),
        "n_patients": int(mural.obs["sample"].nunique()),
        "patient_coverage": sorted(mural.obs["sample"].unique().tolist()),
        "selection_rule": "phase-1 Leiden 0.5 cluster 1; selected by canonical mural markers; PDE8B excluded",
        "primary_descriptive_resolution": 0.4,
        "resolution_grid": resolution_summary,
    }
    (OUT / "gse293580_phase2_mural_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
