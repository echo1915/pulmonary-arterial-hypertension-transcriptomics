from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


RAW = PROJECT_DATA_ROOT / "raw" / "current-data-snapshot" / "pah_scrna"
BASE = RAW / "gse210248" / "Human_PA_Integrated"
OUT = PROJECT_DATA_ROOT / "outputs" / "smc-gene-state-trends-review-v1"
CONTRACTILE = ["ACTA2", "TAGLN", "MYH11", "CNN1", "MYL9", "TPM2", "CALD1"]
SYNTHETIC = ["KLF4", "COL1A1", "COL3A1", "FN1", "VIM"]
TARGETS = ["MYH11", "FN1", "PDE8B"]
COLORS = {"SMC 1": "#E86B6B", "SMC 2": "#42B883"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def moving_average(x: np.ndarray, y: np.ndarray, bins: int = 28) -> tuple[np.ndarray, np.ndarray]:
    edges = np.quantile(x, np.linspace(0, 1, bins + 1))
    edges = np.unique(edges)
    centers, means = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (x >= lo) & (x <= hi if hi == edges[-1] else x < hi)
        if mask.sum() >= 10:
            centers.append(np.median(x[mask]))
            means.append(np.mean(y[mask]))
    return np.asarray(centers), np.asarray(means)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(RAW / "GSE210248_Human_PA_Metadata.csv.gz", sep=";", decimal=",", index_col=0)
    with gzip.open(BASE / "features.tsv.gz", "rt") as handle:
        features = [line.rstrip().split("\t")[1] for line in handle]
    with gzip.open(BASE / "barcodes.tsv.gz", "rt") as handle:
        barcodes = [line.strip() for line in handle]
    if barcodes != list(meta.index):
        raise ValueError("Matrix columns and metadata barcodes are not aligned")

    wanted_genes = set(CONTRACTILE + SYNTHETIC + TARGETS)
    wanted_rows = {i + 1: gene for i, gene in enumerate(features) if gene in wanted_genes}
    smc_mask = meta["Cell_annotation"].isin(["SMC 1", "SMC 2"]).to_numpy()
    smc_cols = np.flatnonzero(smc_mask) + 1
    col_to_local = {int(col): i for i, col in enumerate(smc_cols)}
    counts = {gene: np.zeros(len(smc_cols), dtype=float) for gene in wanted_genes}
    libraries = np.zeros(len(smc_cols), dtype=float)
    with gzip.open(BASE / "matrix.mtx.gz", "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                break
        for line in handle:
            row, col, value = line.split()
            row, col, value = int(row), int(col), float(value)
            local = col_to_local.get(col)
            if local is None:
                continue
            libraries[local] += value
            gene = wanted_rows.get(row)
            if gene is not None:
                counts[gene][local] += value

    cell_meta = meta.loc[smc_mask, ["new.ident", "disease", "Cell_annotation"]].copy()
    cell_meta.columns = ["sample", "disease", "cell_type"]
    expr = pd.DataFrame(index=cell_meta.index)
    for gene in sorted(wanted_genes):
        expr[gene] = np.log1p(10000 * counts[gene] / np.maximum(libraries, 1))
    z = expr.apply(lambda v: (v - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else 0)
    state_score = z[SYNTHETIC].mean(axis=1) - z[CONTRACTILE].mean(axis=1)
    source = cell_meta.assign(state_score=state_score.to_numpy())
    for gene in TARGETS:
        source[gene] = expr[gene].to_numpy()
    source.to_csv(OUT / "Figure_SMC_gene_state_trends_source_data.csv", index_label="barcode")

    fig, axes = plt.subplots(1, 3, figsize=(8.4, 3.15), sharex=True)
    for label, ax, gene in zip("abc", axes, TARGETS):
        order = np.argsort(source["state_score"].to_numpy())
        q = source.iloc[order]
        for cell_type in ["SMC 1", "SMC 2"]:
            part = q[q["cell_type"].eq(cell_type)]
            ax.scatter(part["state_score"], part[gene], s=5.0, alpha=0.24,
                       color=COLORS[cell_type], linewidths=0, rasterized=True,
                       label=cell_type)
        x = q["state_score"].to_numpy()
        y = q[gene].to_numpy()
        bx, by = moving_average(x, y)
        ax.plot(bx, by, color="#D64A61", lw=2.0, zorder=5)
        ax.set_title(gene, fontstyle="italic", pad=8)
        ax.set_xlabel("SMC state score\nSynthetic/ECM  −  Contractile")
        ax.set_ylabel("Expression (log CP10K)")
        ax.grid(False)
        ax.text(-0.14, 1.04, label, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")
    handles, labels = axes[-1].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, markerscale=2.2)
    fig.subplots_adjust(left=0.08, right=0.90, top=0.88, bottom=0.23, wspace=0.32)

    stem = OUT / "Figure_SMC_MYH11_FN1_PDE8B_state_trends_review_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
