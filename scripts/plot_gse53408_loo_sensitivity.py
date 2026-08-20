"""Plot GSE53408 leave-one-sample-out candidate sensitivity results."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


IN_DIR = PROJECT_DATA_ROOT / "outputs" / "gse53408-sensitivity"
OUT = IN_DIR
COLORS = {"PDE8B": "#B84343", "PIEZO2": "#477AA8", "SLC16A12": "#187B7B"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def main() -> None:
    loo = pd.read_csv(IN_DIR / "gse53408_candidate_leave_one_sample_out.csv")
    summary = pd.read_csv(IN_DIR / "gse53408_candidate_loo_summary.csv")
    candidates = ["PDE8B", "PIEZO2", "SLC16A12"]

    fig = plt.figure(figsize=(10.8, 6.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1], height_ratios=[1, 1], hspace=0.46, wspace=0.30)

    ax = fig.add_subplot(gs[:, 0])
    y_base = {gene: 2 - i for i, gene in enumerate(candidates)}
    jitter = np.linspace(-0.18, 0.18, 23)
    for gene in candidates:
        z = loo[loo.symbol.eq(gene)].sort_values("omitted_gsm")
        y = y_base[gene] + jitter
        group_marker = z.omitted_group.map({"PAH": "o", "Control": "s"})
        for marker in ["o", "s"]:
            keep = group_marker.eq(marker).to_numpy()
            ax.scatter(z.loc[keep, "hedges_g"], y[keep], s=22, marker=marker,
                       facecolor=COLORS[gene], edgecolor="white", linewidth=0.35, alpha=0.78)
        full = summary.loc[summary.symbol.eq(gene), "full_hedges_g"].iloc[0]
        ax.vlines(full, y_base[gene] - 0.28, y_base[gene] + 0.28, color=COLORS[gene], lw=2.0)
    ax.axvline(0, color="#7A858F", lw=0.8)
    ax.set_yticks([2, 1, 0], candidates)
    ax.set_xlabel("Hedges g after omitting one sample")
    ax.set_title("A  Effect direction and magnitude remain stable", loc="left", fontweight="bold")
    ax.scatter([], [], marker="o", color="#68747E", label="Omitted PAH sample")
    ax.scatter([], [], marker="s", color="#68747E", label="Omitted control sample")
    ax.plot([], [], color="#68747E", lw=2, label="Full-cohort estimate")
    ax.legend(frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(3)
    ax.barh(y, summary["fdr_lt_0.05_fraction"] * 100,
            color=[COLORS[g] for g in summary.symbol], height=0.55)
    ax.set_yticks(y, summary.symbol)
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel("LOO iterations with FDR < 0.05 (%)")
    ax.set_title("B  Significance persists in all 23 iterations", loc="left", fontweight="bold")
    for i, value in enumerate(summary["fdr_lt_0.05_fraction"] * 100):
        ax.text(value - 2, i, f"{value:.0f}%", ha="right", va="center", color="white", fontweight="bold")

    ax = fig.add_subplot(gs[1, 1])
    z = loo.copy()
    z["abs_delta"] = z.delta_g_from_full.abs()
    top = z.sort_values("abs_delta", ascending=False).groupby("symbol", sort=False).head(1)
    y = np.arange(len(top))
    ax.hlines(y, 0, top.abs_delta, color="#CCD4DA", lw=1.3)
    ax.scatter(top.abs_delta, y, c=[COLORS[g] for g in top.symbol], s=42)
    labels = [f"{row.symbol}: {row.omitted_gsm}\n({row.omitted_group})" for _, row in top.iterrows()]
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Maximum |change in Hedges g|")
    ax.set_title("C  Most influential omitted sample", loc="left", fontweight="bold")
    for i, value in enumerate(top.abs_delta):
        ax.text(value + 0.008, i, f"{value:.3f}", va="center", fontsize=7)

    fig.suptitle("GSE53408 raw-lung validation is not driven by a single sample", x=0.04, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.04, 0.015, "LOO: leave one sample out. Each iteration refits genome-wide limma with robust empirical Bayes; FDR is recalculated across all genes.", fontsize=7, color="#56616B")
    fig.savefig(OUT / "Figure_GSE53408_leave_one_out_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_GSE53408_leave_one_out_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_GSE53408_leave_one_out_draft.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

