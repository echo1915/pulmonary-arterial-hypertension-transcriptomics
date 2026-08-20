"""Build Figure 6: PDE8B-associated cyclic-nucleotide program.

The figure is deliberately associative. Cross-cohort correlations are based on
subject-level bulk-lung expression and do not establish that PDE8B directly
controls cGMP abundance or pathway activity.
"""

from __future__ import annotations

from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import PROJECT_DATA_ROOT


OUT = PROJECT_DATA_ROOT / "outputs" / "figure6-cyclic-nucleotide-program"
COEX = PROJECT_DATA_ROOT / "outputs" / "pde8b-coexpression" / "pde8b_cross_cohort_coexpression_meta.csv"
PATHWAYS = PROJECT_DATA_ROOT / "outputs" / "figure5-program-review-v2" / "Figure_5_panel_b_selected_pathways_source_data.csv"

INK = "#243746"
MUTED = "#687782"
GRID = "#D9E0E4"
NEUTRAL = "#AAB4BA"
CAMP = "#2878A6"
BRIDGE = "#2A8C82"
CGMP = "#D97932"
PALE_BLUE = "#DCECF4"
PALE_TEAL = "#DCEEEB"
PALE_ORANGE = "#F8E7D8"

GENES = [
    "ADCY4", "PRKAR2A", "RAPGEF3", "PRKACA", "CREB1",
    "PDE3B", "PDE5A", "PRKG1", "ABCC4", "PDE3A",
]
BRANCH = {
    "ADCY4": "cAMP-associated", "PRKAR2A": "cAMP-associated",
    "RAPGEF3": "cAMP-associated", "PRKACA": "cAMP-associated",
    "CREB1": "cAMP-associated", "ABCC4": "cAMP-associated",
    "PDE3B": "cAMP–cGMP bridge", "PDE3A": "cAMP–cGMP bridge",
    "PDE5A": "cGMP–PKG", "PRKG1": "cGMP–PKG",
}
COLORS = {"cAMP-associated": CAMP, "cAMP–cGMP bridge": BRIDGE, "cGMP–PKG": CGMP}


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.2,
    "axes.titlesize": 8.4,
    "axes.titleweight": "semibold",
    "axes.labelsize": 7.4,
    "axes.edgecolor": INK,
    "axes.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def panel_label(ax: plt.Axes, label: str, x_fig: float | None = None) -> None:
    fig = ax.figure
    y_fig = fig.transFigure.inverted().transform(ax.transAxes.transform((0, 1.10)))[1]
    if x_fig is None:
        x_fig = .015 if ax.get_position().x0 < .50 else .55
    fig.text(x_fig, y_fig, label.upper(), fontsize=9, fontweight=700,
             color=INK, va="bottom", ha="left")


def adjust_axes(ax: plt.Axes, *, dx_points: float = 0, dy_points: float = 0,
                width_scale: float = 1, height_scale: float = 1) -> None:
    """Reproduce the author's panel translation and proportional resizing."""
    pos = ax.get_position()
    fig = ax.figure
    ax.set_position([
        pos.x0 + dx_points / (72 * fig.get_figwidth()),
        pos.y0 + dy_points / (72 * fig.get_figheight()),
        pos.width * width_scale,
        pos.height * height_scale,
    ])


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    co = pd.read_csv(COEX)
    missing = sorted(set(GENES) - set(co["symbol"]))
    if missing:
        raise ValueError(f"Missing prespecified genes: {missing}")
    selected = co[co.symbol.isin(GENES)].copy()
    selected["branch"] = selected.symbol.map(BRANCH)
    selected["significant"] = selected.adjusted_fdr.lt(0.05)
    selected["display_order"] = selected.symbol.map({g: i for i, g in enumerate(GENES)})
    selected = selected.sort_values("display_order")
    pathways = pd.read_csv(PATHWAYS)
    pathways = pathways[(pathways.query_set == "core_positive")].copy()
    return selected, pathways


def plot_evidence_map(ax: plt.Axes, data: pd.DataFrame) -> None:
    ax.axis("off")
    positions = {
        "ADCY4": (.10, .72), "PDE8B": (.10, .38), "cAMP": (.39, .56),
        "PRKAR2A": (.40, .84), "RAPGEF3": (.60, .78), "PDE3B": (.59, .40),
        "cGMP": (.80, .56), "PDE5A": (.79, .24), "PRKG1": (.93, .78),
    }
    node_style = {
        "ADCY4": (PALE_BLUE, CAMP), "PDE8B": (CAMP, CAMP),
        "cAMP": (PALE_BLUE, CAMP), "PRKAR2A": (PALE_BLUE, CAMP),
        "RAPGEF3": (PALE_BLUE, CAMP), "PDE3B": (PALE_TEAL, BRIDGE),
        "cGMP": (PALE_ORANGE, CGMP), "PDE5A": (PALE_ORANGE, CGMP),
        "PRKG1": (PALE_ORANGE, CGMP),
    }
    node_radius = {
        "PDE8B": .058, "cAMP": .058, "cGMP": .058,
        "ADCY4": .047, "PRKAR2A": .053, "RAPGEF3": .053,
        "PDE3B": .047, "PDE5A": .047, "PRKG1": .047,
    }

    def edge(start: str, end: str, color: str, *, inhibit: bool = False,
             curved: float = 0.0, width: float = 1.5) -> None:
        style = "-[" if inhibit else "-|>"
        ax.annotate("", xy=positions[end], xytext=positions[start], xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle=style, color=color, lw=width,
                                    mutation_scale=9, shrinkA=16, shrinkB=16,
                                    connectionstyle=f"arc3,rad={curved}"))

    edge("ADCY4", "cAMP", CAMP)
    edge("PDE8B", "cAMP", CAMP, inhibit=True, curved=-.08, width=2.1)
    edge("cAMP", "PRKAR2A", CAMP)
    edge("cAMP", "RAPGEF3", CAMP, curved=-.08)
    edge("PDE3B", "cAMP", BRIDGE, inhibit=True, curved=-.08, width=1.8)
    edge("cGMP", "PDE3B", BRIDGE, inhibit=True, curved=-.08, width=1.8)
    edge("PDE5A", "cGMP", CGMP, inhibit=True, curved=.06, width=2.0)
    edge("cGMP", "PRKG1", CGMP, curved=-.06, width=2.0)

    for name, (x, y) in positions.items():
        face, edge_color = node_style[name]
        radius = node_radius[name]
        circle = mpl.patches.Circle((x, y), radius, transform=ax.transAxes,
                                    facecolor=face, edgecolor=edge_color, lw=1.25, zorder=3)
        ax.add_patch(circle)
        if name in {"PDE8B", "cAMP", "cGMP"}:
            node_fontsize = 6.6
        elif len(name) >= 7:
            node_fontsize = 5.0
        else:
            node_fontsize = 5.7
        ax.text(x, y, name, transform=ax.transAxes, ha="center", va="center",
                fontsize=node_fontsize,
                fontweight="bold" if name in {"PDE8B", "cAMP", "cGMP"} else "normal",
                color="white" if name == "PDE8B" else INK, zorder=4)

    ax.text(.08, .08, "cAMP", color=CAMP, transform=ax.transAxes, fontweight="bold", fontsize=6.5)
    ax.text(.24, .08, "PDE3 bridge", color=BRIDGE, transform=ax.transAxes, fontweight="bold", fontsize=6.5)
    ax.text(.50, .08, "cGMP–PKG", color=CGMP, transform=ax.transAxes, fontweight="bold", fontsize=6.5)


def plot_forest(ax: plt.Axes, data: pd.DataFrame) -> None:
    d = data.sort_values("adjusted_meta_r")
    y = np.arange(len(d))
    ax.axvline(0, color=MUTED, lw=.7, ls="--")
    for i, (_, row) in enumerate(d.iterrows()):
        color = COLORS[row.branch] if row.significant else NEUTRAL
        ax.plot([row.adjusted_ci_low_r, row.adjusted_ci_high_r], [i, i], color=color, lw=1.25)
        ax.scatter(row.adjusted_meta_r, i, s=29, color=color, edgecolor="white", lw=.5, zorder=3)
    ax.set_yticks(y, d.symbol)
    ax.set_xlim(-.25, .72)
    ax.set_xlabel("Meta-correlation with PDE8B (95% CI)")
    ax.grid(axis="x", color=GRID, lw=.45)


def plot_heatmap(ax: plt.Axes, data: pd.DataFrame) -> None:
    cohorts = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]
    matrix = data.set_index("symbol")[[f"r_{x}" for x in cohorts]].loc[GENES]
    cmap = LinearSegmentedColormap.from_list("corr", ["#6D89A6", "#F6F7F7", "#D47A38"])
    im = ax.imshow(matrix, cmap=cmap, vmin=-.6, vmax=.6, aspect="auto")
    ax.set_xticks(range(4), cohorts, rotation=28, ha="right")
    ax.set_yticks(range(len(GENES)), GENES)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                    fontsize=5.6, color="white" if abs(value) > .38 else INK)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=.032, pad=.018)
    cbar.ax.set_title("r", fontsize=6.5, pad=3)
    cbar.ax.tick_params(labelsize=5.8)


def plot_pathways(ax: plt.Axes, pathways: pd.DataFrame) -> None:
    p = pathways.sort_values("fdr", ascending=False).copy()
    p["display_name"] = p.term_name.replace({
        "cilium organization": "cilium organization",
        "muscle system process": "muscle process",
        "actin filament-based process": "actin-filament process",
        "cGMP-PKG signaling pathway": "cGMP–PKG",
        "smooth muscle contraction": "SMC contraction",
        "cellular response to calcium ion": "calcium-ion response",
    })
    y = np.arange(len(p))
    size = 28 + 150 * p.recall / p.recall.max()
    color = [CGMP if "cGMP" in x else BRIDGE for x in p.term_name]
    ax.hlines(y, 0, -np.log10(p.fdr), color=GRID, lw=1.0)
    ax.scatter(-np.log10(p.fdr), y, s=size, c=color, edgecolor="white", lw=.6)
    ax.set_yticks(y, p.display_name)
    ax.set_xlabel("−log10(FDR)")
    ax.grid(axis="x", color=GRID, lw=.45)


def save_bundle(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})


def main() -> None:
    data, pathways = load_data()
    fig = plt.figure(figsize=(7.25, 6.65))
    gs = fig.add_gridspec(2, 2, left=.085, right=.985, bottom=.095, top=.965,
                          width_ratios=[1.08, .92], height_ratios=[.92, 1.08],
                          wspace=.64, hspace=.44)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    adjust_axes(axes[0], dx_points=-11.5, dy_points=-31.15,
                width_scale=1.094, height_scale=1.106)
    adjust_axes(axes[1], dx_points=-10.4, dy_points=-17.1)
    adjust_axes(axes[2], dx_points=12.7, dy_points=-.1)
    adjust_axes(axes[3], dx_points=1.6, dy_points=-.1)
    plot_evidence_map(axes[0], data)
    plot_forest(axes[1], data)
    plot_heatmap(axes[2], data)
    plot_pathways(axes[3], pathways)
    for ax, label in zip(axes, "abcd"):
        panel_label(ax, label, x_fig=.01 if label in "ac" else .52)
    stem = OUT / "Figure_6_PDE8B_cyclic_nucleotide_program"
    save_bundle(fig, stem)
    plt.close(fig)

    OUT.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUT / "Figure_6_panels_a-c_source_data.csv", index=False)
    pathways.to_csv(OUT / "Figure_6_panel_d_source_data.csv", index=False)
    legend = """# Fig. 6 | PDE8B marks a cross-cohort cyclic-nucleotide program

**a,** Data-grounded organization of the PDE8B-associated cyclic-nucleotide program. Values are random-effects meta-correlations across four independent bulk-lung cohorts. PDE8B is biochemically cAMP-specific; links to the cGMP–PKG branch represent transcriptional coupling rather than direct cGMP hydrolysis. **b,** Random-effects meta-correlations and 95% confidence intervals for prespecified cAMP-associated, PDE3 bridge and cGMP–PKG genes. Coloured estimates pass Benjamini–Hochberg FDR < 0.05; grey estimates do not. Cohort/subject—not gene—is the biological replication level within each contributing analysis. **c,** Cohort-specific correlations for the same prespecified genes. **d,** Functional enrichment of the positive high-confidence PDE8B core module; point size reflects pathway recall and labels show the number of intersecting genes. Enrichment FDR values were obtained using the existing pathway-analysis pipeline. The figure supports an associative PDE8B-linked cyclic-nucleotide program and does not establish direct regulation of cGMP or causal pathway activity. Source data are provided with the figure.
"""
    (OUT / "Figure_6_legend.md").write_text(legend, encoding="utf-8")
    qa = """# Figure 6 QA notes

- Core claim: PDE8B covaries with a cross-cohort cyclic-nucleotide program spanning cAMP-associated genes, PDE3B and the PDE5A–PRKG1/cGMP–PKG branch.
- Evidence type: subject-level bulk-lung cross-cohort correlation and gene-set enrichment; associative, not causal.
- Prespecified genes are all retained, including non-significant PRKACA, CREB1, PDE3A and ABCC4 results.
- No cell-level inferential statistics are used.
- Panel d reuses the audited positive-core enrichment table from Figure 5; no pathway was recomputed or selectively added.
- Exclusions: none among the ten prespecified cyclic-nucleotide genes; all four primary cohorts are shown.
- Export: editable SVG/PDF, 300-dpi PNG and LZW-compressed 600-dpi TIFF.
"""
    (OUT / "Figure_6_QA_notes.md").write_text(qa, encoding="utf-8")
    print(stem)


if __name__ == "__main__":
    main()
