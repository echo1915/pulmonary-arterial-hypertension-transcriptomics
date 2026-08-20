"""Assemble review-v1 Figure 3: cell atlas, marker validation, and PDE8B SMC disease contrasts."""

from __future__ import annotations

from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


O = PROJECT_DATA_ROOT / "outputs"
UMAP = PROJECT_DATA_ROOT / "processed" / "gse210248_harmony" / "gse210248_raw_vs_harmony_umap.csv.gz"
EXPRESSION = PROJECT_DATA_ROOT / "processed" / "gse210248_umap" / "gse210248_recomputed_umap.csv.gz"
MARKERS = O / "single-cell-panels-review-v1" / "Supplementary_overall_marker_dotplot_source_data.csv"
PDE8B = O / "pde8b-smc-disease-review-v1" / "Figure_PDE8B_SMC1_SMC2_disease_source_data.csv"
STATS = O / "scrna-patient-pseudobulk-audit" / "pde8b_patient_level_effect_summary.csv"
OUT = O / "figure3-localization-review-v1"

ORDER = ["Fibro", "Endo 1", "Endo 2", "Mono / Macs", "T cells", "NK cells",
         "B cells", "DC", "Granulocytes", "Epithelial", "Mast cells", "SMC 1", "SMC 2"]
MARKER_BLOCKS = {
    "Fibro": ["COL1A1", "DCN"], "Endo 1": ["PECAM1", "VWF"],
    "Endo 2": ["CA4", "RGCC"], "Mono / Macs": ["LYZ", "C1QC"],
    "T cells": ["CD3D", "CD3E"], "NK cells": ["NKG7", "GNLY"],
    "B cells": ["CD79A", "MS4A1"], "DC": ["FCER1A", "CD1C"],
    "Granulocytes": ["S100A8", "S100A9"], "Epithelial": ["EPCAM", "KRT19"],
    "Mast cells": ["TPSAB1", "KIT"], "SMC 1": ["TAGLN", "MYH11"],
    "SMC 2": ["COL3A1", "FN1"],
}
PALETTE = {
    "Fibro":"#56B4E9", "Endo 1":"#009E73", "Endo 2":"#CCB974",
    "Mono / Macs":"#0072B2", "T cells":"#E69F00", "NK cells":"#CC79A7",
    "B cells":"#F0E442", "DC":"#AA7BC3", "Granulocytes":"#6C6C9C",
    "Epithelial":"#E76F51", "Mast cells":"#8C6D5A",
    "SMC 1":"#D55E00", "SMC 2":"#0072B2",
}
DONOR, PAH, INK, GRID = "#0072B2", "#D55E00", "#263746", "#E2E6E9"

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":7, "axes.spines.top":False, "axes.spines.right":False,
    "axes.linewidth":.75, "legend.frameon":False, "svg.fonttype":"none", "pdf.fonttype":42,
})


def panel_label(ax, label: str) -> None:
    fig = ax.figure
    y_fig = fig.transFigure.inverted().transform(ax.transAxes.transform((0, 1.10)))[1]
    x_fig = .0015 if ax.get_position().x0 < .50 else .55
    fig.text(x_fig, y_fig, label.upper(), fontsize=10, fontweight=700,
             va="bottom", ha="left")


def fixed_jitter(n: int) -> np.ndarray:
    return np.linspace(-.055, .055, n) if n > 1 else np.zeros(1)


def save(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression":"tiff_lzw"})


def plot_umap(ax, um: pd.DataFrame) -> None:
    for ct in ORDER:
        q = um[um.cell_type.eq(ct)]
        ax.scatter(q.harmony_UMAP1, q.harmony_UMAP2, s=.75, color=PALETTE[ct],
                   alpha=.65, linewidths=0, rasterized=True, label=ct)
    ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .50), ncol=1, fontsize=5.0,
              markerscale=3.0, handletextpad=.2, columnspacing=.55,
              borderaxespad=0, labelspacing=.18)


def plot_pde8b_umap(ax, um: pd.DataFrame, fig) -> None:
    zero = um.PDE8B_log_cp10k.le(0)
    ax.scatter(um.loc[zero, "harmony_UMAP1"], um.loc[zero, "harmony_UMAP2"],
               s=.65, color="#D9DDE0", alpha=.42, linewidths=0, rasterized=True)
    pos = um.loc[~zero].sort_values("PDE8B_log_cp10k")
    sc = ax.scatter(pos.harmony_UMAP1, pos.harmony_UMAP2, c=pos.PDE8B_log_cp10k,
                    s=4.0, cmap="cividis", vmin=0, vmax=pos.PDE8B_log_cp10k.quantile(.98),
                    linewidths=0, rasterized=True)
    ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_title("PDE8B expression", fontstyle="italic", pad=5)
    cb = fig.colorbar(sc, ax=ax, fraction=.040, pad=.025)
    cb.set_label("Expression (log CP10K)", fontsize=5.2)
    cb.ax.tick_params(labelsize=5)


def plot_markers(ax, dot: pd.DataFrame, fig) -> None:
    genes = [g for ct in ORDER for g in MARKER_BLOCKS[ct] if g in dot.gene.unique()]
    starts, cursor = {}, 0
    for ct in ORDER:
        block = [g for g in MARKER_BLOCKS[ct] if g in genes]
        starts[ct] = (cursor, cursor + len(block)); cursor += len(block)
    for k, ct in enumerate(ORDER):
        lo, hi = starts[ct]
        if k % 2: ax.axvspan(lo-.5, hi-.5, color="#EEF0F2", zorder=0)
    for yi, ct in enumerate(ORDER):
        z = dot[dot.cell_type.eq(ct)].set_index("gene").reindex(genes)
        ax.scatter(np.arange(len(genes)), np.full(len(genes), yi),
                   s=5 + 100*z.fraction_expressed, c=z.scaled_expression,
                   cmap="cividis", vmin=0, vmax=1, edgecolor="white", linewidth=.2, zorder=2)
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=5)
    ax.set_yticks(range(len(ORDER)), ORDER, fontsize=5.8); ax.invert_yaxis()
    for ct, (lo, hi) in starts.items():
        ax.text((lo+hi-1)/2, -.85, ct.replace("Mono / Macs", "Mono/Macs"),
                ha="center", va="bottom", fontsize=5.0, rotation=35, color=PALETTE[ct])
    handles = [ax.scatter([], [], s=5+100*x, color="#6F7E58", edgecolor="none") for x in [.1,.3,.5]]
    ax.legend(handles, ["10%", "30%", "50%"], title="Detected", loc="upper left",
              bbox_to_anchor=(1.0, 1.0), fontsize=5.0, title_fontsize=5.1)
    sm = plt.cm.ScalarMappable(cmap="cividis", norm=plt.Normalize(0,1))
    cb = fig.colorbar(sm, ax=ax, fraction=.018, pad=.075)
    cb.set_label("Relative expression", fontsize=5.2)
    cb.set_ticks([0,1]); cb.set_ticklabels(["Low","High"]); cb.ax.tick_params(labelsize=5)


def plot_disease(ax, values: pd.DataFrame, stats: pd.DataFrame, cell_type: str) -> None:
    z = values[values.cell_type.eq(cell_type)]
    for xpos, (group, color) in enumerate([("Donor", DONOR), ("PAH", PAH)]):
        q = z[z.disease.eq(group)].sort_values("sample")
        ax.scatter(xpos + fixed_jitter(len(q)), q.log_cp10k, s=38, color=color,
                   edgecolor="white", linewidth=.7, alpha=.9, zorder=3)
        ax.hlines(q.log_cp10k.mean(), xpos-.19, xpos+.19, color=INK, lw=1.45, zorder=4)
    stat = stats[stats.analysis.eq(f"GSE210248 {cell_type}")].iloc[0]
    ax.text(.04, .96, f"Hedges' g = {stat.hedges_g:.2f}\nExact P = {stat.exact_two_sided_permutation_p:.2f}",
            transform=ax.transAxes, va="top", fontsize=6.6)
    ax.set_xticks([0,1], ["Donor\n(n=3)", "PAH\n(n=3)"])
    ax.set_xlim(-.42,1.42); ax.set_ylim(-.012,.282); ax.grid(axis="y", color=GRID, lw=.6)
    ax.set_title(cell_type.replace(" ", ""), color=PALETTE[cell_type], fontweight="semibold", pad=5)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    um = pd.read_csv(UMAP)
    expression = pd.read_csv(EXPRESSION, usecols=["barcode", "PDE8B_log_cp10k"])
    um = um.merge(expression, on="barcode", how="left", validate="one_to_one")
    if um.PDE8B_log_cp10k.isna().any():
        raise ValueError("PDE8B expression is missing after barcode merge")
    dot = pd.read_csv(MARKERS)
    values = pd.read_csv(PDE8B); stats = pd.read_csv(STATS)
    # The two embeddings are descriptive context, the marker panel validates annotation,
    # and the bottom row contains the patient-level inferential evidence.
    fig = plt.figure(figsize=(7.35, 9.05))
    gs = fig.add_gridspec(3, 1, left=.08, right=.94, top=.97, bottom=.065,
                          height_ratios=[1.0, 1.18, .82], hspace=.53)
    # Reserve a narrow middle column for the cell-type legend so both UMAP
    # plotting areas remain equal in physical size after adding the colorbar.
    top = gs[0].subgridspec(1, 3, width_ratios=[1.0, .27, 1.0], wspace=.10)
    a=fig.add_subplot(top[0,0])
    a_pos=a.get_position(); panel_shift=13/(72*fig.get_figwidth())
    a.set_position([a_pos.x0-panel_shift,a_pos.y0,a_pos.width,a_pos.height])
    plot_umap(a,um); panel_label(a,"a")
    b=fig.add_subplot(top[0,2]); plot_pde8b_umap(b,um,fig); panel_label(b,"b")
    c=fig.add_subplot(gs[1]); plot_markers(c,dot,fig); panel_label(c,"c")
    bottom=gs[2].subgridspec(1,2,wspace=.20)
    d=fig.add_subplot(bottom[0,0]); plot_disease(d,values,stats,"SMC 1"); panel_label(d,"d")
    e=fig.add_subplot(bottom[0,1],sharey=d)
    e_pos=e.get_position()
    e.set_position([e_pos.x0+panel_shift,e_pos.y0,e_pos.width,e_pos.height])
    plot_disease(e,values,stats,"SMC 2"); panel_label(e,"e")
    d.set_ylabel("PDE8B expression (log CP10K)")
    e.tick_params(labelleft=False)
    save(fig, OUT / "Figure_3_cell_atlas_markers_PDE8B_SMC_review_v3")
    plt.close(fig)


if __name__ == "__main__":
    main()
