from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


BLUE = "#4169A1"
PURPLE = "#7652A6"
MAGENTA = "#C51B8A"
RED = "#F21A24"
GREY = "#C7C9C9"


def save_all(fig: plt.Figure, prefix: Path) -> None:
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(prefix.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main() -> None:
    outputs = PROJECT_DATA_ROOT / "outputs"
    continuum_dir = outputs / "smc-state-continuum-v1"
    enrichment_dir = outputs / "smc-symmetric-enrichment-v1"
    out_dir = outputs / "smc-directional-bars-review-v1"
    out_dir.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(continuum_dir / "Figure_SMC1_SMC2_state_continuum_source_data.csv")
    enrich = pd.read_csv(enrichment_dir / "Figure_SMC1_SMC2_symmetric_enrichment_source_data.csv")

    score_cols = ["Contractile", "Mural remodeling", "PDE8B", "Synthetic/ECM"]
    # Restrict the directional summary to PAH patients. Scores were standardized
    # across all samples, so averaging donor and PAH deltas together collapses
    # near zero and obscures the within-PAH SMC2-versus-SMC1 contrast.
    pah_scores = scores[scores["disease"].eq("PAH")].copy()
    wide = pah_scores.pivot(index="sample", columns="cell_type", values=score_cols)
    deltas = pd.DataFrame(
        {
            "program": score_cols,
            "mean_delta_smc2_minus_smc1": [
                (wide[(col, "SMC 2")] - wide[(col, "SMC 1")]).mean() for col in score_cols
            ],
        }
    )
    # Matplotlib barh displays the first category at the bottom. PDE8B is
    # omitted here because its paired SMC2-minus-SMC1 difference is near zero.
    panel_a_order = ["Synthetic/ECM", "Contractile", "Mural remodeling"]
    deltas = deltas[deltas["program"].isin(panel_a_order)].copy()
    deltas["program"] = pd.Categorical(deltas["program"], categories=panel_a_order, ordered=True)
    deltas = deltas.sort_values("program")

    enrich = enrich.copy()
    enrich["minus_log10_fdr"] = -np.log10(enrich["fdr"].astype(float).clip(lower=1e-300))
    selected_terms = [
        "Extracellular matrix organization",
        "collagen fibril organization",
        "ECM proteoglycans",
        "vasculature development",
        "actin cytoskeleton organization",
        "oxidative phosphorylation",
        "Smooth Muscle Contraction",
        "muscle contraction",
    ]
    pathways = enrich[enrich["term"].isin(selected_terms)].copy()
    pathways = pathways.sort_values("minus_log10_fdr", ascending=True)

    fig = plt.figure(figsize=(6.1, 7.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.82, 1.18], hspace=0.50)

    ax1 = fig.add_subplot(gs[0])
    y = np.arange(len(deltas))
    vals = deltas["mean_delta_smc2_minus_smc1"].to_numpy()
    max_abs = max(abs(vals).max(), 0.1)
    norm = mpl.colors.TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("direction", [BLUE, PURPLE, MAGENTA, RED])
    ax1.barh(y, vals, color=cmap(norm(vals)), height=0.68)
    ax1.axvline(0, color="#263746", lw=0.8)
    ax1.set_yticks(y, deltas["program"])
    ax1.set_xlabel("Mean paired score difference in PAH")
    ax1.set_xlim(-max_abs * 1.35, max_abs * 1.35)
    ax1.grid(axis="x", color="#E5E7E9", lw=0.6, zorder=0)
    ax1.tick_params(axis="y", length=0)
    ax1.text(0.22, -0.25, "\u2190 SMC1-associated", transform=ax1.transAxes, color=BLUE,
             ha="center", va="top", fontsize=8)
    ax1.text(0.78, -0.25, "SMC2-associated \u2192", transform=ax1.transAxes, color=RED,
             ha="center", va="top", fontsize=8)
    ax1.text(-0.12, 1.04, "a", transform=ax1.transAxes, fontsize=12, fontweight="bold")

    ax2 = fig.add_subplot(gs[1])
    y2 = np.arange(len(pathways))
    score = pathways["minus_log10_fdr"].to_numpy()
    side_colors = np.where(pathways["side"].eq("SMC2"), PURPLE, GREY)
    ax2.barh(y2, score, color=side_colors, height=0.70)
    ax2.set_yticks(y2, pathways["term"])
    ax2.set_xlabel(r"$-\log_{10}$(FDR)")
    ax2.grid(axis="x", color="#E5E7E9", lw=0.6, zorder=0)
    ax2.tick_params(axis="y", length=0)
    ax2.text(-0.12, 1.04, "b", transform=ax2.transAxes, fontsize=12, fontweight="bold")
    handles = [
        mpl.patches.Patch(color=GREY, label="SMC1-associated"),
        mpl.patches.Patch(color=PURPLE, label="SMC2-associated"),
    ]
    ax2.legend(handles=handles, loc="lower right", frameon=False, fontsize=7)

    fig.subplots_adjust(left=0.43, right=0.96, top=0.97, bottom=0.09)
    prefix = out_dir / "Figure_SMC_directional_bars_review_v1"
    save_all(fig, prefix)

    deltas.to_csv(out_dir / "Figure_SMC_directional_program_source_data.csv", index=False)
    pathways.to_csv(out_dir / "Figure_SMC_directional_pathway_source_data.csv", index=False)
    print(prefix.with_suffix(".png"))


if __name__ == "__main__":
    main()
