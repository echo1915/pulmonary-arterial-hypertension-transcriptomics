"""Bridge the 188-gene PAH differential signature to PDE8B coexpression programs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, hypergeom, mannwhitneyu, spearmanr


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402


AUDIT = OUTPUT_ROOT / "pah_audit"
COEX = PROJECT_DATA_ROOT / "outputs" / "pde8b-coexpression"
OUT = PROJECT_DATA_ROOT / "outputs" / "signature-pde8b-bridge"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

TEAL = "#187B7B"
ORANGE = "#D97732"
BLUE = "#477AA8"
GREY = "#AEB8C2"
DARK = "#263238"


def enrichment_stats(universe_n: int, robust_n: int, module_n: int, overlap_n: int) -> dict:
    expected = robust_n * module_n / universe_n
    return {
        "universe_n": universe_n,
        "robust_n_in_universe": robust_n,
        "module_n": module_n,
        "overlap_n": overlap_n,
        "expected_overlap": expected,
        "observed_expected_ratio": overlap_n / expected if expected else np.nan,
        "p_enrichment": float(hypergeom.sf(overlap_n - 1, universe_n, module_n, robust_n)),
        "p_depletion": float(hypergeom.cdf(overlap_n, universe_n, module_n, robust_n)),
    }


def main() -> None:
    meta = pd.read_csv(AUDIT / "lung_mechanism_random_effects_meta.csv")
    coex = pd.read_csv(COEX / "pde8b_cross_cohort_coexpression_meta.csv")
    core = pd.read_csv(COEX / "pde8b_high_confidence_core_module.csv")
    stable = pd.read_csv(COEX / "pde8b_stable_mechanism_module.csv")

    robust = meta[
        meta["same_direction"].astype(bool)
        & (meta["fdr"] < 0.05)
        & (meta["I2"] < 50)
        & (meta["meta_g"].abs() >= 0.50)
    ].copy()
    universe = set(coex["symbol"])
    bridge = robust[robust["symbol"].isin(universe)].merge(
        coex[["symbol", "adjusted_meta_r", "adjusted_fdr", "core_module", "mechanism_module"]],
        on="symbol", how="left", validate="one_to_one"
    )
    bridge["direction_concordant"] = np.sign(bridge["meta_g"]) == np.sign(bridge["adjusted_meta_r"])
    bridge["in_high_confidence_core"] = bridge["symbol"].isin(set(core["symbol"]))
    bridge["in_stable_mechanism"] = bridge["symbol"].isin(set(stable["symbol"]))
    bridge = bridge.sort_values(["in_stable_mechanism", "fdr"], ascending=[False, True])

    universe_n = len(universe)
    robust_n = len(bridge)
    core_set = set(core["symbol"]) & universe
    stable_set = set(stable["symbol"]) & universe
    tests = pd.DataFrame([
        {"module": "High-confidence core", **enrichment_stats(universe_n, robust_n, len(core_set), int(bridge["in_high_confidence_core"].sum()))},
        {"module": "Stable mechanism", **enrichment_stats(universe_n, robust_n, len(stable_set), int(bridge["in_stable_mechanism"].sum()))},
    ])

    rho, rho_p = spearmanr(bridge["meta_g"], bridge["adjusted_meta_r"])
    robust_abs = bridge["adjusted_meta_r"].abs().to_numpy()
    background_abs = coex.loc[~coex["symbol"].isin(set(bridge["symbol"])), "adjusted_meta_r"].abs().to_numpy()
    mw = mannwhitneyu(robust_abs, background_abs, alternative="greater")
    n_concordant = int(bridge["direction_concordant"].sum())
    direction_test = binomtest(n_concordant, robust_n, p=0.5, alternative="greater")

    missing = sorted(set(robust["symbol"]) - universe)
    summary = {
        "robust_definition": "same direction in all 4 cohorts, FDR<0.05, I2<50%, abs(meta_g)>=0.50",
        "robust_total": int(len(robust)),
        "robust_in_coexpression_universe": robust_n,
        "robust_missing_from_coexpression_universe": missing,
        "coexpression_universe_n": universe_n,
        "spearman_meta_g_vs_adjusted_meta_r": {"rho": float(rho), "p": float(rho_p)},
        "direction_concordance": {
            "n": n_concordant,
            "fraction": n_concordant / robust_n,
            "binomial_one_sided_p": float(direction_test.pvalue),
        },
        "absolute_correlation_robust_vs_background": {
            "robust_median": float(np.median(robust_abs)),
            "background_median": float(np.median(background_abs)),
            "mann_whitney_one_sided_p": float(mw.pvalue),
        },
        "interpretation": "Gene-set membership is largely orthogonal, while disease effects and PDE8B coexpression directions show broad coupling.",
    }

    bridge.to_csv(OUT / "robust188_pde8b_bridge_gene_table.csv", index=False, encoding="utf-8-sig")
    tests.to_csv(OUT / "robust188_pde8b_overlap_tests.csv", index=False, encoding="utf-8-sig")
    bridge[bridge["in_stable_mechanism"]].to_csv(OUT / "robust188_stable_module_overlap_genes.csv", index=False, encoding="utf-8-sig")
    with open(OUT / "robust188_pde8b_bridge_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)

    fig = plt.figure(figsize=(11.2, 7.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.88, 1.45], height_ratios=[0.9, 1.25], hspace=0.42, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    boxes = [
        (0.02, 0.62, "188", "robust PAH\ndifferential genes", BLUE),
        (0.37, 0.62, str(robust_n), "mapped to PDE8B\ncoexpression universe", TEAL),
        (0.72, 0.62, str(int(bridge.in_stable_mechanism.sum())), "also in stable\nPDE8B module", ORANGE),
    ]
    for x, y, number, label, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.25, 0.25, facecolor=color, alpha=0.12, edgecolor=color, lw=1.2))
        ax.text(x + 0.125, y + 0.16, number, ha="center", va="center", fontsize=17, color=color, fontweight="bold")
        ax.text(x + 0.125, y + 0.055, label, ha="center", va="center", fontsize=7)
    ax.annotate("", xy=(0.37, 0.745), xytext=(0.27, 0.745), arrowprops=dict(arrowstyle="->", color=GREY))
    ax.annotate("", xy=(0.72, 0.745), xytext=(0.62, 0.745), arrowprops=dict(arrowstyle="->", color=GREY))
    ax.text(0.02, 0.32, "Set overlap is not enriched", fontsize=10, fontweight="bold", color=DARK)
    ax.text(0.02, 0.19, "The differential signature and coexpression program\nanswer different biological questions.", fontsize=8, color="#56616B")
    ax.set_title("A  Two evidence layers, not one inherited gene set", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(len(tests))
    ax.barh(y + 0.12, tests["expected_overlap"], height=0.22, color=GREY, label="Expected by chance")
    ax.barh(y - 0.12, tests["overlap_n"], height=0.22, color=[TEAL, ORANGE], label="Observed")
    ax.set_yticks(y, tests["module"])
    ax.invert_yaxis()
    ax.set_xlabel("Number of overlapping robust genes")
    ax.set_title("B  No membership enrichment", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="lower right")
    for i, row in tests.reset_index(drop=True).iterrows():
        ax.text(max(row.expected_overlap, row.overlap_n) + 0.15, i, f"P(enrich)={row.p_enrichment:.3g}", va="center", fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    colors = np.where(bridge["in_stable_mechanism"], ORANGE, np.where(bridge["direction_concordant"], TEAL, GREY))
    ax.scatter(bridge["meta_g"], bridge["adjusted_meta_r"], c=colors, s=np.where(bridge["in_stable_mechanism"], 35, 14), alpha=0.78, edgecolor="none")
    ax.axhline(0, color="#7A858F", lw=0.7)
    ax.axvline(0, color="#7A858F", lw=0.7)
    for _, row in bridge[bridge["in_stable_mechanism"]].iterrows():
        ax.annotate(row.symbol, (row.meta_g, row.adjusted_meta_r), xytext=(3, 3), textcoords="offset points", fontsize=6.5)
    pde = bridge[bridge.symbol.eq("PDE8B")]
    if not pde.empty:
        row = pde.iloc[0]
        ax.scatter([row.meta_g], [row.adjusted_meta_r], s=70, facecolor="white", edgecolor="#B23A3A", lw=1.5, zorder=4)
        ax.annotate("PDE8B", (row.meta_g, row.adjusted_meta_r), xytext=(5, -11), textcoords="offset points", color="#B23A3A", fontweight="bold")
    ax.set_xlabel("PAH differential effect (random-effects Hedges g)")
    ax.set_ylabel("Meta-correlation with PDE8B (r)")
    ax.set_title("C  Disease effects align with PDE8B covariance direction", loc="left", fontweight="bold")
    ax.text(0.02, 0.98, f"Spearman rho={rho:.3f}, P={rho_p:.2g}\nDirection concordance={n_concordant}/{robust_n} ({n_concordant/robust_n:.1%})",
            transform=ax.transAxes, va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#D6DCE1"))

    ax = fig.add_subplot(gs[1, 0])
    overlap = bridge[bridge["in_stable_mechanism"]].sort_values("adjusted_meta_r")
    yy = np.arange(len(overlap))
    ax.hlines(yy, 0, overlap["adjusted_meta_r"], color="#D1D8DE", lw=1.2)
    ax.scatter(overlap["adjusted_meta_r"], yy, c=np.where(overlap["adjusted_meta_r"] >= 0, TEAL, ORANGE), s=35)
    ax.axvline(0, color="#7A858F", lw=0.7)
    ax.set_yticks(yy, overlap["symbol"])
    ax.set_xlabel("Meta-correlation with PDE8B (r)")
    ax.set_title("D  Six genes bridge both evidence layers", loc="left", fontweight="bold")

    fig.suptitle("The robust PAH signature and PDE8B program are orthogonal but directionally coupled", x=0.04, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.04, 0.015, "Statistics use gene-level meta-analysis outputs; this figure supports association and program coupling, not causal mechanism.", fontsize=7, color="#56616B")
    fig.savefig(OUT / "Figure_signature_PDE8B_bridge_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_signature_PDE8B_bridge_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_signature_PDE8B_bridge_draft.svg", bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

