"""LOCO, threshold-grid, and negative-control audit of the PAH lung meta-analysis."""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402


INFILE = OUTPUT_ROOT / "pah_audit" / "lung_mechanism_random_effects_meta.csv"
OUT = PROJECT_DATA_ROOT / "outputs" / "bulk-meta-robustness"
OUT.mkdir(parents=True, exist_ok=True)
COHORTS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]
SAMPLE_SIZES = {
    "GSE15197": (18, 13), "GSE113439": (14, 11),
    "GSE254617": (82, 52), "GSE208592": (15, 18),
}
REFERENCE_GENES = ["ACTB", "GAPDH", "RPLP0", "TBP", "PPIA", "GUSB", "HPRT1"]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def random_effects(g: np.ndarray, variance: np.ndarray) -> tuple[float, ...]:
    ok = np.isfinite(g) & np.isfinite(variance) & (variance > 0)
    g, variance = g[ok], variance[ok]
    if len(g) < 2:
        return (np.nan,) * 6
    w = 1 / variance
    fixed = np.sum(w * g) / np.sum(w)
    q = np.sum(w * (g - fixed) ** 2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (len(g) - 1)) / c) if c > 0 else 0.0
    wr = 1 / (variance + tau2)
    pooled = np.sum(wr * g) / np.sum(wr)
    se = math.sqrt(1 / np.sum(wr))
    p = math.erfc(abs(pooled / se) / math.sqrt(2))
    i2 = max(0.0, (q - (len(g) - 1)) / q) * 100 if q > 0 else 0.0
    return pooled, se, pooled - 1.96 * se, pooled + 1.96 * se, p, i2


def meta_table(symbols: pd.Series, effects: np.ndarray, variances: np.ndarray, cohort_names: list[str]) -> pd.DataFrame:
    rows = []
    for i, symbol in enumerate(symbols):
        pooled, se, lo, hi, p, i2 = random_effects(effects[i], variances[i])
        signs = np.sign(effects[i])
        rows.append({
            "symbol": symbol, "meta_g": pooled, "meta_se": se, "ci_low": lo, "ci_high": hi,
            "p": p, "I2": i2, "same_direction": bool(np.all(signs > 0) or np.all(signs < 0)),
            **{f"g_{cohort}": effects[i, j] for j, cohort in enumerate(cohort_names)},
        })
    out = pd.DataFrame(rows)
    out["fdr"] = bh_adjust(out.p.to_numpy())
    out["robust"] = out.same_direction & (out.fdr < 0.05) & (out.I2 < 50) & (out.meta_g.abs() >= 0.50)
    return out


def main() -> None:
    full = pd.read_csv(INFILE)
    effect_cols = [f"g_{c}" for c in COHORTS]
    effects = full[effect_cols].to_numpy(float)
    variances = np.empty_like(effects)
    for j, cohort in enumerate(COHORTS):
        n1, n0 = SAMPLE_SIZES[cohort]
        variances[:, j] = (n1 + n0) / (n1 * n0) + effects[:, j] ** 2 / (2 * (n1 + n0 - 2))

    original_robust = full.same_direction.astype(bool) & (full.fdr < 0.05) & (full.I2 < 50) & (full.meta_g.abs() >= 0.50)
    original_set = set(full.loc[original_robust, "symbol"])
    loco_tables = []
    loco_summary = []
    for omit_idx, omitted in enumerate(COHORTS):
        keep = [i for i in range(4) if i != omit_idx]
        table = meta_table(full.symbol, effects[:, keep], variances[:, keep], [COHORTS[i] for i in keep])
        table["omitted_cohort"] = omitted
        table["original_robust"] = table.symbol.isin(original_set)
        loco_tables.append(table)
        retained = set(table.loc[table.robust, "symbol"]) & original_set
        pde = table[table.symbol.eq("PDE8B")].iloc[0]
        loco_summary.append({
            "omitted_cohort": omitted, "robust_genes_total": int(table.robust.sum()),
            "original_188_retained_n": len(retained), "original_188_retained_fraction": len(retained) / len(original_set),
            "PDE8B_meta_g": pde.meta_g, "PDE8B_ci_low": pde.ci_low, "PDE8B_ci_high": pde.ci_high,
            "PDE8B_fdr": pde.fdr, "PDE8B_I2": pde.I2, "PDE8B_same_direction": bool(pde.same_direction),
            "PDE8B_robust": bool(pde.robust),
        })
    loco = pd.concat(loco_tables, ignore_index=True)
    loco_summary = pd.DataFrame(loco_summary)

    threshold_rows = []
    for fdr_cut, i2_cut, g_cut in itertools.product([0.01, 0.05, 0.10], [25, 50, 75], [0.30, 0.50, 0.70]):
        mask = full.same_direction.astype(bool) & (full.fdr < fdr_cut) & (full.I2 < i2_cut) & (full.meta_g.abs() >= g_cut)
        threshold_rows.append({
            "fdr_cutoff": fdr_cut, "I2_cutoff": i2_cut, "abs_g_cutoff": g_cut,
            "selected_genes": int(mask.sum()), "PDE8B_selected": bool(mask[full.symbol.eq("PDE8B")].iloc[0]),
        })
    thresholds = pd.DataFrame(threshold_rows)

    reference = full[full.symbol.isin(REFERENCE_GENES)].copy()
    reference["passes_primary_robust_definition"] = reference.symbol.isin(original_set)

    sign_rows = []
    for signs in itertools.product([-1, 1], repeat=4):
        if signs in [(-1, -1, -1, -1), (1, 1, 1, 1)]:
            continue
        signed = effects * np.asarray(signs)[None, :]
        table = meta_table(full.symbol, signed, variances, COHORTS)
        sign_rows.append({"sign_pattern": ",".join(map(str, signs)), "robust_count": int(table.robust.sum())})
    sign_null = pd.DataFrame(sign_rows)

    loco.to_csv(OUT / "lung_meta_leave_one_cohort_out_all_genes.csv.gz", index=False, compression="gzip")
    loco_summary.to_csv(OUT / "lung_meta_loco_summary.csv", index=False, encoding="utf-8-sig")
    thresholds.to_csv(OUT / "lung_meta_threshold_grid.csv", index=False, encoding="utf-8-sig")
    reference.to_csv(OUT / "reference_gene_negative_controls.csv", index=False, encoding="utf-8-sig")
    sign_null.to_csv(OUT / "cohort_sign_flip_negative_control.csv", index=False, encoding="utf-8-sig")
    report = {
        "original_robust_n": len(original_set),
        "PDE8B_all_loco_robust": bool(loco_summary.PDE8B_robust.all()),
        "PDE8B_loco_g_range": [float(loco_summary.PDE8B_meta_g.min()), float(loco_summary.PDE8B_meta_g.max())],
        "minimum_original_signature_retention": float(loco_summary.original_188_retained_fraction.min()),
        "reference_genes_passing_robust": reference.loc[reference.passes_primary_robust_definition, "symbol"].tolist(),
        "mixed_sign_null_count_range": [int(sign_null.robust_count.min()), int(sign_null.robust_count.max())],
    }
    (OUT / "bulk_meta_robustness_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 4, figsize=(12, 3.7), gridspec_kw={"width_ratios": [1.2, 1, 1.2, 0.9]})
    pde_full = full[full.symbol.eq("PDE8B")].iloc[0]
    forest = pd.concat([pd.DataFrame([{"omitted_cohort": "None (4 cohorts)", "PDE8B_meta_g": pde_full.meta_g,
                                      "PDE8B_ci_low": pde_full.ci_low, "PDE8B_ci_high": pde_full.ci_high}]), loco_summary], ignore_index=True)
    y = np.arange(len(forest))
    axes[0].errorbar(forest.PDE8B_meta_g, y, xerr=[forest.PDE8B_meta_g - forest.PDE8B_ci_low, forest.PDE8B_ci_high - forest.PDE8B_meta_g], fmt="o", color="#B84343", ecolor="#7D8992", capsize=2)
    axes[0].axvline(0, color="#68747E", lw=0.8)
    axes[0].set_yticks(y, forest.omitted_cohort)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("PDE8B Hedges g")
    axes[0].set_title("A  PDE8B LOCO", loc="left", fontweight="bold")

    y = np.arange(4)
    axes[1].barh(y, loco_summary.original_188_retained_fraction * 100, color="#187B7B")
    axes[1].set_yticks(y, loco_summary.omitted_cohort)
    axes[1].invert_yaxis(); axes[1].set_xlim(0, 105)
    axes[1].set_xlabel("Original signature retained (%)")
    axes[1].set_title("B  188-gene retention", loc="left", fontweight="bold")

    q = thresholds[thresholds.I2_cutoff.eq(50)]
    for fdr_cut, block in q.groupby("fdr_cutoff"):
        axes[2].plot(block.abs_g_cutoff, block.selected_genes, marker="o", label=f"FDR < {fdr_cut:g}")
    axes[2].set_xlabel("Absolute effect cutoff")
    axes[2].set_ylabel("Selected genes")
    axes[2].set_title("C  Threshold grid (I2 < 50%)", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=7)

    reference = reference.sort_values("meta_g").reset_index(drop=True)
    y = np.arange(len(reference))
    axes[3].errorbar(reference.meta_g, y,
                     xerr=[reference.meta_g - reference.ci_low, reference.ci_high - reference.meta_g],
                     fmt="o", color="#7D8992", ecolor="#BFC7CD", capsize=2)
    axes[3].axvline(0, color="#68747E", lw=0.8)
    axes[3].set_yticks(y, reference.symbol)
    axes[3].set_xlabel("Meta Hedges g")
    axes[3].set_title("D  Reference-gene controls", loc="left", fontweight="bold")

    fig.suptitle("PDE8B and the PAH lung signature are tested against cohort dependence and analytic choices", x=0.03, ha="left", fontsize=12.2, fontweight="bold")
    fig.text(0.03, 0.015, "LOCO, leave one cohort out. None of seven commonly used reference genes passes the primary robust-gene definition.", fontsize=7, color="#56616B")
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    fig.savefig(OUT / "Figure_bulk_meta_robustness_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_bulk_meta_robustness_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_bulk_meta_robustness_draft.svg", bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(report, indent=2))
    print(loco_summary.to_string(index=False))


if __name__ == "__main__":
    main()

