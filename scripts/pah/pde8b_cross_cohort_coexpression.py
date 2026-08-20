#!/usr/bin/env python3
"""Cross-cohort PDE8B co-expression meta-analysis in human PAH lung tissue."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(r"C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics")
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "pah"))

from project_paths import PROJECT_DATA_ROOT
from formal_meta_and_confounding import bh_adjust, random_effects
from lung_mechanism_meta import (
    prepare_array_cohort,
    prepare_gse208592,
    prepare_gse254617,
)


OUT = PROJECT_DATA_ROOT / "outputs" / "pde8b-coexpression"
OUT.mkdir(parents=True, exist_ok=True)
COHORTS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]


def load_cohorts():
    cohorts = {}
    for gse in ["GSE15197", "GSE113439"]:
        matrix, case, control, qc = prepare_array_cohort(gse)
        cohorts[gse] = (matrix, case, control, qc)
    matrix, case, control, qc = prepare_gse254617()
    cohorts["GSE254617"] = (matrix, case, control, qc)
    matrix, case, control, qc = prepare_gse208592()
    cohorts["GSE208592"] = (matrix, case, control, qc)
    return cohorts


def rank_residual_correlations(matrix, case, control, target="PDE8B"):
    """Spearman-like correlation after removing the case/control group mean."""
    samples = case + control
    x = matrix[samples].astype(float)
    ranks = x.rank(axis=1, method="average").to_numpy(dtype=float).copy()
    y = ranks[x.index.get_loc(target)].copy()
    group = np.array([1] * len(case) + [0] * len(control), dtype=int)

    residual = ranks.copy()
    y_residual = y.copy()
    for value in [0, 1]:
        idx = group == value
        residual[:, idx] -= residual[:, idx].mean(axis=1, keepdims=True)
        y_residual[idx] -= y_residual[idx].mean()

    numerator = residual @ y_residual
    denominator = np.sqrt((residual**2).sum(axis=1) * (y_residual**2).sum())
    corr = np.divide(numerator, denominator, out=np.full(len(x), np.nan), where=denominator > 0)
    corr[x.index.get_loc(target)] = np.nan
    return pd.Series(corr, index=x.index), len(samples)


def pah_only_correlations(matrix, case, target="PDE8B"):
    x = matrix[case].astype(float)
    ranks = x.rank(axis=1, method="average").to_numpy(dtype=float).copy()
    ranks -= ranks.mean(axis=1, keepdims=True)
    y = ranks[x.index.get_loc(target)]
    numerator = ranks @ y
    denominator = np.sqrt((ranks**2).sum(axis=1) * (y**2).sum())
    corr = np.divide(numerator, denominator, out=np.full(len(x), np.nan), where=denominator > 0)
    corr[x.index.get_loc(target)] = np.nan
    return pd.Series(corr, index=x.index), len(case)


def fisher_effect(r, n, covariates):
    clipped = np.clip(np.asarray(r, dtype=float), -0.999999, 0.999999)
    z = np.arctanh(clipped)
    variance = np.full_like(z, 1.0 / max(1, n - covariates - 3), dtype=float)
    variance[~np.isfinite(r)] = np.nan
    return z, variance


def meta_table(correlations, ns, covariates, prefix):
    common = set.intersection(*(set(x.dropna().index) for x in correlations.values()))
    rows = []
    for gene in sorted(common):
        rs = np.array([correlations[gse].loc[gene] for gse in COHORTS], dtype=float)
        zs = []
        variances = []
        for i, gse in enumerate(COHORTS):
            z, v = fisher_effect([rs[i]], ns[gse], covariates)
            zs.append(z[0])
            variances.append(v[0])
        zs = np.array(zs)
        variances = np.array(variances)
        pooled_z, se, lo_z, hi_z, p, tau2, q, i2 = random_effects(zs, variances)
        rows.append({
            "symbol": gene,
            **{f"r_{gse}": rs[i] for i, gse in enumerate(COHORTS)},
            f"{prefix}_same_direction": bool(np.isfinite(rs).all() and len(set(np.sign(rs))) == 1),
            f"{prefix}_meta_z": pooled_z,
            f"{prefix}_meta_r": math.tanh(pooled_z) if np.isfinite(pooled_z) else np.nan,
            f"{prefix}_ci_low_r": math.tanh(lo_z) if np.isfinite(lo_z) else np.nan,
            f"{prefix}_ci_high_r": math.tanh(hi_z) if np.isfinite(hi_z) else np.nan,
            f"{prefix}_p": p,
            f"{prefix}_tau2_z": tau2,
            f"{prefix}_Q": q,
            f"{prefix}_I2": i2,
        })
    result = pd.DataFrame(rows)
    result[f"{prefix}_fdr"] = bh_adjust(result[f"{prefix}_p"])
    return result.sort_values([f"{prefix}_fdr", f"{prefix}_I2"])


def leave_one_cohort_out(adjusted):
    records = []
    for omitted in COHORTS:
        kept = [x for x in COHORTS if x != omitted]
        for gene, row in adjusted.set_index("symbol").iterrows():
            rs = np.array([row[f"r_{gse}"] for gse in kept], dtype=float)
            zs = np.arctanh(np.clip(rs, -0.999999, 0.999999))
            variances = np.array([1.0 / max(1, NS_ADJUSTED[gse] - 4) for gse in kept])
            pooled, se, lo, hi, p, tau2, q, i2 = random_effects(zs, variances)
            records.append({
                "symbol": gene,
                "omitted_cohort": omitted,
                "meta_r": math.tanh(pooled),
                "ci_low_r": math.tanh(lo),
                "ci_high_r": math.tanh(hi),
                "p": p,
                "I2": i2,
            })
    return pd.DataFrame(records)


NS_ADJUSTED = {}


def main():
    cohorts = load_cohorts()
    adjusted_corr, pah_corr = {}, {}
    n_adjusted, n_pah = {}, {}
    qc = {}
    per_cohort = []
    for gse in COHORTS:
        matrix, case, control, detail = cohorts[gse]
        if "PDE8B" not in matrix.index:
            raise RuntimeError(f"PDE8B missing from {gse}")
        adjusted_corr[gse], n_adjusted[gse] = rank_residual_correlations(matrix, case, control)
        pah_corr[gse], n_pah[gse] = pah_only_correlations(matrix, case)
        qc[gse] = {
            "n_pah": len(case), "n_control": len(control), "genes": len(matrix),
            "pde8b_min": float(matrix.loc["PDE8B"].min()),
            "pde8b_max": float(matrix.loc["PDE8B"].max()),
            **detail,
        }
        for gene, r in adjusted_corr[gse].items():
            per_cohort.append({"symbol": gene, "cohort": gse, "analysis": "group_adjusted", "r": r, "n": n_adjusted[gse]})
        for gene, r in pah_corr[gse].items():
            per_cohort.append({"symbol": gene, "cohort": gse, "analysis": "PAH_only", "r": r, "n": n_pah[gse]})

    global NS_ADJUSTED
    NS_ADJUSTED = n_adjusted
    adjusted = meta_table(adjusted_corr, n_adjusted, covariates=1, prefix="adjusted")
    pah_only = meta_table(pah_corr, n_pah, covariates=0, prefix="pah_only")
    merged = adjusted.merge(pah_only, on="symbol", suffixes=("", "_pah"))
    merged["adjusted_robust"] = (
        merged.adjusted_same_direction & (merged.adjusted_fdr < 0.05) &
        (merged.adjusted_I2 < 50) & (merged.adjusted_meta_r.abs() >= 0.30)
    )
    merged["pah_support"] = (
        merged.pah_only_same_direction &
        (np.sign(merged.pah_only_meta_r) == np.sign(merged.adjusted_meta_r)) &
        (merged.pah_only_fdr < 0.05) & (merged.pah_only_I2 < 50) &
        (merged.pah_only_meta_r.abs() >= 0.20)
    )
    merged["mechanism_module"] = merged.adjusted_robust & merged.pah_support
    merged["core_module"] = merged.mechanism_module & (merged.adjusted_meta_r.abs() >= 0.40)
    merged = merged.sort_values(["mechanism_module", "adjusted_fdr", "adjusted_I2"], ascending=[False, True, True])

    pd.DataFrame(per_cohort).to_csv(OUT / "pde8b_per_cohort_correlations.csv", index=False)
    merged.to_csv(OUT / "pde8b_cross_cohort_coexpression_meta.csv", index=False)
    module = merged[merged.mechanism_module].copy()
    module[["symbol", "adjusted_meta_r", "adjusted_fdr", "adjusted_I2", "pah_only_meta_r", "pah_only_fdr", "pah_only_I2"]].to_csv(
        OUT / "pde8b_stable_mechanism_module.csv", index=False
    )
    merged[merged.core_module].to_csv(OUT / "pde8b_high_confidence_core_module.csv", index=False)
    loco = leave_one_cohort_out(adjusted)
    selected = set(module.symbol) | {"PDE8B"}
    loco[loco.symbol.isin(selected)].to_csv(OUT / "pde8b_module_leave_one_cohort_out.csv", index=False)

    summary = {
        "method": "Within-cohort rank correlation after removal of PAH/control group means, followed by Fisher-z random-effects meta-analysis.",
        "cohort_qc": qc,
        "common_tested_genes": int(len(merged)),
        "adjusted_fdr_lt_0.05": int((merged.adjusted_fdr < 0.05).sum()),
        "adjusted_robust": int(merged.adjusted_robust.sum()),
        "stable_mechanism_module": int(merged.mechanism_module.sum()),
        "high_confidence_core_abs_r_ge_0.40": int(merged.core_module.sum()),
        "positive_module": int((module.adjusted_meta_r > 0).sum()),
        "negative_module": int((module.adjusted_meta_r < 0).sum()),
        "top_positive": module.sort_values("adjusted_meta_r", ascending=False).head(30).to_dict("records"),
        "top_negative": module.sort_values("adjusted_meta_r").head(30).to_dict("records"),
        "important_limit": "Co-expression supports a reproducible molecular program but does not establish PDE8B as the causal regulator of correlated genes.",
    }
    (OUT / "pde8b_coexpression_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()




