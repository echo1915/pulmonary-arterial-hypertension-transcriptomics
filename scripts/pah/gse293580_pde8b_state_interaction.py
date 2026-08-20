"""Patient-level PDE8B cell-state interaction and program association."""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


INTERIM = PROJECT_DATA_ROOT / "interim" / "gse293580_scanpy"
PROCESSED = PROJECT_DATA_ROOT / "processed" / "gse293580_scanpy"
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis"
INPUT = INTERIM / "gse293580_phase3_mural_annotated.h5ad"
PATIENT_INPUT = PROCESSED / "gse293580_pde8b_patient_pseudobulk.csv"
PROGRAM_INPUT = PROCESSED / "gse293580_smc_patient_programs.csv"

MODULES = {
    "contractile": ["ACTA2", "TAGLN", "MYH11", "CNN1", "LMOD1", "SMTN", "MYL9", "TPM2"],
    "matrix_remodeling": ["FN1", "COL1A1", "COL1A2", "COL3A1", "THBS2", "MGP", "AEBP1", "CLU"],
    "proliferation": ["MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA", "STMN1"],
    "cAMP_response": ["CREB1", "CREM", "ATF3", "FOS", "JUN", "NR4A1", "DUSP1"],
}
SMC_LABELS = {"contractile_SMC", "modulated_SMC", "SMC_like"}


def exact_signflip_p(values: np.ndarray) -> float:
    """Two-sided exact paired sign-flip permutation test on mean difference."""
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan
    obs = abs(values.mean())
    diffs = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        diffs.append(abs(np.mean(values * np.asarray(signs))))
    return float(np.mean(np.asarray(diffs) >= obs - 1e-12))


def exact_group_permutation_p(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sided exact label-permutation p for a difference in means."""
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    pooled = np.r_[x, y]
    observed = abs(x.mean() - y.mean())
    exceed = 0
    total = 0
    for idx in itertools.combinations(range(len(pooled)), len(x)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        candidate = abs(pooled[mask].mean() - pooled[~mask].mean())
        exceed += candidate >= observed - 1e-12
        total += 1
    return exceed / total


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    """Bias-corrected standardized mean difference (x minus y)."""
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    pooled_sd = np.sqrt(
        ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1))
        / (len(x) + len(y) - 2)
    )
    if pooled_sd == 0:
        return 0.0 if x.mean() == y.mean() else np.sign(x.mean() - y.mean()) * np.inf
    correction = 1 - 3 / (4 * (len(x) + len(y)) - 9)
    return correction * (x.mean() - y.mean()) / pooled_sd


def exact_spearman_p(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Exact permutation p for Spearman rho; limited to n <= 9."""
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 4:
        return np.nan, np.nan
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    x_centered = x_rank - x_rank.mean()
    y_centered = y_rank - y_rank.mean()
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))
    if denominator == 0:
        return np.nan, np.nan
    rho = float(np.dot(x_centered, y_centered) / denominator)
    if len(x) > 7:
        return rho, float(stats.spearmanr(x, y).pvalue)
    exceed = 0
    total = 0
    # Spearman correlation is Pearson correlation of ranks. Permuting the
    # already ranked response is exactly equivalent and avoids thousands of
    # repeated rank/statistics object constructions.
    for perm in itertools.permutations(y_centered.tolist()):
        candidate = float(np.dot(x_centered, np.asarray(perm)) / denominator)
        exceed += abs(candidate) >= abs(rho) - 1e-12
        total += 1
    return rho, exceed / total


def leave_one_out_spearman_range(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Range of Spearman rho after omitting each patient once."""
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 5:
        return np.nan, np.nan
    values = [float(stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic) for i in range(len(x))]
    return float(np.nanmin(values)), float(np.nanmax(values))


def mean_expression(adata, mask: np.ndarray, genes: list[str]) -> float:
    present = [g for g in genes if g in adata.raw.var_names]
    if not present or not mask.any():
        return np.nan
    block = adata.raw[adata.obs_names[mask], present].X
    return float(block.mean())


def bh_adjust(p: pd.Series) -> pd.Series:
    values = p.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=p.index)


def main() -> None:
    patient = pd.read_csv(PATIENT_INPUT)

    # Within-patient localization contrast: SMC minus pericyte detection rate.
    det = patient.pivot_table(
        index=["sample", "condition"], columns="compartment",
        values=["n_cells", "PDE8B_detection_fraction"], aggfunc="first",
    ).reset_index()
    det.columns = ["_".join([str(x) for x in col if str(x)]) if isinstance(col, tuple) else col for col in det.columns]
    det = det.rename(columns={"sample_": "sample", "condition_": "condition"})
    det["SMC_minus_pericyte_detection"] = (
        det["PDE8B_detection_fraction_strict_SMC"] - det["PDE8B_detection_fraction_pericyte"]
    )
    det["eligible_pair"] = (det["n_cells_strict_SMC"] >= 5) & (det["n_cells_pericyte"] >= 5)
    det.to_csv(PROCESSED / "gse293580_pde8b_smc_pericyte_paired.csv", index=False)

    interaction_rows = []
    for condition in ["Donor", "IPAH", "SSc-PAH", "PAH_combined"]:
        if condition == "PAH_combined":
            frame = det.loc[det["condition"].isin(["IPAH", "SSc-PAH"]) & det["eligible_pair"]]
        else:
            frame = det.loc[(det["condition"] == condition) & det["eligible_pair"]]
        values = frame["SMC_minus_pericyte_detection"].dropna().to_numpy()
        interaction_rows.append({
            "condition": condition,
            "n_paired_patients": len(values),
            "mean_SMC_minus_pericyte_detection": values.mean() if len(values) else np.nan,
            "median_SMC_minus_pericyte_detection": np.median(values) if len(values) else np.nan,
            "exact_signflip_p": exact_signflip_p(values),
        })
    interaction = pd.DataFrame(interaction_rows)
    interaction.to_csv(OUT / "gse293580_pde8b_cellstate_paired_tests.csv", index=False)

    # Disease-by-compartment interaction: compare each disease group's paired
    # SMC-minus-pericyte contrast against donors. The patient remains the unit.
    eligible_det = det.loc[det["eligible_pair"]].copy()
    donor_delta = eligible_det.loc[
        eligible_det["condition"] == "Donor", "SMC_minus_pericyte_detection"
    ].to_numpy(float)
    between_rows = []
    for contrast, conditions in {
        "IPAH_vs_Donor": ["IPAH"],
        "SSc-PAH_vs_Donor": ["SSc-PAH"],
        "PAH_combined_vs_Donor": ["IPAH", "SSc-PAH"],
    }.items():
        case_delta = eligible_det.loc[
            eligible_det["condition"].isin(conditions), "SMC_minus_pericyte_detection"
        ].to_numpy(float)
        between_rows.append({
            "contrast": contrast,
            "n_case": len(case_delta),
            "n_donor": len(donor_delta),
            "case_mean_delta": case_delta.mean() if len(case_delta) else np.nan,
            "donor_mean_delta": donor_delta.mean() if len(donor_delta) else np.nan,
            "difference_in_differences": case_delta.mean() - donor_delta.mean()
            if len(case_delta) and len(donor_delta) else np.nan,
            "hedges_g": hedges_g(case_delta, donor_delta),
            "exact_permutation_p": exact_group_permutation_p(case_delta, donor_delta),
        })
    between_interaction = pd.DataFrame(between_rows)
    between_interaction.to_csv(
        OUT / "gse293580_pde8b_disease_by_cellstate_interaction.csv", index=False
    )

    # Reuse the completed patient-program checkpoint when valid. Rebuilding it
    # requires repeated slicing of the large H5AD and is unnecessary for the
    # downstream interaction/sensitivity continuation.
    required_program_columns = {
        "sample", "condition", "n_SMC", "PDE8B_log1p_CPM", "PDE8B_detection_fraction",
        *[f"module_{module}" for module in MODULES],
    }
    if PROGRAM_INPUT.exists():
        programs = pd.read_csv(PROGRAM_INPUT)
        missing = required_program_columns.difference(programs.columns)
        if missing:
            raise ValueError(f"program checkpoint missing columns: {sorted(missing)}")
        program_source = "reused_existing_patient_program_checkpoint"
    else:
        import scanpy as sc

        adata = sc.read_h5ad(INPUT)
        labels = adata.obs["locked_mural_annotation"].astype(str).to_numpy()
        smc = np.isin(labels, list(SMC_LABELS))
        rows = []
        for (sample, condition), frame in adata.obs.groupby(["sample", "condition"], observed=True):
            mask = smc & (adata.obs["sample"].to_numpy() == sample)
            if mask.sum() < 5:
                continue
            row = {"sample": sample, "condition": condition, "n_SMC": int(mask.sum())}
            for module, genes in MODULES.items():
                row[f"module_{module}"] = mean_expression(adata, mask, genes)
            rows.append(row)
        programs = pd.DataFrame(rows)
        pde = patient.loc[
            patient["compartment"] == "strict_SMC",
            ["sample", "PDE8B_log1p_CPM", "PDE8B_detection_fraction"],
        ]
        programs = programs.merge(pde, on="sample", how="left")
        programs.to_csv(PROGRAM_INPUT, index=False)
        program_source = "rebuilt_from_phase3_h5ad"

    program_rows = []
    donor_programs = programs.loc[programs["condition"] == "Donor"]
    for contrast, conditions in {
        "IPAH_vs_Donor": ["IPAH"],
        "SSc-PAH_vs_Donor": ["SSc-PAH"],
        "PAH_combined_vs_Donor": ["IPAH", "SSc-PAH"],
    }.items():
        case_programs = programs.loc[programs["condition"].isin(conditions)]
        for module in MODULES:
            column = f"module_{module}"
            x = case_programs[column].to_numpy(float)
            y = donor_programs[column].to_numpy(float)
            program_rows.append({
                "contrast": contrast,
                "module": module,
                "n_case": int(np.isfinite(x).sum()),
                "n_donor": int(np.isfinite(y).sum()),
                "case_mean": np.nanmean(x),
                "donor_mean": np.nanmean(y),
                "mean_difference": np.nanmean(x) - np.nanmean(y),
                "hedges_g": hedges_g(x, y),
                "exact_permutation_p": exact_group_permutation_p(x, y),
            })
    program_tests = pd.DataFrame(program_rows)
    for contrast in program_tests["contrast"].unique():
        idx = (program_tests["contrast"] == contrast) & program_tests["exact_permutation_p"].notna()
        program_tests.loc[idx, "BH_q_within_contrast"] = bh_adjust(
            program_tests.loc[idx, "exact_permutation_p"]
        )
    program_tests.to_csv(OUT / "gse293580_smc_program_group_tests.csv", index=False)

    corr_rows = []
    # PAH-only is primary for state association; all-sample is sensitivity.
    for scope, frame in [("PAH_only", programs.loc[programs["condition"] != "Donor"]), ("all", programs)]:
        for module in MODULES:
            x = frame["PDE8B_log1p_CPM"].to_numpy(float)
            y = frame[f"module_{module}"].to_numpy(float)
            rho, pvalue = exact_spearman_p(x, y)
            loo_min, loo_max = leave_one_out_spearman_range(x, y)
            corr_rows.append({
                "scope": scope,
                "module": module,
                "n_patients": int((np.isfinite(x) & np.isfinite(y)).sum()),
                "spearman_rho": rho,
                "loo_rho_min": loo_min,
                "loo_rho_max": loo_max,
                "permutation_or_asymptotic_p": pvalue,
            })
    correlations = pd.DataFrame(corr_rows)
    for scope in correlations["scope"].unique():
        idx = (correlations["scope"] == scope) & correlations["permutation_or_asymptotic_p"].notna()
        correlations.loc[idx, "BH_q_within_scope"] = bh_adjust(correlations.loc[idx, "permutation_or_asymptotic_p"])
    correlations.to_csv(OUT / "gse293580_pde8b_smc_program_correlations.csv", index=False)

    summary = {
        "statistical_unit": "patient/sample",
        "minimum_cells_per_compartment": 5,
        "primary_state_test": "within-patient SMC minus pericyte PDE8B detection fraction",
        "disease_by_state_test": "exact comparison of patient-level SMC-minus-pericyte contrasts versus donors",
        "primary_program_scope": "PAH-only strict SMC",
        "program_sensitivity": "leave-one-patient-out Spearman rho range",
        "patient_program_source": program_source,
        "modules": MODULES,
        "caveat": "Exploratory exact tests; very small patient counts constrain inference.",
    }
    (OUT / "gse293580_phase4_state_interaction_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(interaction.to_string(index=False))
    print(between_interaction.to_string(index=False))
    print(program_tests.to_string(index=False))
    print(correlations.to_string(index=False))


if __name__ == "__main__":
    main()
