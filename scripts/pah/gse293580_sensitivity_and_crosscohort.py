"""Sensitivity audit for GSE293580 and directional replication in GSE210248.

This script continues from existing patient-level checkpoints. It does not
reload single-cell matrices, recluster cells, or change locked annotations.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


PROCESSED = PROJECT_DATA_ROOT / "processed" / "gse293580_scanpy"
GSE293_OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis"
GSE210_STATE = PROJECT_DATA_ROOT / "outputs" / "smc-composition-within-state"
OUT = GSE293_OUT / "sensitivity"

PATIENT_INPUT = PROCESSED / "gse293580_pde8b_patient_pseudobulk.csv"
PROGRAM_INPUT = PROCESSED / "gse293580_smc_patient_programs.csv"
GSE210_EFFECT_INPUT = GSE210_STATE / "gse210248_smc_composition_within_state_effects.csv"

THRESHOLDS = (3, 5, 10, 20)
PROGRAM_THRESHOLDS = (5, 10, 20)
CONTRASTS = {
    "IPAH_vs_Donor": ("IPAH",),
    "SSc-PAH_vs_Donor": ("SSc-PAH",),
    "PAH_combined_vs_Donor": ("IPAH", "SSc-PAH"),
}
PROGRAMS = ("contractile", "matrix_remodeling", "proliferation", "cAMP_response")


def exact_group_permutation_p(case: np.ndarray, donor: np.ndarray) -> float:
    case = case[np.isfinite(case)]
    donor = donor[np.isfinite(donor)]
    if len(case) < 2 or len(donor) < 2:
        return np.nan
    pooled = np.r_[case, donor]
    observed = abs(case.mean() - donor.mean())
    exceed = 0
    total = 0
    for idx in itertools.combinations(range(len(pooled)), len(case)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        candidate = abs(pooled[mask].mean() - pooled[~mask].mean())
        exceed += candidate >= observed - 1e-12
        total += 1
    return exceed / total


def hedges_g(case: np.ndarray, donor: np.ndarray) -> float:
    case = case[np.isfinite(case)]
    donor = donor[np.isfinite(donor)]
    if len(case) < 2 or len(donor) < 2:
        return np.nan
    pooled_sd = np.sqrt(
        ((len(case) - 1) * case.var(ddof=1) + (len(donor) - 1) * donor.var(ddof=1))
        / (len(case) + len(donor) - 2)
    )
    if pooled_sd == 0:
        return 0.0 if case.mean() == donor.mean() else np.sign(case.mean() - donor.mean()) * np.inf
    correction = 1 - 3 / (4 * (len(case) + len(donor)) - 9)
    return correction * (case.mean() - donor.mean()) / pooled_sd


def leave_one_out_direction(case: np.ndarray, donor: np.ndarray) -> tuple[int, int, float, float]:
    case = case[np.isfinite(case)]
    donor = donor[np.isfinite(donor)]
    if not len(case) or not len(donor):
        return 0, 0, np.nan, np.nan
    full_delta = case.mean() - donor.mean()
    full_sign = np.sign(full_delta)
    deltas = []
    if len(case) > 1:
        deltas.extend(np.delete(case, i).mean() - donor.mean() for i in range(len(case)))
    if len(donor) > 1:
        deltas.extend(case.mean() - np.delete(donor, i).mean() for i in range(len(donor)))
    if not deltas:
        return 0, 0, np.nan, np.nan
    values = np.asarray(deltas, dtype=float)
    concordant = int(np.sum(np.sign(values) == full_sign))
    return concordant, len(values), float(values.min()), float(values.max())


def comparison_record(
    case: np.ndarray,
    donor: np.ndarray,
    **labels: object,
) -> dict[str, object]:
    case = case[np.isfinite(case)]
    donor = donor[np.isfinite(donor)]
    concordant, iterations, loo_min, loo_max = leave_one_out_direction(case, donor)
    return {
        **labels,
        "n_case": len(case),
        "n_donor": len(donor),
        "case_mean": case.mean() if len(case) else np.nan,
        "donor_mean": donor.mean() if len(donor) else np.nan,
        "mean_difference": case.mean() - donor.mean() if len(case) and len(donor) else np.nan,
        "hedges_g": hedges_g(case, donor),
        "exact_permutation_p": exact_group_permutation_p(case, donor),
        "loo_direction_concordant_n": concordant,
        "loo_iterations_n": iterations,
        "loo_difference_min": loo_min,
        "loo_difference_max": loo_max,
    }


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return pd.Series(out, index=values.index)


def main() -> None:
    patient = pd.read_csv(PATIENT_INPUT)
    programs = pd.read_csv(PROGRAM_INPUT)
    gse210 = pd.read_csv(GSE210_EFFECT_INPUT)
    OUT.mkdir(parents=True, exist_ok=True)

    coverage_rows = []
    expression_rows = []
    interaction_rows = []
    program_rows = []

    for threshold in THRESHOLDS:
        eligible = patient.loc[patient["n_cells"] >= threshold].copy()
        for (compartment, condition), block in patient.groupby(["compartment", "condition"]):
            coverage_rows.append({
                "analysis": "PDE8B_compartment",
                "threshold": threshold,
                "compartment": compartment,
                "condition": condition,
                "eligible_patients": int((block["n_cells"] >= threshold).sum()),
                "total_patients": int(len(block)),
                "min_cells": int(block["n_cells"].min()),
                "median_cells": float(block["n_cells"].median()),
                "max_cells": int(block["n_cells"].max()),
            })

        for compartment in ("strict_SMC", "pericyte", "broad_mural"):
            block = eligible.loc[eligible["compartment"] == compartment]
            donors = block.loc[block["condition"] == "Donor"]
            for contrast, conditions in CONTRASTS.items():
                cases = block.loc[block["condition"].isin(conditions)]
                for metric in ("PDE8B_log1p_CPM", "PDE8B_detection_fraction"):
                    expression_rows.append(comparison_record(
                        cases[metric].to_numpy(float),
                        donors[metric].to_numpy(float),
                        threshold=threshold,
                        compartment=compartment,
                        contrast=contrast,
                        metric=metric,
                    ))

        paired = patient.pivot_table(
            index=["sample", "condition"],
            columns="compartment",
            values=["n_cells", "PDE8B_detection_fraction"],
            aggfunc="first",
        ).reset_index()
        paired.columns = [
            "_".join(str(x) for x in col if str(x)) if isinstance(col, tuple) else col
            for col in paired.columns
        ]
        paired = paired.rename(columns={"sample_": "sample", "condition_": "condition"})
        paired["delta"] = (
            paired["PDE8B_detection_fraction_strict_SMC"]
            - paired["PDE8B_detection_fraction_pericyte"]
        )
        paired = paired.loc[
            (paired["n_cells_strict_SMC"] >= threshold)
            & (paired["n_cells_pericyte"] >= threshold)
        ]
        donor_delta = paired.loc[paired["condition"] == "Donor", "delta"].to_numpy(float)
        for contrast, conditions in CONTRASTS.items():
            case_delta = paired.loc[paired["condition"].isin(conditions), "delta"].to_numpy(float)
            interaction_rows.append(comparison_record(
                case_delta,
                donor_delta,
                threshold=threshold,
                contrast=contrast,
                metric="SMC_minus_pericyte_detection",
            ))

        # The persisted program checkpoint was originally created with a
        # minimum of five SMCs, so thresholds below five are not recoverable
        # without reopening the H5AD and are deliberately not imputed.
        if threshold in PROGRAM_THRESHOLDS:
            eligible_programs = programs.loc[programs["n_SMC"] >= threshold]
            for condition, block in programs.groupby("condition"):
                coverage_rows.append({
                    "analysis": "SMC_program",
                    "threshold": threshold,
                    "compartment": "strict_SMC",
                    "condition": condition,
                    "eligible_patients": int((block["n_SMC"] >= threshold).sum()),
                    "total_patients": int(len(block)),
                    "min_cells": int(block["n_SMC"].min()),
                    "median_cells": float(block["n_SMC"].median()),
                    "max_cells": int(block["n_SMC"].max()),
                })
            donor_programs = eligible_programs.loc[eligible_programs["condition"] == "Donor"]
            for contrast, conditions in CONTRASTS.items():
                case_programs = eligible_programs.loc[eligible_programs["condition"].isin(conditions)]
                for program in PROGRAMS:
                    column = f"module_{program}"
                    program_rows.append(comparison_record(
                        case_programs[column].to_numpy(float),
                        donor_programs[column].to_numpy(float),
                        threshold=threshold,
                        contrast=contrast,
                        program=program,
                    ))

    coverage = pd.DataFrame(coverage_rows)
    expression = pd.DataFrame(expression_rows)
    interaction = pd.DataFrame(interaction_rows)
    program_tests = pd.DataFrame(program_rows)

    for frame, family in ((expression, ["threshold", "contrast", "metric"]),
                          (program_tests, ["threshold", "contrast"])):
        frame["BH_q_within_family"] = np.nan
        for _, idx in frame.groupby(family).groups.items():
            valid = frame.loc[idx, "exact_permutation_p"].notna()
            valid_idx = frame.loc[idx].index[valid]
            if len(valid_idx):
                frame.loc[valid_idx, "BH_q_within_family"] = bh_adjust(
                    frame.loc[valid_idx, "exact_permutation_p"]
                )

    # Cross-cohort comparison uses standardized effects/directions only because
    # the underlying module-score scales and SMC definitions differ by cohort.
    cross_rows = []
    mapping = {
        "PDE8B": ("PDE8B", "PDE8B"),
        "Contractile": ("Contractile", "contractile"),
        "ECM_remodeling": ("Synthetic/ECM", "matrix_remodeling"),
    }
    gse293_primary = pd.concat([
        expression.loc[
            (expression["threshold"] == 5)
            & (expression["compartment"] == "strict_SMC")
            & (expression["contrast"] == "PAH_combined_vs_Donor")
            & (expression["metric"] == "PDE8B_log1p_CPM")
        ].assign(program="PDE8B"),
        program_tests.loc[
            (program_tests["threshold"] == 5)
            & (program_tests["contrast"] == "PAH_combined_vs_Donor")
        ],
    ], ignore_index=True, sort=False)
    for harmonized, (gse210_feature, gse293_feature) in mapping.items():
        for cell_type in ("SMC 1", "SMC 2"):
            row210 = gse210.loc[
                (gse210["layer"] == "Within-state")
                & (gse210["cell_type"] == cell_type)
                & (gse210["feature"] == gse210_feature)
            ].iloc[0]
            row293 = gse293_primary.loc[gse293_primary["program"] == gse293_feature].iloc[0]
            cross_rows.append({
                "harmonized_feature": harmonized,
                "gse210248_cell_type": cell_type,
                "gse210248_hedges_g": row210["hedges_g"],
                "gse293580_hedges_g": row293["hedges_g"],
                "direction_concordant": bool(
                    np.sign(row210["hedges_g"]) == np.sign(row293["hedges_g"])
                ),
                "gse210248_exact_p": row210["exact_permutation_p"],
                "gse293580_exact_p": row293["exact_permutation_p"],
                "interpretation": "directional_only_noncomparable_score_scales",
            })
    cross = pd.DataFrame(cross_rows)

    coverage.to_csv(OUT / "gse293580_threshold_sample_coverage.csv", index=False)
    expression.to_csv(OUT / "gse293580_pde8b_threshold_sensitivity.csv", index=False)
    interaction.to_csv(OUT / "gse293580_smc_pericyte_interaction_sensitivity.csv", index=False)
    program_tests.to_csv(OUT / "gse293580_smc_program_threshold_sensitivity.csv", index=False)
    cross.to_csv(OUT / "gse210248_gse293580_directional_replication.csv", index=False)

    primary_expression = expression.loc[
        (expression["threshold"] == 5)
        & (expression["compartment"] == "strict_SMC")
        & (expression["contrast"] == "PAH_combined_vs_Donor")
        & (expression["metric"] == "PDE8B_log1p_CPM")
    ].iloc[0]
    primary_interaction = interaction.loc[
        (interaction["threshold"] == 5)
        & (interaction["contrast"] == "PAH_combined_vs_Donor")
    ].iloc[0]
    summary = {
        "statistical_unit": "patient/sample",
        "checkpoint_strategy": "existing patient pseudobulk and program matrices reused; no reclustering",
        "thresholds": list(THRESHOLDS),
        "program_thresholds": list(PROGRAM_THRESHOLDS),
        "program_threshold_floor_reason": "patient program checkpoint was created with >=5 strict SMCs",
        "primary_threshold": 5,
        "primary_combined_PAH_SMC_PDE8B_hedges_g": float(primary_expression["hedges_g"]),
        "primary_combined_PAH_SMC_PDE8B_exact_p": float(primary_expression["exact_permutation_p"]),
        "primary_disease_by_cellstate_hedges_g": float(primary_interaction["hedges_g"]),
        "primary_disease_by_cellstate_exact_p": float(primary_interaction["exact_permutation_p"]),
        "cross_cohort_direction_concordance": cross.groupby("harmonized_feature")["direction_concordant"].all().to_dict(),
        "inference_boundary": "Exploratory small-n sensitivity analysis; direction agreement is not causal evidence.",
    }
    (OUT / "gse293580_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print("\nPrimary threshold sensitivity:")
    print(expression.loc[
        (expression["compartment"] == "strict_SMC")
        & (expression["contrast"] == "PAH_combined_vs_Donor")
        & (expression["metric"] == "PDE8B_log1p_CPM")
    ].to_string(index=False))
    print("\nCross-cohort directional replication:")
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
