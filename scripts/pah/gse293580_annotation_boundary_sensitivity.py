"""Patient-level PDE8B sensitivity to locked mural annotation boundaries."""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


INPUT = PROJECT_DATA_ROOT / "interim" / "gse293580_scanpy" / "gse293580_phase3_mural_annotated.h5ad"
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis" / "sensitivity"
MIN_CELLS = 5
BOUNDARIES = {
    "core_SMC": {"contractile_SMC", "modulated_SMC"},
    "strict_SMC_with_like": {"contractile_SMC", "modulated_SMC", "SMC_like"},
    "contractile_SMC": {"contractile_SMC"},
    "modulated_SMC": {"modulated_SMC"},
    "SMC_like": {"SMC_like"},
    "pericyte": {"pericyte"},
}
CONTRASTS = {
    "IPAH_vs_Donor": ("IPAH",),
    "SSc-PAH_vs_Donor": ("SSc-PAH",),
    "PAH_combined_vs_Donor": ("IPAH", "SSc-PAH"),
}


def exact_p(case: np.ndarray, donor: np.ndarray) -> float:
    if len(case) < 2 or len(donor) < 2:
        return np.nan
    pooled = np.r_[case, donor]
    observed = abs(case.mean() - donor.mean())
    values = []
    for idx in itertools.combinations(range(len(pooled)), len(case)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        values.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
    return float(np.mean(np.asarray(values) >= observed - 1e-12))


def hedges_g(case: np.ndarray, donor: np.ndarray) -> float:
    if len(case) < 2 or len(donor) < 2:
        return np.nan
    pooled_sd = np.sqrt(
        ((len(case) - 1) * case.var(ddof=1) + (len(donor) - 1) * donor.var(ddof=1))
        / (len(case) + len(donor) - 2)
    )
    if pooled_sd == 0:
        return 0.0 if case.mean() == donor.mean() else np.sign(case.mean() - donor.mean()) * np.inf
    return (case.mean() - donor.mean()) / pooled_sd * (1 - 3 / (4 * (len(case) + len(donor)) - 9))


def loo(case: np.ndarray, donor: np.ndarray) -> tuple[int, int, float, float]:
    if not len(case) or not len(donor):
        return 0, 0, np.nan, np.nan
    sign = np.sign(case.mean() - donor.mean())
    deltas = []
    if len(case) > 1:
        deltas.extend(np.delete(case, i).mean() - donor.mean() for i in range(len(case)))
    if len(donor) > 1:
        deltas.extend(case.mean() - np.delete(donor, i).mean() for i in range(len(donor)))
    values = np.asarray(deltas)
    return int(np.sum(np.sign(values) == sign)), len(values), float(values.min()), float(values.max())


def main() -> None:
    import anndata as ad

    adata = ad.read_h5ad(INPUT, backed="r")
    obs = adata.obs[["sample", "condition", "locked_mural_annotation"]].copy()
    pde8b_idx = adata.var_names.get_loc("PDE8B")
    counts = adata.layers["counts"][:]
    adata.file.close()
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    library_per_cell = np.asarray(counts.sum(axis=1)).ravel()
    pde8b_per_cell = np.asarray(counts[:, pde8b_idx].toarray()).ravel()

    patient_rows = []
    samples = obs[["sample", "condition"]].drop_duplicates()
    labels = obs["locked_mural_annotation"].astype(str).to_numpy()
    sample_values = obs["sample"].to_numpy()
    for boundary, included in BOUNDARIES.items():
        boundary_mask = np.isin(labels, list(included))
        for sample, condition in samples.itertuples(index=False):
            take = boundary_mask & (sample_values == sample)
            n_cells = int(take.sum())
            library = float(library_per_cell[take].sum())
            pde_counts = float(pde8b_per_cell[take].sum())
            detected = int((pde8b_per_cell[take] > 0).sum())
            cpm = 1e6 * pde_counts / library if library > 0 else np.nan
            patient_rows.append({
                "boundary": boundary,
                "sample": sample,
                "condition": condition,
                "n_cells": n_cells,
                "library_counts": library,
                "PDE8B_counts": pde_counts,
                "PDE8B_detected_cells": detected,
                "PDE8B_detection_fraction": detected / n_cells if n_cells else np.nan,
                "PDE8B_log1p_CPM": np.log1p(cpm) if np.isfinite(cpm) else np.nan,
            })
    patient = pd.DataFrame(patient_rows)

    test_rows = []
    eligible = patient.loc[patient["n_cells"] >= MIN_CELLS]
    for boundary in BOUNDARIES:
        block = eligible.loc[eligible["boundary"] == boundary]
        donors = block.loc[block["condition"] == "Donor"]
        for contrast, conditions in CONTRASTS.items():
            cases = block.loc[block["condition"].isin(conditions)]
            for metric in ("PDE8B_log1p_CPM", "PDE8B_detection_fraction"):
                case = cases[metric].dropna().to_numpy(float)
                donor = donors[metric].dropna().to_numpy(float)
                concordant, iterations, loo_min, loo_max = loo(case, donor)
                test_rows.append({
                    "boundary": boundary,
                    "contrast": contrast,
                    "metric": metric,
                    "minimum_cells": MIN_CELLS,
                    "n_case": len(case),
                    "n_donor": len(donor),
                    "mean_difference": case.mean() - donor.mean() if len(case) and len(donor) else np.nan,
                    "hedges_g": hedges_g(case, donor),
                    "exact_permutation_p": exact_p(case, donor),
                    "loo_direction_concordant_n": concordant,
                    "loo_iterations_n": iterations,
                    "loo_difference_min": loo_min,
                    "loo_difference_max": loo_max,
                })
    tests = pd.DataFrame(test_rows)
    OUT.mkdir(parents=True, exist_ok=True)
    patient.to_csv(OUT / "gse293580_annotation_boundary_patient_values.csv", index=False)
    tests.to_csv(OUT / "gse293580_annotation_boundary_tests.csv", index=False)

    primary = tests.loc[
        (tests["contrast"] == "PAH_combined_vs_Donor")
        & (tests["metric"] == "PDE8B_log1p_CPM")
        & tests["boundary"].isin(["core_SMC", "strict_SMC_with_like"])
    ].set_index("boundary")
    summary = {
        "statistical_unit": "patient/sample",
        "minimum_cells": MIN_CELLS,
        "locked_boundaries": {key: sorted(value) for key, value in BOUNDARIES.items()},
        "core_SMC_hedges_g": float(primary.loc["core_SMC", "hedges_g"]),
        "strict_SMC_with_like_hedges_g": float(primary.loc["strict_SMC_with_like", "hedges_g"]),
        "core_and_strict_direction_concordant": bool(
            np.sign(primary.loc["core_SMC", "hedges_g"])
            == np.sign(primary.loc["strict_SMC_with_like", "hedges_g"])
        ),
        "boundary": "Sensitivity uses locked labels only; no reannotation or reclustering was performed.",
    }
    (OUT / "gse293580_annotation_boundary_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(primary.to_string())


if __name__ == "__main__":
    main()
