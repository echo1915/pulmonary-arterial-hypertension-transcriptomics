"""Patient-level PDE8B validation after locking GSE293580 mural annotations."""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


INTERIM = PROJECT_DATA_ROOT / "interim" / "gse293580_scanpy"
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis"
PROCESSED = PROJECT_DATA_ROOT / "processed" / "gse293580_scanpy"
INPUT = INTERIM / "gse293580_phase2_mural_clustered.h5ad"

ANNOTATION = {
    "0": "non_mural_contaminant",
    "1": "low_quality_contaminant",
    "2": "pericyte",
    "3": "contractile_SMC",
    "4": "modulated_SMC",
    "5": "SMC_like",
}
SMC_LABELS = {"contractile_SMC", "modulated_SMC", "SMC_like"}


def exact_label_permutation(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sided exact permutation p for mean difference."""
    pooled = np.r_[x, y]
    nx = len(x)
    observed = abs(x.mean() - y.mean())
    exceed = 0
    total = 0
    for idx in itertools.combinations(range(len(pooled)), nx):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        diff = abs(pooled[mask].mean() - pooled[~mask].mean())
        exceed += diff >= observed - 1e-12
        total += 1
    return exceed / total


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    pooled_sd = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    if pooled_sd == 0:
        return 0.0 if x.mean() == y.mean() else np.sign(x.mean() - y.mean()) * np.inf
    d = (x.mean() - y.mean()) / pooled_sd
    correction = 1 - 3 / (4 * (nx + ny) - 9)
    return d * correction


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return pd.Series(out, index=values.index)


def aggregate(adata, mask: np.ndarray, compartment: str) -> pd.DataFrame:
    gene_idx = adata.var_names.get_loc("PDE8B")
    rows = []
    all_samples = adata.obs[["sample", "condition"]].drop_duplicates()
    for sample, condition in all_samples.itertuples(index=False):
        take = mask & (adata.obs["sample"].to_numpy() == sample)
        n_cells = int(take.sum())
        if n_cells:
            block = adata.layers["counts"][take]
            total = float(block.sum())
            gene = block[:, gene_idx]
            pde8b_sum = float(gene.sum())
            gene_values = gene.toarray().ravel() if sparse.issparse(gene) else np.asarray(gene).ravel()
            detected = int((gene_values > 0).sum())
        else:
            total = pde8b_sum = 0.0
            detected = 0
        cpm = 1e6 * pde8b_sum / total if total > 0 else np.nan
        rows.append({
            "compartment": compartment,
            "sample": sample,
            "condition": condition,
            "n_cells": n_cells,
            "library_counts": total,
            "PDE8B_counts": pde8b_sum,
            "PDE8B_detected_cells": detected,
            "PDE8B_detection_fraction": detected / n_cells if n_cells else np.nan,
            "PDE8B_CPM": cpm,
            "PDE8B_log1p_CPM": np.log1p(cpm) if np.isfinite(cpm) else np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    adata = sc.read_h5ad(INPUT)
    counts = adata.layers["counts"]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    integer_max_error = float(np.max(np.abs(values - np.rint(values)))) if values.size else 0.0
    if integer_max_error > 1e-8:
        raise ValueError(f"counts layer is not integer-like; max error={integer_max_error}")

    cluster = adata.obs["mural_leiden_0.4"].astype(str)
    adata.obs["locked_mural_annotation"] = cluster.map(ANNOTATION).astype("category")
    labels = adata.obs["locked_mural_annotation"].astype(str)
    strict_smc = labels.isin(SMC_LABELS).to_numpy()
    pericyte = (labels == "pericyte").to_numpy()
    broad_mural = (strict_smc | pericyte)

    tables = [
        aggregate(adata, strict_smc, "strict_SMC"),
        aggregate(adata, pericyte, "pericyte"),
        aggregate(adata, broad_mural, "broad_mural"),
    ]
    patient = pd.concat(tables, ignore_index=True)
    patient.to_csv(PROCESSED / "gse293580_pde8b_patient_pseudobulk.csv", index=False)

    # Primary expression rule is fixed before testing: at least 5 cells in the
    # compartment. This is reported, not tuned to obtain significance.
    eligible = patient.loc[patient["n_cells"] >= 5].copy()
    test_rows = []
    for compartment in ["strict_SMC", "pericyte", "broad_mural"]:
        frame = eligible.loc[eligible["compartment"] == compartment]
        for case in ["IPAH", "SSc-PAH"]:
            x = frame.loc[frame["condition"] == case, "PDE8B_log1p_CPM"].dropna().to_numpy()
            y = frame.loc[frame["condition"] == "Donor", "PDE8B_log1p_CPM"].dropna().to_numpy()
            test_rows.append({
                "compartment": compartment,
                "contrast": f"{case}_vs_Donor",
                "n_case": len(x),
                "n_donor": len(y),
                "case_mean_log1p_CPM": x.mean() if len(x) else np.nan,
                "donor_mean_log1p_CPM": y.mean() if len(y) else np.nan,
                "mean_difference": x.mean() - y.mean() if len(x) and len(y) else np.nan,
                "hedges_g": hedges_g(x, y),
                "exact_permutation_p": exact_label_permutation(x, y) if len(x) >= 2 and len(y) >= 2 else np.nan,
                "interpretability": "formal_exploratory" if len(x) >= 3 and len(y) >= 3 else "descriptive_low_n",
            })
    tests = pd.DataFrame(test_rows)
    valid = tests["exact_permutation_p"].notna()
    tests.loc[valid, "BH_q"] = bh_adjust(tests.loc[valid, "exact_permutation_p"])
    tests.to_csv(OUT / "gse293580_pde8b_patient_level_tests.csv", index=False)

    annotation_counts = (
        adata.obs.groupby(["locked_mural_annotation", "condition", "sample"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    annotation_counts.to_csv(OUT / "gse293580_locked_annotation_sample_counts.csv", index=False)
    adata.write_h5ad(INTERIM / "gse293580_phase3_mural_annotated.h5ad", compression="gzip")

    summary = {
        "integer_count_layer_max_rounding_error": integer_max_error,
        "locked_annotation": ANNOTATION,
        "strict_smc_cells": int(strict_smc.sum()),
        "pericyte_cells": int(pericyte.sum()),
        "broad_mural_cells": int(broad_mural.sum()),
        "minimum_cells_for_expression_test": 5,
        "statistical_unit": "patient/sample",
        "note": "PDE8B was not used for clustering or annotation; tests are exploratory owing to small patient counts.",
    }
    (OUT / "gse293580_phase3_pde8b_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
