"""Patient-level SMC-pericyte ligand-receptor expression compatibility audit.

The resulting scores are expression-based communication proxies, not spatial
evidence, causal effects, or CellChat/NicheNet probabilities.
"""

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
OUT = PROJECT_DATA_ROOT / "outputs" / "gse293580-reanalysis" / "ligand-receptor-proxy"
MIN_CELLS = 5

SMC_BOUNDARIES = {
    "core_SMC": {"contractile_SMC", "modulated_SMC"},
    "strict_SMC_with_like": {"contractile_SMC", "modulated_SMC", "SMC_like"},
}
AXES = {
    "PDGF": {"ligands": ["PDGFB"], "receptors": ["PDGFRB"]},
    "TGFb": {"ligands": ["TGFB1", "TGFB2"], "receptors": ["TGFBR1", "TGFBR2"]},
    "JAG_NOTCH3": {"ligands": ["JAG1", "JAG2"], "receptors": ["NOTCH3"]},
    "EDN": {"ligands": ["EDN1"], "receptors": ["EDNRA", "EDNRB"]},
    "CXCL12_CXCR4": {"ligands": ["CXCL12"], "receptors": ["CXCR4"]},
    "ECM_integrin": {
        "ligands": ["COL1A1", "COL1A2", "COL3A1", "FN1"],
        "receptors": ["ITGA1", "ITGA2", "ITGA5", "ITGB1"],
    },
    "THBS1_CD47": {"ligands": ["THBS1"], "receptors": ["CD47"]},
}
CONTRASTS = {
    "IPAH_vs_Donor": ("IPAH",),
    "SSc-PAH_vs_Donor": ("SSc-PAH",),
    "PAH_combined_vs_Donor": ("IPAH", "SSc-PAH"),
}


def exact_p(case: np.ndarray, donor: np.ndarray) -> float:
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
    return (case.mean() - donor.mean()) / pooled_sd * (1 - 3 / (4 * (len(case) + len(donor)) - 9))


def loo(case: np.ndarray, donor: np.ndarray) -> tuple[int, int, float, float]:
    case = case[np.isfinite(case)]
    donor = donor[np.isfinite(donor)]
    if not len(case) or not len(donor):
        return 0, 0, np.nan, np.nan
    full_sign = np.sign(case.mean() - donor.mean())
    deltas = []
    if len(case) > 1:
        deltas.extend(np.delete(case, i).mean() - donor.mean() for i in range(len(case)))
    if len(donor) > 1:
        deltas.extend(case.mean() - np.delete(donor, i).mean() for i in range(len(donor)))
    values = np.asarray(deltas, dtype=float)
    return (
        int(np.sum(np.sign(values) == full_sign)),
        len(values),
        float(values.min()),
        float(values.max()),
    )


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return pd.Series(out, index=values.index)


def main() -> None:
    import anndata as ad

    genes = sorted({gene for axis in AXES.values() for key in ("ligands", "receptors") for gene in axis[key]})
    adata = ad.read_h5ad(INPUT, backed="r")
    missing = sorted(set(genes).difference(adata.var_names))
    if missing:
        raise ValueError(f"missing prespecified ligand-receptor genes: {missing}")
    obs = adata.obs[["sample", "condition", "locked_mural_annotation"]].copy()
    gene_idx = [adata.var_names.get_loc(gene) for gene in genes]
    counts = adata.layers["counts"][:]
    adata.file.close()
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    library_per_cell = np.asarray(counts.sum(axis=1)).ravel()
    gene_counts = counts[:, gene_idx].toarray()

    labels = obs["locked_mural_annotation"].astype(str).to_numpy()
    sample_values = obs["sample"].to_numpy()
    compartments = {"pericyte": {"pericyte"}, **SMC_BOUNDARIES}
    samples = obs[["sample", "condition"]].drop_duplicates()
    expression_rows = []
    for compartment, included_labels in compartments.items():
        compartment_mask = np.isin(labels, list(included_labels))
        for sample, condition in samples.itertuples(index=False):
            take = compartment_mask & (sample_values == sample)
            n_cells = int(take.sum())
            library = float(library_per_cell[take].sum())
            sums = gene_counts[take].sum(axis=0)
            detected = (gene_counts[take] > 0).sum(axis=0)
            for gene, count, detected_cells in zip(genes, sums, detected):
                cpm = 1e6 * float(count) / library if library > 0 else np.nan
                expression_rows.append({
                    "compartment": compartment,
                    "sample": sample,
                    "condition": condition,
                    "n_cells": n_cells,
                    "library_counts": library,
                    "gene": gene,
                    "counts": float(count),
                    "detected_cells": int(detected_cells),
                    "detection_fraction": float(detected_cells) / n_cells if n_cells else np.nan,
                    "log1p_CPM": np.log1p(cpm) if np.isfinite(cpm) else np.nan,
                })
    expression = pd.DataFrame(expression_rows)

    score_rows = []
    for boundary in SMC_BOUNDARIES:
        for sample, condition in samples.itertuples(index=False):
            for direction, sender, receiver in (
                ("SMC_to_pericyte", boundary, "pericyte"),
                ("pericyte_to_SMC", "pericyte", boundary),
            ):
                sender_block = expression.loc[
                    (expression["compartment"] == sender) & (expression["sample"] == sample)
                ].set_index("gene")
                receiver_block = expression.loc[
                    (expression["compartment"] == receiver) & (expression["sample"] == sample)
                ].set_index("gene")
                sender_cells = int(sender_block["n_cells"].iloc[0])
                receiver_cells = int(receiver_block["n_cells"].iloc[0])
                for axis, members in AXES.items():
                    ligand_score = float(sender_block.loc[members["ligands"], "log1p_CPM"].mean())
                    receptor_score = float(receiver_block.loc[members["receptors"], "log1p_CPM"].mean())
                    proxy = np.sqrt(ligand_score * receptor_score) if np.isfinite(ligand_score * receptor_score) else np.nan
                    score_rows.append({
                        "smc_boundary": boundary,
                        "direction": direction,
                        "axis": axis,
                        "sample": sample,
                        "condition": condition,
                        "sender_compartment": sender,
                        "receiver_compartment": receiver,
                        "sender_cells": sender_cells,
                        "receiver_cells": receiver_cells,
                        "eligible_pair": sender_cells >= MIN_CELLS and receiver_cells >= MIN_CELLS,
                        "ligand_module_log1p_CPM": ligand_score,
                        "receptor_module_log1p_CPM": receptor_score,
                        "compatibility_proxy": proxy,
                    })
    scores = pd.DataFrame(score_rows)

    test_rows = []
    eligible = scores.loc[scores["eligible_pair"]]
    for boundary in SMC_BOUNDARIES:
        for direction in ("SMC_to_pericyte", "pericyte_to_SMC"):
            for axis in AXES:
                block = eligible.loc[
                    (eligible["smc_boundary"] == boundary)
                    & (eligible["direction"] == direction)
                    & (eligible["axis"] == axis)
                ]
                donor = block.loc[block["condition"] == "Donor", "compatibility_proxy"].to_numpy(float)
                for contrast, conditions in CONTRASTS.items():
                    case = block.loc[
                        block["condition"].isin(conditions), "compatibility_proxy"
                    ].to_numpy(float)
                    concordant, iterations, loo_min, loo_max = loo(case, donor)
                    test_rows.append({
                        "smc_boundary": boundary,
                        "direction": direction,
                        "axis": axis,
                        "contrast": contrast,
                        "minimum_cells_each_compartment": MIN_CELLS,
                        "n_case": int(np.isfinite(case).sum()),
                        "n_donor": int(np.isfinite(donor).sum()),
                        "case_mean": np.nanmean(case),
                        "donor_mean": np.nanmean(donor),
                        "mean_difference": np.nanmean(case) - np.nanmean(donor),
                        "hedges_g": hedges_g(case, donor),
                        "exact_permutation_p": exact_p(case, donor),
                        "loo_direction_concordant_n": concordant,
                        "loo_iterations_n": iterations,
                        "loo_difference_min": loo_min,
                        "loo_difference_max": loo_max,
                    })
    tests = pd.DataFrame(test_rows)
    tests["BH_q_within_boundary_direction_contrast"] = np.nan
    families = ["smc_boundary", "direction", "contrast"]
    for _, idx in tests.groupby(families).groups.items():
        valid_idx = tests.loc[idx].index[tests.loc[idx, "exact_permutation_p"].notna()]
        if len(valid_idx):
            tests.loc[valid_idx, "BH_q_within_boundary_direction_contrast"] = bh_adjust(
                tests.loc[valid_idx, "exact_permutation_p"]
            )

    component_rows = []
    for boundary in SMC_BOUNDARIES:
        for direction in ("SMC_to_pericyte", "pericyte_to_SMC"):
            for axis in AXES:
                block = eligible.loc[
                    (eligible["smc_boundary"] == boundary)
                    & (eligible["direction"] == direction)
                    & (eligible["axis"] == axis)
                ]
                for contrast, conditions in CONTRASTS.items():
                    for component, column in (
                        ("sender_ligand", "ligand_module_log1p_CPM"),
                        ("receiver_receptor", "receptor_module_log1p_CPM"),
                    ):
                        case = block.loc[block["condition"].isin(conditions), column].to_numpy(float)
                        donor = block.loc[block["condition"] == "Donor", column].to_numpy(float)
                        concordant, iterations, loo_min, loo_max = loo(case, donor)
                        component_rows.append({
                            "smc_boundary": boundary,
                            "direction": direction,
                            "axis": axis,
                            "contrast": contrast,
                            "component": component,
                            "n_case": int(np.isfinite(case).sum()),
                            "n_donor": int(np.isfinite(donor).sum()),
                            "mean_difference": np.nanmean(case) - np.nanmean(donor),
                            "hedges_g": hedges_g(case, donor),
                            "exact_permutation_p": exact_p(case, donor),
                            "loo_direction_concordant_n": concordant,
                            "loo_iterations_n": iterations,
                            "loo_difference_min": loo_min,
                            "loo_difference_max": loo_max,
                        })
    components = pd.DataFrame(component_rows)
    components["BH_q_within_boundary_direction_contrast_component"] = np.nan
    component_families = ["smc_boundary", "direction", "contrast", "component"]
    for _, idx in components.groupby(component_families).groups.items():
        valid_idx = components.loc[idx].index[components.loc[idx, "exact_permutation_p"].notna()]
        if len(valid_idx):
            components.loc[valid_idx, "BH_q_within_boundary_direction_contrast_component"] = bh_adjust(
                components.loc[valid_idx, "exact_permutation_p"]
            )

    OUT.mkdir(parents=True, exist_ok=True)
    expression.to_csv(OUT / "gse293580_mural_lr_gene_patient_pseudobulk.csv", index=False)
    scores.to_csv(OUT / "gse293580_mural_lr_patient_proxy_scores.csv", index=False)
    tests.to_csv(OUT / "gse293580_mural_lr_patient_tests.csv", index=False)
    components.to_csv(OUT / "gse293580_mural_lr_component_tests.csv", index=False)

    primary = tests.loc[
        (tests["smc_boundary"] == "core_SMC")
        & (tests["contrast"] == "PAH_combined_vs_Donor")
    ].sort_values(["direction", "hedges_g"], ascending=[True, False])
    summary = {
        "statistical_unit": "patient/sample",
        "minimum_cells_each_compartment": MIN_CELLS,
        "axes": AXES,
        "smc_boundaries": {key: sorted(value) for key, value in SMC_BOUNDARIES.items()},
        "score_definition": "geometric mean of sender ligand-module and receiver receptor-module log1p CPM",
        "primary_scope": "core SMC, combined PAH versus donor",
        "inference_boundary": "Expression compatibility proxy only; no spatial, protein, or causal communication inference.",
    }
    (OUT / "gse293580_mural_lr_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(primary[[
        "direction", "axis", "n_case", "n_donor", "mean_difference", "hedges_g",
        "exact_permutation_p", "BH_q_within_boundary_direction_contrast",
        "loo_direction_concordant_n", "loo_iterations_n",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
