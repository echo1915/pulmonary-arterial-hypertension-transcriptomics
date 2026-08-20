"""Replicate GSE293580 ligand-receptor model SMC-side components in GSE210248.

GSE210248 has no independent pericyte annotation. This script therefore tests
only patient-level SMC ligand/receptor module directions and does not claim
cross-cohort replication of SMC-pericyte communication.
"""

from __future__ import annotations

import gzip
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import DATA_ROOT, PROJECT_DATA_ROOT  # noqa: E402


RAW = DATA_ROOT / "pah_scrna"
MATRIX_DIR = RAW / "gse210248" / "Human_PA_Integrated"
METADATA = RAW / "GSE210248_Human_PA_Metadata.csv.gz"
GSE293_COMPONENTS = (
    PROJECT_DATA_ROOT
    / "outputs"
    / "gse293580-reanalysis"
    / "ligand-receptor-proxy"
    / "gse293580_mural_lr_component_tests.csv"
)
OUT = PROJECT_DATA_ROOT / "outputs" / "gse210248-lr-smc-component-replication"

AXES = {
    "PDGF": {"ligand": ["PDGFB"], "receptor": ["PDGFRB"]},
    "TGFb": {"ligand": ["TGFB1", "TGFB2"], "receptor": ["TGFBR1", "TGFBR2"]},
    "JAG_NOTCH3": {"ligand": ["JAG1", "JAG2"], "receptor": ["NOTCH3"]},
    "EDN": {"ligand": ["EDN1"], "receptor": ["EDNRA", "EDNRB"]},
    "CXCL12_CXCR4": {"ligand": ["CXCL12"], "receptor": ["CXCR4"]},
    "ECM_integrin": {
        "ligand": ["COL1A1", "COL1A2", "COL3A1", "FN1"],
        "receptor": ["ITGA1", "ITGA2", "ITGA5", "ITGB1"],
    },
    "THBS1_CD47": {"ligand": ["THBS1"], "receptor": ["CD47"]},
}
CELL_TYPES = ("SMC 1", "SMC 2")


def exact_p(case: np.ndarray, donor: np.ndarray) -> float:
    pooled = np.r_[case, donor]
    observed = abs(case.mean() - donor.mean())
    values = []
    for idx in itertools.combinations(range(len(pooled)), len(case)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        values.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
    return float(np.mean(np.asarray(values) >= observed - 1e-12))


def hedges_g(case: np.ndarray, donor: np.ndarray) -> float:
    pooled_sd = np.sqrt(
        ((len(case) - 1) * case.var(ddof=1) + (len(donor) - 1) * donor.var(ddof=1))
        / (len(case) + len(donor) - 2)
    )
    if pooled_sd == 0:
        return 0.0 if case.mean() == donor.mean() else np.sign(case.mean() - donor.mean()) * np.inf
    return (case.mean() - donor.mean()) / pooled_sd * (1 - 3 / (4 * (len(case) + len(donor)) - 9))


def loo(case: np.ndarray, donor: np.ndarray) -> tuple[int, int, float, float]:
    sign = np.sign(case.mean() - donor.mean())
    deltas = [np.delete(case, i).mean() - donor.mean() for i in range(len(case))]
    deltas.extend(case.mean() - np.delete(donor, i).mean() for i in range(len(donor)))
    values = np.asarray(deltas)
    return int(np.sum(np.sign(values) == sign)), len(values), float(values.min()), float(values.max())


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return pd.Series(out, index=values.index)


def main() -> None:
    metadata = pd.read_csv(METADATA, sep=";", decimal=",", index_col=0)
    with gzip.open(MATRIX_DIR / "features.tsv.gz", "rt") as handle:
        features = [line.rstrip().split("\t")[1] for line in handle]
    with gzip.open(MATRIX_DIR / "barcodes.tsv.gz", "rt") as handle:
        barcodes = [line.strip() for line in handle]
    if barcodes != list(metadata.index):
        raise ValueError("GSE210248 barcode-metadata alignment failed")

    genes = sorted({gene for axis in AXES.values() for role in ("ligand", "receptor") for gene in axis[role]})
    row_to_gene = {i + 1: gene for i, gene in enumerate(features) if gene in genes}
    missing = sorted(set(genes).difference(row_to_gene.values()))
    if missing:
        raise ValueError(f"missing prespecified genes: {missing}")

    samples = metadata["new.ident"].astype(str).to_numpy()
    cell_types = metadata["Cell_annotation"].astype(str).to_numpy()
    groups = list(zip(samples, cell_types))
    library = defaultdict(float)
    counts = defaultdict(float)
    detected_cells: dict[tuple[str, str, str], set[int]] = defaultdict(set)

    with gzip.open(MATRIX_DIR / "matrix.mtx.gz", "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            break
        for line in handle:
            row, column, value = line.split()
            row, column, value = int(row), int(column), float(value)
            sample, cell_type = groups[column - 1]
            if cell_type not in CELL_TYPES:
                continue
            library[(sample, cell_type)] += value
            if row in row_to_gene:
                gene = row_to_gene[row]
                counts[(gene, sample, cell_type)] += value
                if value > 0:
                    detected_cells[(gene, sample, cell_type)].add(column)

    cell_coverage = (
        metadata.loc[metadata["Cell_annotation"].isin(CELL_TYPES)]
        .groupby(["new.ident", "disease", "Cell_annotation"], observed=True)
        .size().rename("n_cells").reset_index()
        .rename(columns={"new.ident": "sample", "Cell_annotation": "cell_type"})
    )
    records = []
    for sample in sorted(metadata["new.ident"].astype(str).unique()):
        disease = "PAH" if sample.startswith("PAH") else "Donor"
        for cell_type in CELL_TYPES:
            n_cells = int(cell_coverage.loc[
                (cell_coverage["sample"] == sample) & (cell_coverage["cell_type"] == cell_type),
                "n_cells",
            ].iloc[0])
            lib = library[(sample, cell_type)]
            for gene in genes:
                count = counts[(gene, sample, cell_type)]
                cpm = 1e6 * count / lib if lib else np.nan
                records.append({
                    "gene": gene,
                    "sample": sample,
                    "disease": disease,
                    "cell_type": cell_type,
                    "n_cells": n_cells,
                    "library_counts": lib,
                    "counts": count,
                    "detected_cells": len(detected_cells[(gene, sample, cell_type)]),
                    "log1p_CPM": np.log1p(cpm),
                })
    expression = pd.DataFrame(records)

    score_rows = []
    for axis, roles in AXES.items():
        for role, members in roles.items():
            block = expression.loc[expression["gene"].isin(members)]
            scores = (
                block.groupby(["sample", "disease", "cell_type"], as_index=False)
                .agg(module_log1p_CPM=("log1p_CPM", "mean"), genes_observed=("gene", "nunique"), n_cells=("n_cells", "first"))
            )
            scores["axis"] = axis
            scores["role"] = role
            score_rows.append(scores)
    scores = pd.concat(score_rows, ignore_index=True)

    effect_rows = []
    for (cell_type, axis, role), block in scores.groupby(["cell_type", "axis", "role"]):
        case = block.loc[block["disease"] == "PAH", "module_log1p_CPM"].to_numpy(float)
        donor = block.loc[block["disease"] == "Donor", "module_log1p_CPM"].to_numpy(float)
        concordant, iterations, loo_min, loo_max = loo(case, donor)
        effect_rows.append({
            "cell_type": cell_type,
            "axis": axis,
            "role": role,
            "n_PAH": len(case),
            "n_Donor": len(donor),
            "mean_difference": case.mean() - donor.mean(),
            "hedges_g": hedges_g(case, donor),
            "exact_permutation_p": exact_p(case, donor),
            "loo_direction_concordant_n": concordant,
            "loo_iterations_n": iterations,
            "loo_difference_min": loo_min,
            "loo_difference_max": loo_max,
        })
    effects = pd.DataFrame(effect_rows)
    effects["BH_q_within_cell_type_role"] = np.nan
    for _, idx in effects.groupby(["cell_type", "role"]).groups.items():
        effects.loc[idx, "BH_q_within_cell_type_role"] = bh_adjust(effects.loc[idx, "exact_permutation_p"])

    # Compare only the molecular component located in SMC in GSE293580:
    # ligand for SMC->pericyte and receptor for pericyte->SMC.
    gse293 = pd.read_csv(GSE293_COMPONENTS)
    gse293 = gse293.loc[
        (gse293["smc_boundary"] == "core_SMC")
        & (gse293["contrast"] == "PAH_combined_vs_Donor")
        & (
            ((gse293["direction"] == "SMC_to_pericyte") & (gse293["component"] == "sender_ligand"))
            | ((gse293["direction"] == "pericyte_to_SMC") & (gse293["component"] == "receiver_receptor"))
        )
    ].copy()
    gse293["role"] = np.where(gse293["direction"] == "SMC_to_pericyte", "ligand", "receptor")

    cross_rows = []
    for _, row293 in gse293.iterrows():
        for cell_type in CELL_TYPES:
            row210 = effects.loc[
                (effects["cell_type"] == cell_type)
                & (effects["axis"] == row293["axis"])
                & (effects["role"] == row293["role"])
            ].iloc[0]
            cross_rows.append({
                "axis": row293["axis"],
                "smc_role": row293["role"],
                "gse293580_direction_context": row293["direction"],
                "gse210248_cell_type": cell_type,
                "gse293580_hedges_g": row293["hedges_g"],
                "gse210248_hedges_g": row210["hedges_g"],
                "direction_concordant": bool(np.sign(row293["hedges_g"]) == np.sign(row210["hedges_g"])),
                "gse293580_exact_p": row293["exact_permutation_p"],
                "gse210248_exact_p": row210["exact_permutation_p"],
                "inference": "SMC_component_direction_only_not_communication_replication",
            })
    cross = pd.DataFrame(cross_rows)

    OUT.mkdir(parents=True, exist_ok=True)
    expression.to_csv(OUT / "gse210248_lr_gene_smc_patient_pseudobulk.csv", index=False)
    scores.to_csv(OUT / "gse210248_lr_smc_patient_module_scores.csv", index=False)
    effects.to_csv(OUT / "gse210248_lr_smc_component_effects.csv", index=False)
    cross.to_csv(OUT / "gse210248_gse293580_lr_smc_component_concordance.csv", index=False)
    summary = {
        "statistical_unit": "patient/sample",
        "gse210248_groups": "3 PAH versus 3 donors",
        "cell_types": list(CELL_TYPES),
        "pericyte_annotation_available": False,
        "scope": "SMC-side ligand/receptor component replication only",
        "axes": AXES,
        "overall_direction_concordance_fraction": float(cross["direction_concordant"].mean()),
        "inference_boundary": "No cross-cohort SMC-pericyte communication claim is possible without a pericyte annotation in GSE210248.",
    }
    (OUT / "gse210248_lr_smc_component_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(effects.sort_values(["cell_type", "role", "hedges_g"], ascending=[True, True, False]).to_string(index=False))
    print("\nCross-cohort SMC-component concordance:")
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
