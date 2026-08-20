#!/usr/bin/env python3
"""Audit GSE113439 PAH subtype mixing with an IPAH-only sensitivity analysis."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from project_paths import OUTPUT_ROOT
from formal_meta_and_confounding import hedges_g, prepare_gene_matrix
from lung_mechanism_meta import classify_array, read_geo_metadata


OUT = OUTPUT_ROOT / "public-data-audit"
CANDIDATES = ["PDE8B", "PIEZO2", "SLC16A12"]


def subtype(row: dict[str, object]) -> str:
    text = str(row["text"])
    if "ctep h" in text or "cteph" in text:
        return "CTEPH"
    if "idiopathic pah" in text:
        return "IPAH"
    if "pah and ctd" in text:
        return "CTD-PAH"
    if "pah and chd" in text:
        return "CHD-PAH"
    if "normal control" in text:
        return "Control"
    return "Unresolved"


def candidate_effects(matrix, cases, controls):
    effects = hedges_g(matrix[cases], matrix[controls])
    return {
        gene: float(effects.loc[gene]) if gene in effects.index else None
        for gene in CANDIDATES
    }


def pca_diagnostics(matrix, samples):
    x = matrix[samples].T.to_numpy(dtype=float, copy=True)
    x -= x.mean(axis=0, keepdims=True)
    scale = x.std(axis=0, ddof=1)
    keep = np.isfinite(scale) & (scale > 0)
    x = x[:, keep] / scale[keep]
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    variance = s**2 / np.sum(s**2)
    return {
        "variance_explained_pc1": float(variance[0]),
        "variance_explained_pc2": float(variance[1]),
        "scores": [
            {"sample": sample, "pc1": float(scores[i, 0]), "pc2": float(scores[i, 1])}
            for i, sample in enumerate(samples)
        ],
    }


def main():
    matrix, scale_qc = prepare_gene_matrix("GSE113439", "GPL6244")
    rows = read_geo_metadata("GSE113439")
    labels = {row["gsm"]: subtype(row) for row in rows}
    pipeline = {row["gsm"]: classify_array("GSE113439", row) for row in rows}
    controls = [gsm for gsm, label in labels.items() if label == "Control" and gsm in matrix.columns]
    ipah = [gsm for gsm, label in labels.items() if label == "IPAH" and gsm in matrix.columns]
    all_pah = [gsm for gsm, label in pipeline.items() if label == "PAH" and gsm in matrix.columns]
    cteph = [gsm for gsm, label in labels.items() if label == "CTEPH" and gsm in matrix.columns]
    selected = all_pah + controls

    result = {
        "dataset": "GSE113439",
        "purpose": "etiology-mixing audit; not a replacement primary analysis",
        "subtype_counts": {
            label: sum(value == label for value in labels.values())
            for label in sorted(set(labels.values()))
        },
        "pipeline_counts": {
            label: sum(value == label for value in pipeline.values())
            for label in sorted(set(pipeline.values()))
        },
        "cteph_excluded_samples": cteph,
        "all_pah_vs_control": {
            "n_case": len(all_pah),
            "n_control": len(controls),
            "candidate_hedges_g": candidate_effects(matrix, all_pah, controls),
        },
        "ipah_only_vs_control": {
            "n_case": len(ipah),
            "n_control": len(controls),
            "candidate_hedges_g": candidate_effects(matrix, ipah, controls),
        },
        "pca_all_pah_and_control": pca_diagnostics(matrix, selected),
        "scale_qc": scale_qc,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "GSE113439_etiology_sensitivity.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if not k.startswith("pca_")}, indent=2))
    print(path)


if __name__ == "__main__":
    main()
