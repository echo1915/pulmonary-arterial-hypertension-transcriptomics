"""Prepare audited GSE210248 SMC matrices for scTenifoldKnk feasibility tests."""

from __future__ import annotations

import csv
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite

from project_paths import PROJECT_DATA_ROOT


RAW = PROJECT_DATA_ROOT / "raw" / "current-data-snapshot" / "pah_scrna"
INTERIM = PROJECT_DATA_ROOT / "interim" / "sctenifoldknk-pde8b"
ARCHIVE = RAW / "GSE210248_Human_PA_Integrated.tar.gz"
METADATA = RAW / "GSE210248_Human_PA_Metadata.csv.gz"


def extract_inputs() -> Path:
    target = INTERIM / "Human_PA_Integrated"
    target.mkdir(parents=True, exist_ok=True)
    expected = [target / "barcodes.tsv.gz", target / "features.tsv.gz", target / "matrix.mtx.gz"]
    if not all(p.exists() for p in expected):
        with tarfile.open(ARCHIVE, "r:gz") as handle:
            members = [m for m in handle.getmembers() if Path(m.name).name in {p.name for p in expected}]
            handle.extractall(INTERIM, members=members, filter="data")
    return target


def read_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def main() -> None:
    matrix_dir = extract_inputs()
    barcodes = read_lines(matrix_dir / "barcodes.tsv.gz")
    feature_rows = [row.split("\t") for row in read_lines(matrix_dir / "features.tsv.gz")]
    genes = [row[1] if len(row) > 1 else row[0] for row in feature_rows]
    meta = pd.read_csv(METADATA, sep=";")
    matrix = mmread(matrix_dir / "matrix.mtx.gz").tocsr()

    barcode_candidates = [c for c in meta.columns if c.lower() in {"barcode", "cell", "cell_id", "cellid", "unnamed: 0"}]
    if not barcode_candidates:
        overlap = {c: len(set(meta[c].astype(str)).intersection(barcodes)) for c in meta.columns}
        barcode_col = max(overlap, key=overlap.get)
    else:
        barcode_col = barcode_candidates[0]
    meta[barcode_col] = meta[barcode_col].astype(str)
    meta = meta.set_index(barcode_col).reindex(barcodes)

    cell_type_col = "Cell_annotation"
    sample_col = "new.ident"
    disease_col = next(c for c in meta.columns if c.lower() in {"disease", "condition", "group"})
    if matrix.shape == (len(barcodes), len(genes)):
        matrix = matrix.T.tocsr()
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(f"Matrix shape {matrix.shape} does not match genes x cells {(len(genes), len(barcodes))}")

    INTERIM.mkdir(parents=True, exist_ok=True)
    pde_rows = np.flatnonzero(np.asarray(genes) == "PDE8B")
    if len(pde_rows) != 1:
        raise ValueError(f"Expected exactly one PDE8B feature, found {len(pde_rows)}")

    audit_rows = []
    for subset_name, allowed in {
        "SMC_all": {"SMC 1", "SMC 2"},
        "SMC1": {"SMC 1"},
        "SMC2": {"SMC 2"},
        "SMC_all_PAH": {"SMC 1", "SMC 2"},
    }.items():
        mask = meta[cell_type_col].isin(allowed).fillna(False).to_numpy(copy=True)
        if subset_name.endswith("_PAH"):
            mask &= meta[disease_col].astype(str).eq("PAH").fillna(False).to_numpy()
        idx = np.flatnonzero(mask)
        sub = matrix[:, idx].tocsc()
        pde = np.asarray(sub[pde_rows[0], :].todense()).ravel()
        libs = np.asarray(sub.sum(axis=0)).ravel()
        detected_genes = np.asarray((sub > 0).sum(axis=0)).ravel()
        audit_rows.append({
            "subset": subset_name,
            "cells": len(idx),
            "subjects": int(meta.iloc[idx][sample_col].nunique()),
            "pde8b_nonzero_cells": int((pde > 0).sum()),
            "pde8b_detected_fraction": float((pde > 0).mean()) if len(idx) else 0.0,
            "median_library_size": float(np.median(libs)) if len(idx) else 0.0,
            "median_detected_genes": float(np.median(detected_genes)) if len(idx) else 0.0,
        })
        out = INTERIM / subset_name
        out.mkdir(exist_ok=True)
        mmwrite(out / "counts.mtx", sub)
        pd.Series(genes).to_csv(out / "genes.tsv", index=False, header=False)
        meta.iloc[idx].reset_index().to_csv(out / "cells.csv", index=False)

    pd.DataFrame(audit_rows).to_csv(INTERIM / "input_feasibility_audit.csv", index=False)
    (INTERIM / "input_manifest.json").write_text(json.dumps({
        "archive": str(ARCHIVE), "metadata": str(METADATA),
        "matrix_shape_genes_cells": list(matrix.shape),
        "barcode_column": barcode_col, "cell_type_column": cell_type_col,
        "sample_column": sample_col, "disease_column": disease_col,
    }, indent=2), encoding="utf-8")
    print(pd.DataFrame(audit_rows).to_string(index=False))


if __name__ == "__main__":
    main()
