import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT, OUTPUT_ROOT

import pandas as pd


DATA = DATA_ROOT / "pah_scrna"
OUT = OUTPUT_ROOT / "pah_audit"


def hedges_g(a, b):
    a, b = list(a), list(b)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    sp = math.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    if sp == 0:
        return 0.0
    d = (ma - mb) / sp
    return d * (1 - 3 / (4 * (len(a) + len(b)) - 9))


def main():
    base = DATA / "gse210248" / "Human_PA_Integrated"
    meta = pd.read_csv(
        DATA / "GSE210248_Human_PA_Metadata.csv.gz",
        sep=";", decimal=",", index_col=0,
    )
    prioritized = pd.read_csv(OUT / "lung_mechanism_candidate_prioritization.csv")
    candidates = set(prioritized.loc[prioritized["priority_pass"], "symbol"].astype(str))

    with gzip.open(base / "features.tsv.gz", "rt") as handle:
        features = [line.rstrip("\n").split("\t")[1] for line in handle]
    with gzip.open(base / "barcodes.tsv.gz", "rt") as handle:
        barcodes = [line.strip() for line in handle]
    if barcodes != list(meta.index):
        raise ValueError("Matrix barcodes and metadata rows are not aligned")

    wanted_rows = {i + 1: gene for i, gene in enumerate(features) if gene in candidates}
    groups = list(zip(meta["new.ident"].astype(str), meta["Cell_annotation"].astype(str)))
    group_cells = defaultdict(int)
    for group in groups:
        group_cells[group] += 1

    library = defaultdict(float)
    gene_counts = defaultdict(float)
    detected = defaultdict(int)
    matrix_path = base / "matrix.mtx.gz"
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3:
                nrow, ncol, nnz = map(int, parts)
                break
        for line in handle:
            row, col, value = line.split()
            row, col, value = int(row), int(col), float(value)
            group = groups[col - 1]
            library[group] += value
            if row in wanted_rows:
                gene = wanted_rows[row]
                gene_counts[(gene,) + group] += value
                detected[(gene,) + group] += 1

    records = []
    for gene in sorted(candidates):
        for (sample, celltype), ncells in sorted(group_cells.items()):
            count = gene_counts[(gene, sample, celltype)]
            lib = library[(sample, celltype)]
            records.append({
                "gene": gene,
                "sample": sample,
                "disease": "PAH" if sample.startswith("PAH") else "Donor",
                "cell_type": celltype,
                "n_cells": ncells,
                "count": count,
                "log_cp10k": math.log1p(10000 * count / lib) if lib else 0.0,
                "detection_fraction": detected[(gene, sample, celltype)] / ncells,
            })
    pb = pd.DataFrame(records)
    pb.to_csv(OUT / "scrna_gse210248_candidate_pseudobulk.csv", index=False)

    effects = []
    for (gene, celltype), block in pb.groupby(["gene", "cell_type"]):
        pah = block.loc[block.disease == "PAH", "log_cp10k"]
        donor = block.loc[block.disease == "Donor", "log_cp10k"]
        if len(pah) != 3 or len(donor) != 3:
            continue
        effects.append({
            "gene": gene,
            "cell_type": celltype,
            "pah_mean_log_cp10k": pah.mean(),
            "donor_mean_log_cp10k": donor.mean(),
            "delta_log_cp10k": pah.mean() - donor.mean(),
            "hedges_g": hedges_g(pah, donor),
            "pah_detection": block.loc[block.disease == "PAH", "detection_fraction"].mean(),
            "donor_detection": block.loc[block.disease == "Donor", "detection_fraction"].mean(),
        })
    effects = pd.DataFrame(effects)
    effects.to_csv(OUT / "scrna_gse210248_candidate_celltype_effects.csv", index=False)

    top = (effects.assign(abs_delta=lambda x: x.delta_log_cp10k.abs())
           .sort_values(["gene", "abs_delta"], ascending=[True, False])
           .groupby("gene").head(1))
    summary = {
        "cells": len(meta),
        "samples": meta["new.ident"].value_counts().to_dict(),
        "cell_types": meta["Cell_annotation"].value_counts().to_dict(),
        "candidate_genes_requested": len(candidates),
        "candidate_genes_present": len(wanted_rows),
        "top_localization": top[["gene", "cell_type", "delta_log_cp10k", "hedges_g"]].to_dict("records"),
    }
    (OUT / "scrna_gse210248_candidate_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "top_localization"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()



