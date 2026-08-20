import gzip
import math
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT, OUTPUT_ROOT

import pandas as pd


DATA = DATA_ROOT / "pah_scrna" / "gse293580"
OUT = OUTPUT_ROOT / "pah_audit"


GROUPS = {
    "GSM8885793": "Donor", "GSM8885794": "Donor", "GSM8885795": "Donor", "GSM8885796": "Donor",
    "GSM8885797": "SSc-PAH", "GSM8885798": "SSc-PAH", "GSM8885799": "SSc-PAH", "GSM8885800": "SSc-PAH",
    "GSM8885801": "IPAH", "GSM8885802": "IPAH", "GSM8885803": "IPAH", "GSM8885804": "IPAH",
}


def hedges_g(a, b):
    a, b = list(a), list(b)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    sp = math.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    if sp == 0:
        return 0.0
    return ((ma - mb) / sp) * (1 - 3 / (4 * (len(a) + len(b)) - 9))


def main():
    pri = pd.read_csv(OUT / "lung_mechanism_candidate_prioritization.csv")
    candidates = set(pri.loc[pri.priority_pass, "symbol"].astype(str))
    records = []
    for matrix in sorted(DATA.glob("GSM*_matrix.mtx.gz")):
        acc = matrix.name.split("_")[0]
        feature_file = Path(str(matrix).replace("matrix.mtx.gz", "features.tsv.gz"))
        barcode_file = Path(str(matrix).replace("matrix.mtx.gz", "barcodes.tsv.gz"))
        with gzip.open(feature_file, "rt") as f:
            features = [line.rstrip().split("\t")[1] for line in f]
        with gzip.open(barcode_file, "rt") as f:
            n_cells = sum(1 for _ in f)
        wanted = {i + 1: g for i, g in enumerate(features) if g in candidates}
        counts = defaultdict(float)
        library = 0.0
        with gzip.open(matrix, "rt") as f:
            for line in f:
                if line.startswith("%"):
                    continue
                nrow, ncol, nnz = map(int, line.split())
                break
            for line in f:
                row, col, val = line.split()
                row, val = int(row), float(val)
                library += val
                if row in wanted:
                    counts[wanted[row]] += val
        for gene in candidates:
            records.append({"accession": acc, "sample": matrix.name.split("_")[2], "condition": GROUPS[acc],
                            "n_cells": n_cells, "gene": gene,
                            "log_cp10k": math.log1p(10000 * counts[gene] / library)})
    values = pd.DataFrame(records)
    values.to_csv(OUT / "scrna_gse293580_candidate_pseudobulk.csv", index=False)

    meta = pri.set_index("symbol")["meta_g"].to_dict()
    effects = []
    for gene, block in values.groupby("gene"):
        donor = block.loc[block.condition.eq("Donor"), "log_cp10k"]
        for condition in ["IPAH", "SSc-PAH"]:
            case = block.loc[block.condition.eq(condition), "log_cp10k"]
            g = hedges_g(case, donor)
            effects.append({"gene": gene, "contrast": f"{condition}_vs_Donor", "hedges_g": g,
                            "delta_log_cp10k": case.mean() - donor.mean(), "bulk_meta_g": meta.get(gene),
                            "direction_concordant_with_bulk": (g > 0) == (meta.get(gene, 0) > 0)})
    effects = pd.DataFrame(effects)
    effects.to_csv(OUT / "scrna_gse293580_candidate_effects.csv", index=False)
    focus = ["PDE8B", "PIEZO2", "DCBLD1", "NCAM1", "SELP", "ABCG2"]
    print(values.groupby("condition").n_cells.first())
    print(effects[effects.gene.isin(focus)].sort_values(["gene", "contrast"]).to_string(index=False))


if __name__ == "__main__":
    main()



