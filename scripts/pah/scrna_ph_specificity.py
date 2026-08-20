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


def read_features(path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip().split("\t")[1] for line in handle]


def read_mtx_candidate_counts(path, features, candidates, column_groups=None):
    wanted = {i + 1: g for i, g in enumerate(features) if g in candidates}
    library = defaultdict(float)
    counts = defaultdict(float)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            nrow, ncol, nnz = map(int, line.split())
            break
        if column_groups is None:
            column_groups = ["sample"] * ncol
        for line in handle:
            row, col, value = line.split()
            row, col, value = int(row), int(col), float(value)
            group = column_groups[col - 1]
            library[group] += value
            if row in wanted:
                counts[(group, wanted[row])] += value
    return library, counts


def main():
    pri = pd.read_csv(OUT / "lung_mechanism_candidate_prioritization.csv")
    candidates = set(pri.loc[pri.priority_pass, "symbol"].astype(str))
    records = []

    # PAH versus donor, integrated pulmonary-artery matrix; aggregate by donor ID.
    base = DATA / "gse210248" / "Human_PA_Integrated"
    meta = pd.read_csv(DATA / "GSE210248_Human_PA_Metadata.csv.gz", sep=";", decimal=",", index_col=0)
    groups = meta["new.ident"].astype(str).tolist()
    lib, cnt = read_mtx_candidate_counts(base / "matrix.mtx.gz", read_features(base / "features.tsv.gz"), candidates, groups)
    for sample in sorted(set(groups)):
        for gene in candidates:
            records.append({"cohort": "GSE210248", "condition": "PAH" if sample.startswith("PAH") else "Donor",
                            "sample": sample, "gene": gene,
                            "log_cp10k": math.log1p(10000 * cnt[(sample, gene)] / lib[sample])})

    # PHPF-PH versus donor, one raw 10x matrix per subject.
    base = DATA / "gse228643"
    for matrix in sorted(base.glob("*_matrix.mtx.gz")):
        stem = matrix.name.replace("_matrix.mtx.gz", "")
        label = stem.split("_", 1)[1]
        feature_path = base / matrix.name.replace("matrix.mtx.gz", "features.tsv.gz")
        lib, cnt = read_mtx_candidate_counts(matrix, read_features(feature_path), candidates)
        for gene in candidates:
            records.append({"cohort": "GSE228643", "condition": "PHPF-PH" if label.startswith("PHPF") else "Donor",
                            "sample": label, "gene": gene,
                            "log_cp10k": math.log1p(10000 * cnt[("sample", gene)] / lib["sample"])})

    values = pd.DataFrame(records)
    values.to_csv(OUT / "scrna_pah_phpf_whole_sample_pseudobulk.csv", index=False)
    rows = []
    for cohort, disease in [("GSE210248", "PAH"), ("GSE228643", "PHPF-PH")]:
        x = values[values.cohort.eq(cohort)]
        for gene, block in x.groupby("gene"):
            case = block.loc[block.condition.eq(disease), "log_cp10k"]
            ctrl = block.loc[block.condition.eq("Donor"), "log_cp10k"]
            rows.append({"gene": gene, "cohort": cohort, "contrast": f"{disease}_vs_Donor",
                         "case_mean": case.mean(), "donor_mean": ctrl.mean(),
                         "delta_log_cp10k": case.mean() - ctrl.mean(),
                         "case_n": len(case), "donor_n": len(ctrl)})
    effects = pd.DataFrame(rows)
    wide = effects.pivot(index="gene", columns="contrast", values="delta_log_cp10k").reset_index()
    wide["pah_direction"] = wide["PAH_vs_Donor"].apply(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))
    wide["phpf_direction"] = wide["PHPF-PH_vs_Donor"].apply(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))
    wide["direction_specific_to_pah"] = wide.pah_direction.ne(wide.phpf_direction)
    wide.to_csv(OUT / "scrna_pah_vs_phpf_candidate_specificity.csv", index=False)
    focus = ["PDE8B", "PIEZO2", "DCBLD1", "NCAM1", "SELP", "ABCG2"]
    print(wide[wide.gene.isin(focus)].to_string(index=False))


if __name__ == "__main__":
    main()



