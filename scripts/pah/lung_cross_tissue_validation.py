import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from formal_meta_and_confounding import MODULES, OUT, auc_rank, prepare_gene_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT

DATA = DATA_ROOT / "pah_geo"


DATASETS = ["GSE53408", "GSE113439", "GSE117261"]


def sample_metadata(gse):
    path = DATA / f"{gse}_series_matrix.txt.gz"
    metadata = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                metadata.setdefault(parts[0][1:], []).append(parts[1:])
            if line.startswith("!series_matrix_table_begin"):
                break
    accessions = metadata["Sample_geo_accession"][0]
    titles = metadata["Sample_title"][0]
    sources = metadata["Sample_source_name_ch1"][0]
    chars = metadata.get("Sample_characteristics_ch1", [])
    records = []
    for i, gsm in enumerate(accessions):
        text = " | ".join([titles[i], sources[i], *[row[i] for row in chars]]).lower()
        if "failed donor" in text or "lung_fd" in text or "normal control" in text or "control_lung" in text:
            group = "Control"
        elif "cteph" in text or "who 4" in text:
            group = "Exclude-CTEPH"
        elif "pah" in text:
            group = "PAH"
        else:
            group = "Unclassified"
        records.append((gsm, group, titles[i], text))
    return records


def hedges(case, control):
    n1, n0 = len(case), len(control)
    pooled = np.sqrt(((n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / (n1 + n0 - 2))
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(((np.mean(case) - np.mean(control)) / pooled) * correction) if pooled > 0 else None


def main():
    result = {}
    for gse in DATASETS:
        records = sample_metadata(gse)
        groups = pd.Series({gsm: group for gsm, group, _, _ in records})
        case = groups[groups == "PAH"].index.tolist()
        control = groups[groups == "Control"].index.tolist()
        matrix, qc = prepare_gene_matrix(gse, "GPL6244")
        case = [x for x in case if x in matrix.columns]
        control = [x for x in control if x in matrix.columns]
        selected = matrix[case + control]
        y = np.array([1] * len(case) + [0] * len(control))
        z = selected.sub(selected.mean(axis=1), axis=0).div(selected.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        features = {}
        for gene in ["NRCAM", "HBQ1", "HBM", "HBB"]:
            if gene in selected.index:
                values = selected.loc[gene].to_numpy(dtype=float)
                features[gene] = {"present": True, "g": hedges(values[:len(case)], values[len(case):]), "auc_case_high": auc_rank(y, values)}
            else:
                features[gene] = {"present": False}
        for name, genes in MODULES.items():
            present = [gene for gene in genes if gene in z.index]
            if present:
                values = z.loc[present].mean(axis=0).to_numpy(dtype=float)
                features[f"module_{name}"] = {"genes_used": present, "g": hedges(values[:len(case)], values[len(case):]), "auc_case_high": auc_rank(y, values)}
        result[gse] = {
            "n_pah": len(case),
            "n_control": len(control),
            "excluded_or_unclassified": int(len(records) - len(case) - len(control)),
            "matrix_genes": int(len(matrix)),
            "scale_qc": qc,
            "features": features,
        }
    (OUT / "lung_cross_tissue_validation.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()


