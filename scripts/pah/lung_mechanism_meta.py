import gzip
import json
import math
import re
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT

import numpy as np
import pandas as pd

from formal_meta_and_confounding import OUT, bh_adjust, hedges_g, prepare_gene_matrix, random_effects


DATA = DATA_ROOT / "pah_geo"


def read_geo_metadata(gse):
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
    rows = []
    for i, gsm in enumerate(accessions):
        values = [row[i] for row in chars]
        text = " | ".join([titles[i], sources[i], *values]).lower()
        rows.append({"gsm": gsm, "title": titles[i], "text": text, "characteristics": values})
    return rows


def classify_array(gse, row):
    text = row["text"]
    if gse == "GSE15197":
        if "ipf" in text:
            return "IPF-PH"
        if "pulmonary arterial hypertension" in text or re.search(r"\bp[ah]{2},", text):
            return "PAH"
        if "control" in text or "normal lung" in text:
            return "Control"
    else:
        if "cteph" in text or "who 4" in text:
            return "Exclude"
        if "failed donor" in text or "lung_fd" in text or "normal control" in text or "control_lung" in text:
            return "Control"
        if "pah" in text:
            return "PAH"
    return "Exclude"


def prepare_array_cohort(gse):
    gpl = "GPL6480" if gse == "GSE15197" else "GPL6244"
    matrix, qc = prepare_gene_matrix(gse, gpl)
    rows = read_geo_metadata(gse)
    groups = {row["gsm"]: classify_array(gse, row) for row in rows}
    case = [gsm for gsm, group in groups.items() if group == "PAH" and gsm in matrix.columns]
    control = [gsm for gsm, group in groups.items() if group == "Control" and gsm in matrix.columns]
    other = {group: sum(v == group for v in groups.values()) for group in set(groups.values()) - {"PAH", "Control"}}
    return matrix[case + control], case, control, {**qc, "other_groups": other}


def prepare_gse254617():
    counts = pd.read_csv(DATA / "GSE254617_txiCounts.csv.gz", index_col=0)
    counts.index = counts.index.astype(str).str.upper()
    counts = counts.groupby(level=0).sum()
    subject_columns = {}
    for column in counts.columns:
        match = re.match(r"(PHBI_\d+)", column)
        if match:
            subject_columns.setdefault(match.group(1), []).append(column)
    subject_counts = pd.DataFrame({subject: counts[cols].sum(axis=1) for subject, cols in subject_columns.items()})
    cpm = subject_counts.div(subject_counts.sum(axis=0), axis=1) * 1_000_000
    matrix = np.log2(cpm + 1)

    rows = read_geo_metadata("GSE254617")
    groups = {}
    subtypes = {}
    for row in rows:
        condition = next((x.split(":", 1)[1].strip() for x in row["characteristics"] if x.lower().startswith("condition:")), "")
        subtype = next((x.split(":", 1)[1].strip() for x in row["characteristics"] if x.lower().startswith("ph subtype:")), "")
        subject = row["title"]
        subtypes[subject] = subtype
        if condition == "FD":
            groups[subject] = "Control"
        elif subtype in {"IPAH", "APAH", "HPAH"}:
            groups[subject] = "PAH"
        else:
            groups[subject] = f"Exclude-{subtype or 'Unknown'}"
    case = [x for x, group in groups.items() if group == "PAH" and x in matrix.columns]
    control = [x for x, group in groups.items() if group == "Control" and x in matrix.columns]
    other = {group: sum(v == group for v in groups.values()) for group in set(groups.values()) - {"PAH", "Control"}}
    qc = {
        "raw_count_columns": int(counts.shape[1]),
        "unique_subject_columns": int(matrix.shape[1]),
        "genes": int(matrix.shape[0]),
        "other_groups": other,
    }
    return matrix[case + control], case, control, qc


def ensembl_to_symbol():
    gene_to_symbol = {}
    with gzip.open(DATA / "Homo_sapiens.gene_info.gz", "rt", encoding="utf-8", errors="replace") as handle:
        header = next(handle).lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            gene_to_symbol[row["GeneID"]] = row["Symbol"].upper()
    mapping = {}
    with gzip.open(DATA / "gene2ensembl.gz", "rt", encoding="utf-8", errors="replace") as handle:
        header = next(handle).lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "9606":
                continue
            row = dict(zip(header, parts))
            ensembl = row.get("Ensembl_gene_identifier", "").split(".")[0]
            symbol = gene_to_symbol.get(row.get("GeneID", ""))
            if ensembl and ensembl != "-" and symbol:
                mapping[ensembl] = symbol
    return mapping


def prepare_gse208592():
    counts = pd.read_csv(DATA / "GSE208592_Counts_Lung.csv.gz", index_col=0)
    counts.index = counts.index.astype(str).str.split(".").str[0]
    mapping = ensembl_to_symbol()
    symbols = pd.Series(counts.index.map(mapping), index=counts.index)
    counts["symbol"] = symbols.values
    counts = counts.dropna(subset=["symbol"]).groupby("symbol").sum(numeric_only=True)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1_000_000
    matrix = np.log2(cpm + 1)
    case = [column for column in matrix.columns if column.startswith("PAH_Lung_")]
    control = [column for column in matrix.columns if column.startswith("Control_Lung_")]
    qc = {"mapped_genes": int(len(matrix)), "raw_samples": int(len(case) + len(control))}
    return matrix[case + control], case, control, qc


def cohort_effect(matrix, case, control):
    g = hedges_g(matrix[case], matrix[control])
    n1, n0 = len(case), len(control)
    variance = (n1 + n0) / (n1 * n0) + (g**2) / (2 * (n1 + n0 - 2))
    return g, variance


def main():
    cohorts = {}
    qc = {}
    for gse in ["GSE15197", "GSE113439"]:
        matrix, case, control, detail = prepare_array_cohort(gse)
        g, variance = cohort_effect(matrix, case, control)
        cohorts[gse] = {"g": g, "variance": variance}
        qc[gse] = {"n_pah": len(case), "n_control": len(control), "genes": len(matrix), **detail}
    matrix, case, control, detail = prepare_gse254617()
    g, variance = cohort_effect(matrix, case, control)
    cohorts["GSE254617"] = {"g": g, "variance": variance}
    qc["GSE254617"] = {"n_pah": len(case), "n_control": len(control), **detail}
    matrix, case, control, detail = prepare_gse208592()
    g, variance = cohort_effect(matrix, case, control)
    cohorts["GSE208592"] = {"g": g, "variance": variance}
    qc["GSE208592"] = {"n_pah": len(case), "n_control": len(control), **detail}

    common = set.intersection(*(set(value["g"].dropna().index) for value in cohorts.values()))
    rows = []
    for gene in sorted(common):
        gs = np.array([cohorts[gse]["g"].loc[gene] for gse in cohorts], dtype=float)
        vs = np.array([cohorts[gse]["variance"].loc[gene] for gse in cohorts], dtype=float)
        pooled, se, lo, hi, p, tau2, q, i2 = random_effects(gs, vs)
        rows.append({
            "symbol": gene,
            **{f"g_{gse}": gs[i] for i, gse in enumerate(cohorts)},
            "n_effects": int(np.isfinite(gs).sum()),
            "same_direction": bool(np.isfinite(gs).all() and len(set(np.sign(gs))) == 1),
            "meta_g": pooled, "meta_se": se, "ci_low": lo, "ci_high": hi,
            "p": p, "tau2": tau2, "Q": q, "I2": i2,
        })
    meta = pd.DataFrame(rows)
    meta["fdr"] = bh_adjust(meta["p"])
    meta = meta.sort_values(["fdr", "I2", "meta_g"], ascending=[True, True, False])
    meta.to_csv(OUT / "lung_mechanism_random_effects_meta.csv", index=False, encoding="utf-8-sig")
    robust = meta[(meta["same_direction"]) & (meta["fdr"] < 0.05) & (meta["I2"] < 50) & (meta["meta_g"].abs() >= 0.50)]

    validation_matrix, validation_case, validation_control, validation_qc = prepare_array_cohort("GSE53408")
    validation_g, _ = cohort_effect(validation_matrix, validation_case, validation_control)
    validation = []
    for _, row in robust.iterrows():
        gene = row["symbol"]
        observed = float(validation_g.loc[gene]) if gene in validation_g.index and np.isfinite(validation_g.loc[gene]) else None
        validation.append({
            "symbol": gene,
            "meta_g": float(row["meta_g"]),
            "present_in_GSE53408_filtered_matrix": observed is not None,
            "g_GSE53408": observed,
            "direction_concordant": bool(observed is not None and np.sign(observed) == np.sign(row["meta_g"])),
        })
    summary = {
        "cohort_qc": qc,
        "independence_note": "GSE117261 excluded because its PHBI subjects overlap the expanded GSE254617 RNA-seq cohort. GSE53408 excluded from primary meta-analysis because its public matrix contains only ~2,882 prefiltered genes.",
        "common_genes": int(len(meta)),
        "fdr_lt_0.05": int((meta["fdr"] < 0.05).sum()),
        "robust_definition": "same direction in all 4 cohorts, FDR<0.05, I2<50%, abs(meta_g)>=0.50",
        "robust_count": int(len(robust)),
        "top_robust": robust.head(50).to_dict("records"),
        "secondary_validation_GSE53408": {"qc": validation_qc, "candidates": validation},
    }
    (OUT / "lung_mechanism_meta_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()



