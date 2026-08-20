import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT, OUTPUT_ROOT


DATA = DATA_ROOT / "pah_geo"
OUT = OUTPUT_ROOT / "pah_audit"
OUT.mkdir(parents=True, exist_ok=True)


def clean(value):
    return value.strip().strip('"')


def parse_matrix(path):
    sample_meta = {}
    n_rows = 0
    sample_ids = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0][1:]
                sample_meta.setdefault(key, []).append([clean(x) for x in parts[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif in_table and line.startswith('"ID_REF"'):
                sample_ids = [clean(x) for x in line.rstrip("\n").split("\t")[1:]]
            elif in_table and line.startswith("!series_matrix_table_end"):
                break
            elif in_table and line.strip():
                n_rows += 1

    def first(key):
        values = sample_meta.get(key, [[]])
        return values[0]

    titles = first("Sample_title")
    accessions = first("Sample_geo_accession")
    sources = first("Sample_source_name_ch1")
    characteristics = sample_meta.get("Sample_characteristics_ch1", [])
    platforms = first("Sample_platform_id")
    if not sample_ids:
        sample_ids = accessions

    records = []
    for i, gsm in enumerate(accessions):
        char_values = [row[i] for row in characteristics if i < len(row)]
        text = " | ".join([titles[i], sources[i], *char_values]).lower()
        if "ssc-ph-ild" in text or "ssc-ild-ph" in text or ("interstitial" in text and "hypertension" in text):
            group = "SSc-PH-ILD"
        elif ("lssc_pah" in text or "ssc-pah" in text or "ssc_ph" in text
              or "ssc with pah" in text or "ssc-ph" in text
              or "scleroderma-associated pulmonary" in text):
            group = "SSc-PAH"
        elif ("lssc_nopah" in text or "ssc_nopah" in text or "ssc-no" in text
              or "ssc without pah" in text or "without pulmonary" in text):
            group = "SSc-noPH"
        elif "idiopathic pulmonary arterial hypertension" in text or "ipah" in text:
            group = "IPAH"
        elif "normal" in text or "control" in text or "healthy" in text:
            group = "Healthy"
        else:
            group = "Unclassified"

        title = titles[i]
        subject = re.sub(r"_(R|1yr)$", "", title, flags=re.IGNORECASE)
        subject = re.sub(r"^(LSSC_(?:NoPAH|PAH)|Normal_NoPAH)_", "", subject, flags=re.IGNORECASE)
        is_technical_repeat = bool(re.search(r"_R$", title, flags=re.IGNORECASE))
        is_followup = bool(re.search(r"_1yr$", title, flags=re.IGNORECASE))
        records.append({
            "gsm": gsm,
            "title": title,
            "source": sources[i],
            "platform": platforms[i] if i < len(platforms) else "",
            "group": group,
            "subject_key": subject,
            "technical_repeat": is_technical_repeat,
            "followup": is_followup,
            "characteristics": " | ".join(char_values),
        })
    return records, n_rows, len(sample_ids)


def parse_annotation(path):
    gene_ids = set()
    gene_symbols = set()
    header = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("#", "!", "^")) or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if header is None:
                header = row
                continue
            rec = dict(zip(header, row))
            for gene_id in re.split(r"\s*///\s*|\s*//\s*|[,;]\s*", rec.get("Gene ID", "")):
                if gene_id and gene_id != "---" and gene_id.isdigit():
                    gene_ids.add(gene_id)
            symbol_value = rec.get("Gene Symbol", rec.get("Gene symbol", ""))
            for symbol in re.split(r"\s*///\s*|\s*//\s*|[,;]\s*", symbol_value):
                if symbol and symbol != "---":
                    gene_symbols.add(symbol.upper())
    return gene_ids, gene_symbols


def main():
    datasets = {}
    all_records = []
    for path in sorted(DATA.glob("GSE*_series_matrix.txt.gz")):
        gse = path.name.split("_")[0]
        records, n_features, n_columns = parse_matrix(path)
        all_records.extend({"dataset": gse, **record} for record in records)
        baseline_candidates = [r for r in records if not r["followup"]]
        baseline_by_subject = {}
        for record in baseline_candidates:
            key = (record["group"], record["subject_key"].lower())
            current = baseline_by_subject.get(key)
            if current is None or (current["technical_repeat"] and not record["technical_repeat"]):
                baseline_by_subject[key] = record
        baseline = list(baseline_by_subject.values())
        datasets[gse] = {
            "matrix_samples": len(records),
            "matrix_columns": n_columns,
            "features": n_features,
            "platforms": sorted(set(r["platform"] for r in records)),
            "all_group_counts": dict(Counter(r["group"] for r in records)),
            "baseline_unique_counts": dict(Counter(r["group"] for r in baseline)),
            "technical_repeats": sum(r["technical_repeat"] for r in records),
            "followups": sum(r["followup"] for r in records),
            "unclassified_titles": [r["title"] for r in records if r["group"] == "Unclassified"],
        }

    annotations = {}
    for path in sorted(DATA.glob("GPL*.annot.gz")):
        gpl = path.name.split(".")[0]
        ids, symbols = parse_annotation(path)
        annotations[gpl] = {"gene_ids": len(ids), "gene_symbols": len(symbols), "symbol_set": symbols}

    symbol_sets = [value["symbol_set"] for value in annotations.values()]
    intersection = set.intersection(*symbol_sets) if symbol_sets else set()
    union = set.union(*symbol_sets) if symbol_sets else set()
    annotation_summary = {
        gpl: {k: v for k, v in value.items() if k != "symbol_set"}
        for gpl, value in annotations.items()
    }
    annotation_summary["all_platforms"] = {
        "common_gene_symbols": len(intersection),
        "union_gene_symbols": len(union),
    }

    with (OUT / "core_sample_manifest.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_records[0].keys()))
        writer.writeheader()
        writer.writerows(all_records)

    result = {"datasets": datasets, "annotations": annotation_summary}
    (OUT / "core_feasibility_audit.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()



