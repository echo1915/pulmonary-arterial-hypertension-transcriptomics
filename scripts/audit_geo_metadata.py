#!/usr/bin/env python3
"""Parse GEO series-matrix sample metadata without loading expression values."""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

from project_paths import DATA_ROOT, OUTPUT_ROOT


INPUT = DATA_ROOT / "geo_metadata"
OUTPUT = OUTPUT_ROOT / "audit"


def clean(value: str) -> str:
    return value.strip().strip('"')


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def parse(path: Path):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                rows.append(next(csv.reader([line.rstrip("\n")], delimiter="\t")))

    accessions = next(r[1:] for r in rows if r[0] == "!Sample_geo_accession")
    samples = [{"geo_accession": clean(x)} for x in accessions]
    seen = Counter()
    for row in rows:
        base = row[0].removeprefix("!Sample_")
        if base == "geo_accession":
            continue
        seen[base] += 1
        column = slug(base if seen[base] == 1 else f"{base}_{seen[base]}")
        for sample, value in zip(samples, row[1:]):
            sample[column] = clean(value)

    # Promote every key:value characteristic/description entry to a named column.
    for sample in samples:
        for key, value in list(sample.items()):
            if key.startswith("characteristics_ch1") or key.startswith("description"):
                if ":" in value:
                    label, content = value.split(":", 1)
                    sample[slug(label)] = content.strip()
    return samples


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report = {}
    for path in sorted(INPUT.glob("GSE*_series_matrix.txt.gz")):
        gse = path.name.split("_")[0]
        samples = parse(path)
        columns = sorted({k for s in samples for k in s})
        with (OUTPUT / f"{gse}_sample_metadata.csv").open("w", newline="", encoding="utf-8-sig") as out:
            writer = csv.DictWriter(out, fieldnames=columns)
            writer.writeheader()
            writer.writerows(samples)

        categorical = {}
        for column in columns:
            vals = [s.get(column, "") for s in samples]
            counts = Counter(v for v in vals if v != "")
            if 0 < len(counts) <= 20:
                categorical[column] = dict(counts.most_common())
        report[gse] = {
            "n_samples": len(samples),
            "columns": columns,
            "categorical_counts": categorical,
        }

    (OUTPUT / "metadata_audit.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("Audited: " + ", ".join(f"{gse}={item['n_samples']} samples" for gse, item in report.items()))


if __name__ == "__main__":
    main()



