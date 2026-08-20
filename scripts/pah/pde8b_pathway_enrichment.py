import csv
import json
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(r"G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\outputs\pde8b-coexpression")
META = ROOT / "pde8b_cross_cohort_coexpression_meta.csv"
OUT = ROOT / "enrichment"
URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "REAC", "KEGG", "WP"]


def query_gprofiler(label, genes, background):
    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background,
        "no_evidences": False,
    }
    request = urllib.request.Request(
        URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "PAH-PDE8B-reproducible-analysis/1.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        result = json.loads(response.read().decode("utf-8"))
    (OUT / f"{label}_gprofiler_raw.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rows = result.get("result", [])
    query_genes = result.get("meta", {}).get("query_metadata", {}).get("queries", {}).get("query_1", genes)
    normalized = []
    for row in rows:
        evidence = row.get("intersections", [])
        hits = [query_genes[i] for i, sources in enumerate(evidence) if sources]
        normalized.append({
            "query_set": label,
            "source": row.get("source"),
            "term_id": row.get("native"),
            "term_name": row.get("name"),
            "fdr": row.get("p_value"),
            "term_size": row.get("term_size"),
            "query_size": row.get("query_size"),
            "intersection_size": row.get("intersection_size"),
            "precision": row.get("precision"),
            "recall": row.get("recall"),
            "intersection": ";".join(hits),
        })
    return normalized


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(META)
    background = sorted(meta.symbol.dropna().astype(str).unique())
    core = meta[meta.core_module.astype(bool)].copy()
    sets = {
        "core_positive": sorted(core.loc[core.adjusted_meta_r > 0, "symbol"].astype(str).unique()),
        "core_negative": sorted(core.loc[core.adjusted_meta_r < 0, "symbol"].astype(str).unique()),
    }
    all_rows = []
    for label, genes in sets.items():
        all_rows.extend(query_gprofiler(label, genes, background))
    pd.DataFrame(all_rows).sort_values(["query_set", "fdr"]).to_csv(
        OUT / "pde8b_core_pathway_enrichment.csv", index=False, quoting=csv.QUOTE_MINIMAL
    )
    summary = {
        "service": "g:Profiler gost API",
        "organism": "hsapiens",
        "sources": SOURCES,
        "correction": "FDR",
        "custom_background_n": len(background),
        "query_sizes": {key: len(value) for key, value in sets.items()},
        "significant_terms": {
            key: sum(row["query_set"] == key for row in all_rows) for key in sets
        },
    }
    (OUT / "enrichment_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
