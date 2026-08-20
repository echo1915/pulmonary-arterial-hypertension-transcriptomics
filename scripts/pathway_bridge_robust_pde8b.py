"""Compare pathway enrichment of the robust PAH signature and PDE8B core program."""

from __future__ import annotations

import csv
import json
import ssl
import sys
import urllib.request
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import certifi


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402


META = OUTPUT_ROOT / "pah_audit" / "lung_mechanism_random_effects_meta.csv"
COEX_DIR = PROJECT_DATA_ROOT / "outputs" / "pde8b-coexpression"
EXISTING = COEX_DIR / "enrichment" / "pde8b_core_pathway_enrichment.csv"
OUT = PROJECT_DATA_ROOT / "outputs" / "pathway-bridge"
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "REAC", "KEGG", "WP"]

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def query(label: str, genes: list[str], background: list[str]) -> list[dict]:
    payload = {
        "organism": "hsapiens", "query": genes, "sources": SOURCES,
        "user_threshold": 0.05, "significance_threshold_method": "fdr",
        "domain_scope": "custom", "background": background, "no_evidences": False,
    }
    request = urllib.request.Request(URL, data=json.dumps(payload).encode("utf-8"),
                                     headers={"Content-Type": "application/json", "User-Agent": "PAH-PDE8B-pathway-bridge/1.0"}, method="POST")
    ssl_context = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(request, timeout=120, context=ssl_context) as response:
        result = json.loads(response.read().decode("utf-8"))
    (OUT / f"{label}_gprofiler_raw.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    query_genes = result.get("meta", {}).get("query_metadata", {}).get("queries", {}).get("query_1", genes)
    rows = []
    for row in result.get("result", []):
        hits = [query_genes[i] for i, evidence in enumerate(row.get("intersections", [])) if evidence]
        rows.append({
            "query_set": label, "source": row.get("source"), "term_id": row.get("native"),
            "term_name": row.get("name"), "fdr": row.get("p_value"), "term_size": row.get("term_size"),
            "query_size": row.get("query_size"), "intersection_size": row.get("intersection_size"),
            "precision": row.get("precision"), "recall": row.get("recall"), "intersection": ";".join(hits),
        })
    return rows


def main() -> None:
    meta = pd.read_csv(META)
    coex = pd.read_csv(COEX_DIR / "pde8b_cross_cohort_coexpression_meta.csv")
    background_series = coex.symbol
    background_missing_n = int(background_series.isna().sum())
    background = sorted(set(background_series[background_series.notna()].astype(str)))
    robust = meta[meta.same_direction.astype(bool) & (meta.fdr < 0.05) & (meta.I2 < 50) & (meta.meta_g.abs() >= 0.50)]
    robust = robust[robust.symbol.isin(background)]
    sets = {
        "robust_positive": sorted(robust.loc[robust.meta_g > 0, "symbol"].astype(str)),
        "robust_negative": sorted(robust.loc[robust.meta_g < 0, "symbol"].astype(str)),
    }
    rows = []
    for label, genes in sets.items():
        rows.extend(query(label, genes, background))
    robust_enrichment = pd.DataFrame(rows)
    robust_enrichment.to_csv(OUT / "robust_signature_pathway_enrichment.csv", index=False, quoting=csv.QUOTE_MINIMAL)

    core = pd.read_csv(EXISTING)
    combined = pd.concat([robust_enrichment, core], ignore_index=True)
    combined.to_csv(OUT / "combined_signature_core_pathway_enrichment.csv", index=False, quoting=csv.QUOTE_MINIMAL)

    comparisons = [("positive", "robust_positive", "core_positive"), ("negative", "robust_negative", "core_negative")]
    shared_rows = []
    summary_rows = []
    for direction, robust_label, core_label in comparisons:
        a = robust_enrichment[robust_enrichment.query_set.eq(robust_label)].drop_duplicates("term_id")
        b = core[core.query_set.eq(core_label)].drop_duplicates("term_id")
        shared = a.merge(b, on=["source", "term_id"], suffixes=("_robust", "_core"))
        shared["direction"] = direction
        shared["bridge_score"] = -np.log10(shared.fdr_robust.clip(lower=1e-300)) - np.log10(shared.fdr_core.clip(lower=1e-300))
        shared_rows.append(shared)
        union_n = len(set(a.term_id) | set(b.term_id))
        summary_rows.append({
            "direction": direction, "robust_terms": len(a), "core_terms": len(b),
            "shared_exact_terms": len(shared), "term_jaccard": len(shared) / union_n if union_n else 0,
        })
    shared = pd.concat(shared_rows, ignore_index=True) if shared_rows else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    shared.to_csv(OUT / "shared_exact_pathway_terms.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "pathway_bridge_summary.csv", index=False, encoding="utf-8-sig")

    themes = {
        "Muscle/cytoskeleton": r"muscle|contract|actin|cytoskeleton|myofibril|sarcomere",
        "ECM/adhesion": r"extracellular matrix|collagen|adhesion|focal adhesion",
        "Cilium/projection": r"cili|cell projection|microtubule",
        "Vascular signaling": r"cAMP|cGMP|PKG|calcium|vascular|smooth muscle",
        "Stress/protein": r"stress|translation|ribosome|protein",
    }
    theme_rows = []
    for query_set, block in combined.groupby("query_set"):
        names = block.term_name.fillna("")
        for theme, pattern in themes.items():
            hit = block[names.str.contains(pattern, case=False, regex=True)]
            theme_rows.append({
                "query_set": query_set, "theme": theme, "terms_n": len(hit),
                "best_fdr": hit.fdr.min() if len(hit) else np.nan,
                "best_term": hit.sort_values("fdr").term_name.iloc[0] if len(hit) else "",
            })
    theme = pd.DataFrame(theme_rows)
    theme.to_csv(OUT / "pathway_theme_bridge_matrix.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(11, 4), gridspec_kw={"width_ratios": [0.85, 1.45, 1.15]})
    qsets = ["robust_positive", "core_positive", "robust_negative", "core_negative"]
    labels = ["Robust +", "PDE8B core +", "Robust -", "PDE8B core -"]
    counts = [int((combined.query_set == q).sum()) for q in qsets]
    axes[0].barh(np.arange(4), counts, color=["#477AA8", "#D97732", "#7D8992", "#B89A78"])
    axes[0].set_yticks(np.arange(4), labels); axes[0].invert_yaxis()
    axes[0].set_xlabel("Significant enriched terms")
    axes[0].set_title("A  Enrichment breadth", loc="left", fontweight="bold")

    matrix = theme.pivot(index="theme", columns="query_set", values="best_fdr").reindex(columns=qsets)
    values = -np.log10(matrix.clip(lower=1e-12)).fillna(0)
    image = axes[1].imshow(values, cmap="Blues", aspect="auto", vmin=0)
    axes[1].set_xticks(np.arange(4), labels, rotation=25, ha="right")
    axes[1].set_yticks(np.arange(len(values)), values.index)
    axes[1].set_title("B  Convergence at biological-theme level", loc="left", fontweight="bold")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            axes[1].text(j, i, f"{values.iloc[i, j]:.1f}" if values.iloc[i, j] else "-", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=axes[1], fraction=0.04, pad=0.03, label="-log10(best FDR)")

    top = shared.sort_values("bridge_score", ascending=False).head(10) if len(shared) else shared
    if len(top):
        y = np.arange(len(top))
        axes[2].barh(y, top.bridge_score, color=np.where(top.direction.eq("positive"), "#D97732", "#7D8992"))
        axes[2].set_yticks(y, top.term_name_robust.str.slice(0, 42))
        axes[2].invert_yaxis(); axes[2].set_xlabel("Combined enrichment evidence")
    else:
        axes[2].text(0.5, 0.5, "No exact shared terms", ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_xticks([]); axes[2].set_yticks([])
    axes[2].set_title("C  Exact shared pathway terms", loc="left", fontweight="bold")

    fig.suptitle("Disease differential genes and the PDE8B program are compared at pathway level", x=0.03, ha="left", fontsize=12.2, fontweight="bold")
    fig.text(0.03, 0.015, "g:Profiler FDR enrichment uses the same 15,015-gene custom background for all query sets. Theme matching is descriptive; exact term overlap is reported separately.", fontsize=7, color="#56616B")
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    fig.savefig(OUT / "Figure_pathway_bridge_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_pathway_bridge_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_pathway_bridge_draft.svg", bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({"query_sizes": {k: len(v) for k, v in sets.items()}, "summary": summary_rows}, indent=2))


if __name__ == "__main__":
    main()


