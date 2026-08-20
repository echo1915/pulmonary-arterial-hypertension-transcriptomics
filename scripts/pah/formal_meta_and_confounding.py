import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from pilot_effect_consistency import DATA, OUT, read_annotation, read_matrix, choose_baseline_samples, hedges_g


PLATFORM_BY_GSE = {"GSE19617": "GPL6480", "GSE22356": "GPL570", "GSE33463": "GPL6947"}

MODULES = {
    "erythroid": ["AHSP", "ALAS2", "ANK1", "BPGM", "CA1", "DMTN", "EPB41", "EPB42", "FECH", "GYPA", "GYPB", "HBA1", "HBA2", "HBB", "HBD", "HBM", "HBQ1", "KLF1", "SLC4A1", "SLC25A37", "SLC25A39", "TMOD1", "TRIM58"],
    "platelet": ["PPBP", "PF4", "NRGN", "GNG11", "RGS18", "GP9", "ITGA2B", "GP1BA", "TUBB1", "TREML1"],
    "monocyte": ["CTSS", "FCN1", "LYZ", "S100A8", "S100A9", "CTSD", "LILRB1", "TYMP", "FCGR3A", "MS4A7"],
    "t_cell": ["CD3D", "CD3E", "CD247", "TRAC", "LCK", "MAL", "IL7R", "LTB", "TRBC1", "TRBC2"],
    "b_cell": ["CD79A", "CD79B", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD22", "BANK1", "BLK", "CD19"],
    "interferon": ["IFI27", "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "OAS3", "RSAD2", "XAF1"],
}


def auc_rank(labels, scores):
    labels = np.asarray(labels)
    ranks = pd.Series(np.asarray(scores)).rank(method="average").to_numpy()
    n1 = int((labels == 1).sum())
    n0 = int((labels == 0).sum())
    return float((ranks[labels == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def random_effects(g, variance):
    ok = np.isfinite(g) & np.isfinite(variance) & (variance > 0)
    g, variance = g[ok], variance[ok]
    k = len(g)
    if k < 2:
        return (np.nan,) * 8
    w = 1 / variance
    fixed = np.sum(w * g) / np.sum(w)
    q = np.sum(w * (g - fixed) ** 2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    wr = 1 / (variance + tau2)
    pooled = np.sum(wr * g) / np.sum(wr)
    se = math.sqrt(1 / np.sum(wr))
    z = pooled / se
    p = math.erfc(abs(z) / math.sqrt(2))
    ci_low, ci_high = pooled - 1.96 * se, pooled + 1.96 * se
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    return pooled, se, ci_low, ci_high, p, tau2, q, i2


def prepare_gene_matrix(gse, gpl):
    matrix = read_matrix(gse)
    raw_quantiles = matrix.stack().quantile([0, 0.01, 0.5, 0.99, 1]).to_dict()
    transformed = False
    if raw_quantiles[0.99] > 100:
        matrix = np.log2(matrix.clip(lower=0) + 1)
        transformed = True
    annotation = read_annotation(gpl)
    common_ids = matrix.index.intersection(annotation["ID"])
    matrix = matrix.loc[common_ids]
    mapping = annotation.set_index("ID").loc[common_ids, "symbol"]
    matrix["symbol"] = mapping.values
    gene_matrix = matrix.groupby("symbol").median(numeric_only=True)
    return gene_matrix, {"raw_quantiles": {str(k): float(v) for k, v in raw_quantiles.items()}, "log2_applied": transformed}


def module_scores(matrix):
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    scores = {}
    used = {}
    for name, genes in MODULES.items():
        present = [g for g in genes if g in z.index]
        used[name] = present
        scores[name] = z.loc[present].mean(axis=0) if present else pd.Series(np.nan, index=z.columns)
    return pd.DataFrame(scores), used


def main():
    with (OUT / "core_sample_manifest.csv").open(encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle))

    effects = {}
    variances = {}
    matrices = {}
    labels = {}
    qc = {}
    module_audits = {}

    for gse, gpl in PLATFORM_BY_GSE.items():
        matrix, scale_qc = prepare_gene_matrix(gse, gpl)
        selected = choose_baseline_samples(manifest, gse)
        case = [r["gsm"] for r in selected if r["group"] == "SSc-PAH"]
        control = [r["gsm"] for r in selected if r["group"] == "SSc-noPH"]
        matrix = matrix[case + control]
        y = pd.Series([1] * len(case) + [0] * len(control), index=case + control)
        g = hedges_g(matrix[case], matrix[control])
        n1, n0 = len(case), len(control)
        var = (n1 + n0) / (n1 * n0) + (g**2) / (2 * (n1 + n0 - 2))
        effects[gse], variances[gse] = g, var
        matrices[gse], labels[gse] = matrix, y
        scores, used = module_scores(matrix)
        module_audits[gse] = {
            name: {"genes_used": used[name], "auc": auc_rank(y, scores[name])}
            for name in MODULES
        }
        qc[gse] = {**scale_qc, "n_case": n1, "n_control": n0, "gene_count": len(matrix)}

    common = set.intersection(*(set(series.index) for series in effects.values()))
    rows = []
    for gene in sorted(common):
        gs = np.array([effects[gse].loc[gene] for gse in PLATFORM_BY_GSE], dtype=float)
        vs = np.array([variances[gse].loc[gene] for gse in PLATFORM_BY_GSE], dtype=float)
        pooled, se, lo, hi, p, tau2, q, i2 = random_effects(gs, vs)
        rows.append({
            "symbol": gene,
            **{f"g_{gse}": gs[i] for i, gse in enumerate(PLATFORM_BY_GSE)},
            "n_effects": int(np.isfinite(gs).sum()),
            "same_direction": bool(np.isfinite(gs).all() and len(set(np.sign(gs))) == 1),
            "meta_g": pooled,
            "meta_se": se,
            "ci_low": lo,
            "ci_high": hi,
            "p": p,
            "tau2": tau2,
            "Q": q,
            "I2": i2,
        })
    meta = pd.DataFrame(rows)
    meta["fdr"] = bh_adjust(meta["p"])
    meta["is_erythroid_gene"] = meta["symbol"].isin(MODULES["erythroid"])
    meta = meta.sort_values(["fdr", "I2", "meta_g"], ascending=[True, True, False])
    meta.to_csv(OUT / "formal_random_effects_meta.csv", index=False, encoding="utf-8-sig")

    robust = meta[(meta["n_effects"] == 3) & (meta["same_direction"]) & (meta["fdr"] < 0.05) & (meta["I2"] < 50) & (meta["meta_g"].abs() >= 0.30)]
    robust_nonery = robust[~robust["is_erythroid_gene"]]
    summary = {
        "qc": qc,
        "module_group_discrimination_auc": module_audits,
        "common_genes": int(len(meta)),
        "fdr_lt_0.05": int((meta["fdr"] < 0.05).sum()),
        "same_direction_fdr_lt_0.05_i2_lt_50_abs_g_ge_0.30": int(len(robust)),
        "robust_after_named_erythroid_gene_exclusion": int(len(robust_nonery)),
        "top_robust_nonerythroid": robust_nonery.head(30).to_dict("records"),
    }
    (OUT / "formal_meta_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

