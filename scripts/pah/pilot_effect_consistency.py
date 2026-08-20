import csv
import gzip
import json
import math
import re
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_paths import DATA_ROOT, OUTPUT_ROOT

import numpy as np
import pandas as pd


DATA = DATA_ROOT / "pah_geo"
OUT = OUTPUT_ROOT / "pah_audit"


def read_annotation(gpl):
    path = DATA / f"{gpl}.annot.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = [line for line in handle if not line.startswith(("#", "!", "^"))]
    from io import StringIO
    annotation = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False)
    symbol_col = "Gene Symbol" if "Gene Symbol" in annotation.columns else "Gene symbol"
    annotation = annotation[["ID", symbol_col]].rename(columns={symbol_col: "symbol"}).dropna()
    annotation = annotation[~annotation["symbol"].str.contains("///|//|,|;", regex=True)]
    annotation = annotation[annotation["symbol"] != "---"]
    annotation["symbol"] = annotation["symbol"].str.upper()
    return annotation.drop_duplicates("ID")


def read_matrix(gse):
    path = DATA / f"{gse}_series_matrix.txt.gz"
    rows = []
    header = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = parts
            else:
                rows.append(parts)
    frame = pd.DataFrame(rows, columns=header).set_index("ID_REF")
    return frame.apply(pd.to_numeric, errors="coerce")


def choose_baseline_samples(manifest, gse):
    records = [row for row in manifest if row["dataset"] == gse and row["group"] in {"SSc-PAH", "SSc-noPH"}]
    records = [row for row in records if row["followup"].lower() != "true"]
    selected = {}
    for row in records:
        key = (row["group"], row["subject_key"].lower())
        if key not in selected or (selected[key]["technical_repeat"].lower() == "true" and row["technical_repeat"].lower() == "false"):
            selected[key] = row
    return list(selected.values())


def hedges_g(case, control):
    n1, n0 = case.shape[1], control.shape[1]
    m1, m0 = case.mean(axis=1), control.mean(axis=1)
    v1, v0 = case.var(axis=1, ddof=1), control.var(axis=1, ddof=1)
    pooled = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
    d = (m1 - m0) / pooled.replace(0, np.nan)
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return d * correction


def main():
    with (OUT / "core_sample_manifest.csv").open(encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle))
    platform_by_gse = {"GSE19617": "GPL6480", "GSE22356": "GPL570", "GSE33463": "GPL6947"}
    effects = {}
    sample_counts = {}
    cohort_matrices = {}
    cohort_groups = {}

    for gse, gpl in platform_by_gse.items():
        matrix = read_matrix(gse)
        annotation = read_annotation(gpl)
        common_ids = matrix.index.intersection(annotation["ID"])
        matrix = matrix.loc[common_ids]
        mapping = annotation.set_index("ID").loc[common_ids, "symbol"]
        matrix["symbol"] = mapping.values
        gene_matrix = matrix.groupby("symbol").median(numeric_only=True)

        selected = choose_baseline_samples(manifest, gse)
        case_ids = [row["gsm"] for row in selected if row["group"] == "SSc-PAH"]
        control_ids = [row["gsm"] for row in selected if row["group"] == "SSc-noPH"]
        case_ids = [x for x in case_ids if x in gene_matrix.columns]
        control_ids = [x for x in control_ids if x in gene_matrix.columns]
        case, control = gene_matrix[case_ids], gene_matrix[control_ids]
        g = hedges_g(case, control)
        effects[gse] = pd.DataFrame({f"g_{gse}": g}, index=gene_matrix.index)
        sample_counts[gse] = {"SSc-PAH": len(case_ids), "SSc-noPH": len(control_ids), "genes": len(gene_matrix)}
        cohort_matrices[gse] = gene_matrix[case_ids + control_ids]
        cohort_groups[gse] = pd.Series([1] * len(case_ids) + [0] * len(control_ids), index=case_ids + control_ids)

    merged = pd.concat(effects.values(), axis=1, join="inner").dropna()
    gcols = [f"g_{gse}" for gse in platform_by_gse]
    signs = np.sign(merged[gcols])
    same_direction = (signs.nunique(axis=1) == 1)
    moderate_all = same_direction & (merged[gcols].abs().min(axis=1) >= 0.30)
    strong_two = same_direction & ((merged[gcols].abs() >= 0.50).sum(axis=1) >= 2)

    correlations = {}
    for i, left in enumerate(gcols):
        for right in gcols[i + 1:]:
            left_rank = merged[left].rank().to_numpy()
            right_rank = merged[right].rank().to_numpy()
            rho = np.corrcoef(left_rank, right_rank)[0, 1]
            correlations[f"{left}_vs_{right}"] = {"rho": float(rho)}

    def auc_rank(labels, scores):
        labels = np.asarray(labels)
        ranks = pd.Series(scores).rank(method="average").to_numpy()
        n1 = int((labels == 1).sum())
        n0 = int((labels == 0).sum())
        return float((ranks[labels == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

    erythroid = {
        "AHSP", "ALAS2", "ANK1", "BPGM", "CA1", "EPB41", "EPB42", "FECH", "GYPA", "GYPB",
        "HBA1", "HBA2", "HBB", "HBD", "HBM", "HBQ1", "KLF1", "SLC4A1", "SLC25A37",
        "SLC25A39", "TMOD1", "TRIM10", "TRIM58", "DMTN",
    }
    loco = {}
    for heldout in platform_by_gse:
        training = [gse for gse in platform_by_gse if gse != heldout]
        train = pd.concat([effects[gse] for gse in training], axis=1, join="inner").dropna()
        train_cols = [f"g_{gse}" for gse in training]
        eligible = train[(np.sign(train[train_cols]).nunique(axis=1) == 1) & (train[train_cols].abs().min(axis=1) >= 0.30)].copy()
        eligible["min_abs_g"] = eligible[train_cols].abs().min(axis=1)
        eligible["direction"] = np.sign(eligible[train_cols].mean(axis=1))
        eligible = eligible.sort_values("min_abs_g", ascending=False)
        hold = cohort_matrices[heldout]
        labels = cohort_groups[heldout]
        results = {}
        for remove_erythroid in (False, True):
            pool = eligible[~eligible.index.isin(erythroid)] if remove_erythroid else eligible
            scenario = "erythroid_excluded" if remove_erythroid else "all_genes"
            results[scenario] = {"eligible_genes": int(len(pool)), "models": {}}
            for k in (3, 5, 10, 20):
                genes = [gene for gene in pool.index if gene in hold.index][:k]
                if len(genes) < k:
                    continue
                z = hold.loc[genes].sub(hold.loc[genes].mean(axis=1), axis=0)
                sd = hold.loc[genes].std(axis=1, ddof=1).replace(0, np.nan)
                z = z.div(sd, axis=0).fillna(0)
                directions = pool.loc[genes, "direction"]
                score = z.mul(directions, axis=0).mean(axis=0)
                results[scenario]["models"][str(k)] = {
                    "auc": auc_rank(labels.loc[score.index].to_numpy(), score.to_numpy()),
                    "genes": genes,
                }
        loco[heldout] = results

    merged["mean_g"] = merged[gcols].mean(axis=1)
    merged["min_abs_g"] = merged[gcols].abs().min(axis=1)
    merged["same_direction"] = same_direction
    candidates = merged.loc[moderate_all | strong_two].sort_values("min_abs_g", ascending=False)
    merged.to_csv(OUT / "pilot_all_common_gene_effects.csv", encoding="utf-8-sig")
    candidates.to_csv(OUT / "pilot_consistent_candidates.csv", encoding="utf-8-sig")

    summary = {
        "sample_counts": sample_counts,
        "common_tested_genes": int(len(merged)),
        "same_direction_genes": int(same_direction.sum()),
        "same_direction_fraction": float(same_direction.mean()),
        "moderate_in_all_three_abs_g_ge_0.30": int(moderate_all.sum()),
        "same_direction_abs_g_ge_0.50_in_at_least_two": int(strong_two.sum()),
        "pairwise_effect_correlations": correlations,
        "leave_one_cohort_out_scores": loco,
        "top_candidates_by_minimum_effect": candidates.head(20).reset_index().to_dict("records"),
    }
    (OUT / "pilot_effect_consistency.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()



