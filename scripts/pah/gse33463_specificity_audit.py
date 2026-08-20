import csv
import json

import numpy as np
import pandas as pd

from formal_meta_and_confounding import MODULES, OUT, auc_rank, module_scores, prepare_gene_matrix


def main():
    with (OUT / "core_sample_manifest.csv").open(encoding="utf-8-sig") as handle:
        records = [row for row in csv.DictReader(handle) if row["dataset"] == "GSE33463"]
    groups = pd.Series({row["gsm"]: row["group"] for row in records})
    matrix, _ = prepare_gene_matrix("GSE33463", "GPL6947")
    matrix = matrix[groups.index]
    module, used = module_scores(matrix)

    features = pd.DataFrame(index=groups.index)
    for gene in ["NRCAM", "HBQ1", "HBM", "HBB"]:
        values = matrix.loc[gene]
        features[gene] = (values - values.mean()) / values.std(ddof=1)
    features["erythroid_module"] = module["erythroid"]
    features["platelet_module"] = module["platelet"]
    features["monocyte_module"] = module["monocyte"]
    features["interferon_module"] = module["interferon"]
    features["group"] = groups

    comparisons = [
        ("SSc-PAH", "SSc-noPH"),
        ("IPAH", "Healthy"),
        ("SSc-PAH", "Healthy"),
        ("SSc-noPH", "Healthy"),
        ("SSc-PH-ILD", "SSc-PAH"),
    ]
    result = {"group_counts": groups.value_counts().to_dict(), "comparisons": {}}
    for case, control in comparisons:
        subset = features[features["group"].isin([case, control])]
        y = (subset["group"] == case).astype(int)
        key = f"{case}_vs_{control}"
        result["comparisons"][key] = {}
        for feature in features.columns.drop("group"):
            case_mean = float(subset.loc[y == 1, feature].mean())
            control_mean = float(subset.loc[y == 0, feature].mean())
            result["comparisons"][key][feature] = {
                "auc_case_high": auc_rank(y, subset[feature]),
                "case_mean_z": case_mean,
                "control_mean_z": control_mean,
                "difference": case_mean - control_mean,
            }
    features.to_csv(OUT / "GSE33463_specificity_features.csv", encoding="utf-8-sig")
    (OUT / "GSE33463_specificity_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

