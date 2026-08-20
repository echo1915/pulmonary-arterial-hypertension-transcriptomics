import json

import numpy as np
import pandas as pd

from formal_meta_and_confounding import MODULES, OUT, hedges_g
from lung_mechanism_meta import prepare_array_cohort, prepare_gene_matrix, read_geo_metadata


SEX_LINKED_EXCLUSIONS = {
    "SRY", "UTY", "KDM5D", "DDX3Y", "EIF1AY", "RPS4Y1", "ZFY", "USP9Y", "TMSB4Y",
    "NLGN4Y", "PRKY", "TTTY14", "XIST", "TSIX",
}


def main():
    meta = pd.read_csv(OUT / "lung_mechanism_random_effects_meta.csv")
    robust = meta[(meta["same_direction"]) & (meta["fdr"] < 0.05) & (meta["I2"] < 50) & (meta["meta_g"].abs() >= 0.50)].copy()

    matrix, _ = prepare_gene_matrix("GSE15197", "GPL6480")
    rows = read_geo_metadata("GSE15197")
    pah = [r["gsm"] for r in rows if "pulmonary arterial hypertension" in r["text"]]
    ipf_ph = [r["gsm"] for r in rows if "ipf" in r["text"]]
    g_specificity = hedges_g(matrix[pah], matrix[ipf_ph])
    robust["g_PAH_vs_IPF_PH_GSE15197"] = robust["symbol"].map(g_specificity)
    robust["specificity_same_direction"] = np.sign(robust["meta_g"]) == np.sign(robust["g_PAH_vs_IPF_PH_GSE15197"])
    robust["specificity_abs_g_ge_0.30"] = robust["g_PAH_vs_IPF_PH_GSE15197"].abs() >= 0.30

    validation_matrix, validation_case, validation_control, _ = prepare_array_cohort("GSE53408")
    validation_g = hedges_g(validation_matrix[validation_case], validation_matrix[validation_control])
    robust["g_GSE53408_secondary"] = robust["symbol"].map(validation_g)
    robust["GSE53408_present"] = robust["g_GSE53408_secondary"].notna()
    robust["GSE53408_direction_concordant"] = np.sign(robust["meta_g"]) == np.sign(robust["g_GSE53408_secondary"])

    composition = set().union(*[set(genes) for genes in MODULES.values()])
    robust["named_blood_composition_gene"] = robust["symbol"].isin(composition)
    robust["sex_linked_exclusion"] = robust["symbol"].isin(SEX_LINKED_EXCLUSIONS)
    robust["priority_pass"] = (
        robust["specificity_same_direction"]
        & robust["specificity_abs_g_ge_0.30"]
        & ~robust["named_blood_composition_gene"]
        & ~robust["sex_linked_exclusion"]
    )
    robust["secondary_support"] = robust["GSE53408_present"] & robust["GSE53408_direction_concordant"]
    robust = robust.sort_values(["priority_pass", "secondary_support", "fdr"], ascending=[False, False, True])
    robust.to_csv(OUT / "lung_mechanism_candidate_prioritization.csv", index=False, encoding="utf-8-sig")

    priority = robust[robust["priority_pass"]]
    summary = {
        "starting_robust_candidates": int(len(robust)),
        "pass_PAH_vs_IPF_PH_specificity_and_basic_confounder_filters": int(len(priority)),
        "priority_with_GSE53408_secondary_support": int(priority["secondary_support"].sum()),
        "top_priority": priority.head(40).to_dict("records"),
    }
    (OUT / "lung_mechanism_candidate_prioritization.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

