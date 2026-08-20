"""Audit patient-level PDE8B pseudobulk evidence in two PAH scRNA-seq cohorts."""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402


AUDIT = OUTPUT_ROOT / "pah_audit"
OUT = PROJECT_DATA_ROOT / "outputs" / "scrna-patient-pseudobulk-audit"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"Donor": "#477AA8", "PAH": "#B84343", "IPAH": "#B84343", "SSc-PAH": "#D97732"}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def hedges_g(case: np.ndarray, control: np.ndarray) -> float:
    n1, n0 = len(case), len(control)
    pooled = math.sqrt(((n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / (n1 + n0 - 2))
    if not np.isfinite(pooled) or pooled == 0:
        return 0.0
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return correction * (np.mean(case) - np.mean(control)) / pooled


def exact_permutation_p(case: np.ndarray, control: np.ndarray) -> float:
    values = np.r_[case, control]
    observed = abs(np.mean(case) - np.mean(control))
    stats = []
    for idx in itertools.combinations(range(len(values)), len(case)):
        idx = set(idx)
        a = np.array([values[i] for i in range(len(values)) if i in idx])
        b = np.array([values[i] for i in range(len(values)) if i not in idx])
        stats.append(abs(np.mean(a) - np.mean(b)))
    return float(np.mean(np.asarray(stats) >= observed - 1e-12))


def leave_one_direction(case: np.ndarray, control: np.ndarray) -> tuple[int, int]:
    full_sign = np.sign(np.mean(case) - np.mean(control))
    signs = []
    for i in range(len(case)):
        signs.append(np.sign(np.mean(np.delete(case, i)) - np.mean(control)))
    for i in range(len(control)):
        signs.append(np.sign(np.mean(case) - np.mean(np.delete(control, i))))
    return int(np.sum(np.asarray(signs) == full_sign)), len(signs)


def summarize(block: pd.DataFrame, group_col: str, case_name: str, control_name: str, label: str) -> dict:
    case = block.loc[block[group_col].eq(case_name), "log_cp10k"].to_numpy()
    control = block.loc[block[group_col].eq(control_name), "log_cp10k"].to_numpy()
    concordant, total = leave_one_direction(case, control)
    return {
        "analysis": label,
        "case": case_name,
        "control": control_name,
        "n_case": len(case),
        "n_control": len(control),
        "delta_log_cp10k": float(np.mean(case) - np.mean(control)),
        "hedges_g": float(hedges_g(case, control)),
        "exact_two_sided_permutation_p": exact_permutation_p(case, control),
        "loo_direction_concordant_n": concordant,
        "loo_iterations_n": total,
    }


def main() -> None:
    gse210 = pd.read_csv(AUDIT / "scrna_gse210248_candidate_pseudobulk.csv")
    gse210 = gse210[gse210.gene.eq("PDE8B") & gse210.cell_type.isin(["SMC 1", "SMC 2"])].copy()
    gse293 = pd.read_csv(AUDIT / "scrna_gse293580_candidate_pseudobulk.csv")
    gse293 = gse293[gse293.gene.eq("PDE8B")].copy()

    rows = []
    for cell_type in ["SMC 1", "SMC 2"]:
        rows.append(summarize(gse210[gse210.cell_type.eq(cell_type)], "disease", "PAH", "Donor", f"GSE210248 {cell_type}"))
    rows.append(summarize(gse293, "condition", "IPAH", "Donor", "GSE293580 whole-sample IPAH"))
    rows.append(summarize(gse293, "condition", "SSc-PAH", "Donor", "GSE293580 whole-sample SSc-PAH"))
    summary = pd.DataFrame(rows)

    carrier = (gse210.assign(detected=lambda x: x["count"] > 0)
               .groupby(["cell_type", "disease"], as_index=False)
               .agg(subjects=("sample", "nunique"), detected_subjects=("detected", "sum"),
                    min_cells=("n_cells", "min"), median_cells=("n_cells", "median"), max_cells=("n_cells", "max")))
    carrier["detected_subject_fraction"] = carrier.detected_subjects / carrier.subjects

    summary.to_csv(OUT / "pde8b_patient_level_effect_summary.csv", index=False, encoding="utf-8-sig")
    gse210.to_csv(OUT / "gse210248_pde8b_smc_patient_values.csv", index=False, encoding="utf-8-sig")
    gse293.to_csv(OUT / "gse293580_pde8b_whole_sample_values.csv", index=False, encoding="utf-8-sig")
    carrier.to_csv(OUT / "gse210248_pde8b_smc_detection_audit.csv", index=False, encoding="utf-8-sig")
    report = {
        "replicate_unit": "patient/sample pseudobulk",
        "gse210248_limitation": "Only 3 PAH and 3 donor subjects; SMC2 signal is sparse and concentrated in PAH_1.",
        "gse293580_scope": "Whole-sample scRNA-seq pseudobulk; no cell-type localization is inferred from this table.",
        "inference_boundary": "Patient points support direction/replication; small n and sparse expression limit formal significance and mechanistic claims.",
    }
    (OUT / "pde8b_patient_pseudobulk_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7), gridspec_kw={"width_ratios": [1.05, 1.05, 1.35]})
    def fixed_jitter(n: int) -> np.ndarray:
        return np.linspace(-0.055, 0.055, n) if n > 1 else np.zeros(1)
    for ax, cell_type in zip(axes[:2], ["SMC 1", "SMC 2"]):
        z = gse210[gse210.cell_type.eq(cell_type)]
        for xpos, condition in enumerate(["Donor", "PAH"]):
            q = z[z.disease.eq(condition)]
            jitter = fixed_jitter(len(q))
            sizes = 28 + 90 * np.sqrt(q.n_cells / q.n_cells.max())
            ax.scatter(xpos + jitter, q.log_cp10k, s=sizes, color=COLORS[condition], alpha=0.78, edgecolor="white", lw=0.5)
            for xx, (_, row) in zip(xpos + jitter, q.iterrows()):
                ax.annotate(row["sample"], (xx, row.log_cp10k), xytext=(2, 3), textcoords="offset points", fontsize=6)
            ax.hlines(q.log_cp10k.mean(), xpos - 0.18, xpos + 0.18, color="#27323A", lw=1.5)
        stat = summary[summary.analysis.eq(f"GSE210248 {cell_type}")].iloc[0]
        ax.set_xticks([0, 1], ["Donor", "PAH"])
        ax.set_ylabel("Normalized PDE8B abundance")
        ax.set_title(f"{'A' if cell_type == 'SMC 1' else 'B'}  GSE210248 {cell_type}", loc="left", fontweight="bold")
        ax.text(0.03, 0.97, f"g={stat.hedges_g:.2f}\nexact P={stat.exact_two_sided_permutation_p:.2f}", transform=ax.transAxes, va="top", fontsize=7)
    axes[1].text(0.5, -0.23, "Point size reflects cells in each patient pseudobulk", transform=axes[1].transAxes, ha="center", fontsize=6.5, color="#56616B")

    ax = axes[2]
    conditions = ["Donor", "IPAH", "SSc-PAH"]
    for xpos, condition in enumerate(conditions):
        q = gse293[gse293.condition.eq(condition)]
        jitter = fixed_jitter(len(q))
        sizes = 28 + 90 * np.sqrt(q.n_cells / gse293.n_cells.max())
        ax.scatter(xpos + jitter, q.log_cp10k, s=sizes, color=COLORS[condition], alpha=0.78, edgecolor="white", lw=0.5)
        for xx, (_, row) in zip(xpos + jitter, q.iterrows()):
            ax.annotate(row["sample"], (xx, row.log_cp10k), xytext=(2, 3), textcoords="offset points", fontsize=6)
        ax.hlines(q.log_cp10k.mean(), xpos - 0.18, xpos + 0.18, color="#27323A", lw=1.5)
    ax.set_xticks(range(3), conditions)
    ax.set_ylabel("Normalized PDE8B abundance")
    ax.set_title("C  GSE293580 whole-sample pseudobulk", loc="left", fontweight="bold")
    ipah = summary[summary.analysis.str.contains("whole-sample IPAH")].iloc[0]
    ssc = summary[summary.analysis.str.contains("whole-sample SSc")].iloc[0]
    ax.text(0.03, 0.97, f"IPAH: g={ipah.hedges_g:.2f}, exact P={ipah.exact_two_sided_permutation_p:.3f}\nSSc-PAH: g={ssc.hedges_g:.2f}, exact P={ssc.exact_two_sided_permutation_p:.3f}", transform=ax.transAxes, va="top", fontsize=7)

    fig.suptitle("Patient-level pseudobulk supports replication but weakens SMC2 localization certainty", x=0.04, ha="left", fontsize=12.5, fontweight="bold")
    fig.text(0.04, 0.015, "Each point is one patient/sample. Exact permutation tests are descriptive because cohort sizes are small; cells are not treated as independent replicates.", fontsize=7, color="#56616B")
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    fig.savefig(OUT / "Figure_PDE8B_patient_pseudobulk_audit_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_PDE8B_patient_pseudobulk_audit_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_PDE8B_patient_pseudobulk_audit_draft.svg", bbox_inches="tight")
    plt.close(fig)
    print(summary.to_string(index=False))
    print(carrier.to_string(index=False))


if __name__ == "__main__":
    main()

