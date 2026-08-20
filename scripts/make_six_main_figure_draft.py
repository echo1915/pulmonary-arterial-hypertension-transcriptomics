"""Build the six-figure PDE8B/SMC story draft from existing project outputs.

All inferential panels use cohort- or patient-level effect estimates.  The figures
are an editorial v1: they establish the evidence order and visual grammar before
final panel-level polishing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402
from figure_style import COLORS, figure_title, panel_label, save_figure  # noqa: E402


ROOT = PROJECT_DATA_ROOT
OUTPUTS = ROOT / "outputs"
OUT = OUTPUTS / "manuscript-figures-sci-v1"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = COLORS["ink"]
BLUE = COLORS["smc1"]
TEAL = COLORS["teal"]
ORANGE = COLORS["smc2"]
RED = COLORS["signal"]
GREY = COLORS["muted"]
PALE = COLORS["pale"]
LIGHT_BLUE = COLORS["smc1_light"]

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def panel(ax, label: str) -> None:
    panel_label(ax, label)


def title(fig, number: int, text: str) -> None:
    figure_title(fig, text)


def save(fig, stem: str) -> None:
    fig.subplots_adjust(top=0.90)
    save_figure(fig, OUT, stem)


def box(ax, xy, width, height, headline, detail, edge=BLUE, fill="white"):
    patch = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.012,rounding_size=0.008",
                           ec=edge, fc=fill, lw=1.0)
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height * 0.65, headline,
            ha="center", va="center", fontweight="bold", color=edge, fontsize=8.5)
    ax.text(xy[0] + width / 2, xy[1] + height * 0.31, detail,
            ha="center", va="center", color=GREY, fontsize=6.8, linespacing=1.2)


def fig1_design():
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    title(fig, 1, "Study architecture and evidence-preserving analysis design")
    cohorts = [
        ("GSE15197", "18 PAH / 13 donor"), ("GSE113439", "14 / 11"),
        ("GSE254617", "82 / 52"), ("GSE208592", "15 / 18"),
    ]
    for i, (name, n) in enumerate(cohorts):
        box(ax, (0.035 + i * 0.235, 0.70), 0.19, 0.13, name, n,
            edge=BLUE, fill=LIGHT_BLUE)
    ax.text(0.035, 0.88, "Independent lung bulk cohorts", color=NAVY,
            fontweight="bold", fontsize=10)
    stages = [
        ("Cross-cohort\nmeta-analysis", "subject-level effects"),
        ("Stable mechanism\nmodule", "543 genes"),
        ("High-confidence\ncore", "201 genes; |r| ≥ 0.40"),
        ("PDE8B\nprioritization", "specificity + replication"),
        ("Single-cell\nvalidation", "patient pseudobulk"),
    ]
    for i, (h, d) in enumerate(stages):
        x = 0.025 + i * 0.195
        box(ax, (x, 0.39), 0.165, 0.14, h, d, edge=TEAL if i < 3 else ORANGE)
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch((x + 0.168, 0.46), (x + 0.19, 0.46),
                                         arrowstyle="-|>", mutation_scale=11,
                                         color=NAVY, lw=1))
    box(ax, (0.06, 0.10), 0.25, 0.14, "GSE210248", "SMC1/SMC2 localization\nn = 3 PAH + 3 donor",
        edge=BLUE, fill=LIGHT_BLUE)
    box(ax, (0.375, 0.10), 0.25, 0.14, "GSE293580", "SMC boundary + programs\nIPAH / SSc-PAH / donor",
        edge=ORANGE, fill="#FFF3E7")
    box(ax, (0.69, 0.10), 0.25, 0.14, "Inference boundary", "association + localization\ncausality not established",
        edge=RED, fill="#FBEDED")
    ax.text(0.02, 0.015,
            "Audit rule: GSE117261 overlaps PHBI/GSE254617 and is not counted as an independent primary cohort; GSE53408 is secondary validation.",
            color=GREY, fontsize=7.5)
    save(fig, "Figure_1_study_design_and_cohort_audit")


def fig2_core_signature():
    meta = pd.read_csv(OUTPUTS / "current-results" / "pah_audit" /
                       "lung_mechanism_random_effects_meta.csv")
    p = meta.loc[meta.symbol.eq("PDE8B")].iloc[0]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1.35, 0.9]})
    title(fig, 2, "A reproducible PAH lung program nominates a stable disease core")
    names = ["GSE15197", "GSE113439", "GSE254617", "GSE208592", "Random effects"]
    vals = [p.g_GSE15197, p.g_GSE113439, p.g_GSE254617, p.g_GSE208592, p.meta_g]
    y = np.arange(5)[::-1]
    axs[0].axvline(0, color=GREY, lw=0.8)
    axs[0].scatter(vals[:4], y[:4], s=45, c=BLUE, zorder=3)
    axs[0].errorbar(p.meta_g, y[4], xerr=[[p.meta_g-p.ci_low], [p.ci_high-p.meta_g]],
                    fmt="D", color=RED, capsize=3)
    axs[0].set_yticks(y, names); axs[0].set_xlabel("Hedges g (PAH vs donor)")
    axs[0].set_title("PDE8B across lung cohorts", fontweight="bold"); panel(axs[0], "A")
    m = meta.copy(); m["mlog"] = -np.log10(m.fdr.clip(lower=1e-300))
    robust = m.same_direction & (m.fdr < 0.05) & (m.I2 < 50)
    axs[1].scatter(m.loc[~robust, "meta_g"], m.loc[~robust, "mlog"], s=3,
                   c="#CDD6DC", alpha=0.45, rasterized=True)
    axs[1].scatter(m.loc[robust, "meta_g"], m.loc[robust, "mlog"], s=6,
                   c=TEAL, alpha=0.7, rasterized=True)
    axs[1].scatter(p.meta_g, -np.log10(p.fdr), s=60, c=RED, zorder=4)
    axs[1].annotate("PDE8B", (p.meta_g, -np.log10(p.fdr)), xytext=(6, 6),
                    textcoords="offset points", color=RED, fontweight="bold")
    axs[1].axvline(0, color=GREY, lw=0.7)
    axs[1].set(xlabel="Random-effects Hedges g", ylabel="−log10(FDR)")
    axs[1].set_title("Cross-cohort meta-analysis", fontweight="bold"); panel(axs[1], "B")
    axs[2].axis("off"); panel(axs[2], "C")
    funnel = [("15,015", "common genes"), ("2,413", "FDR < 0.05"),
              ("817", "robust"), ("543", "stable module"), ("201", "high-confidence core")]
    for i, (n, lab) in enumerate(funnel):
        w = 0.82 - i * 0.11; x = 0.5 - w / 2; yy = 0.82 - i * 0.16
        axs[2].add_patch(FancyBboxPatch((x, yy), w, 0.105, boxstyle="round,pad=.01",
                                        fc=PALE if i < 4 else "#FFF0E1",
                                        ec=TEAL if i < 4 else ORANGE, lw=1.2))
        axs[2].text(0.5, yy + 0.064, n, ha="center", fontweight="bold",
                    color=TEAL if i < 4 else ORANGE, fontsize=10)
        axs[2].text(0.5, yy + 0.025, lab, ha="center", color=GREY, fontsize=7)
    save(fig, "Figure_2_cross_cohort_core_signature")


def fig3_prioritization():
    comp = pd.read_csv(OUTPUTS / "candidate-comparison" / "three_candidate_evidence_matrix.csv")
    coex_summary = json.loads((OUTPUTS / "pde8b-coexpression" /
                               "pde8b_coexpression_summary.json").read_text(encoding="utf-8"))
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.05, 1.1, 1]})
    title(fig, 3, "Transparent prioritization identifies PDE8B as the SMC-centered lead candidate")
    metrics = ["main_meta_g", "PAH_vs_IPF_PH_g", "GSE53408_raw_g", "GSE210248_top_celltype_g"]
    labels = ["Primary meta", "PAH vs IPF-PH", "GSE53408", "Top cell type"]
    x = np.arange(len(metrics)); width = 0.23
    colors = [ORANGE, BLUE, TEAL]
    for j, (_, row) in enumerate(comp.iterrows()):
        axs[0].bar(x + (j-1)*width, [row[m] for m in metrics], width,
                   color=colors[j], label=row.gene)
    axs[0].axhline(0, color=GREY, lw=0.8); axs[0].set_xticks(x, labels, rotation=25, ha="right")
    axs[0].set_ylabel("Hedges g"); axs[0].legend(frameon=False)
    axs[0].set_title("Parallel candidate comparison", fontweight="bold"); panel(axs[0], "A")
    evidence = np.array([[1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=float)
    # The last column is not a significance flag: it records SMC as the top localized compartment.
    evidence[1, 3] = 0
    axs[1].imshow(evidence, cmap=matplotlib.colors.ListedColormap(["#E2E7EA", TEAL]),
                  vmin=0, vmax=1, aspect="auto")
    axs[1].set_yticks(range(3), comp.gene)
    axs[1].set_xticks(range(4), ["4-cohort\ndirection", "disease\nspecificity",
                                "secondary\nreplication", "SMC\nlocalization"], rotation=25)
    for i in range(3):
        for j in range(4): axs[1].text(j, i, "●" if evidence[i,j] else "—", ha="center", va="center", color="white" if evidence[i,j] else GREY)
    axs[1].set_title("Evidence dimensions", fontweight="bold"); panel(axs[1], "B")
    axs[2].axis("off"); panel(axs[2], "C")
    counts = [
        ("4", "independent bulk cohorts"),
        (f"{coex_summary['high_confidence_core_abs_r_ge_0.40']}", "high-confidence coexpressed genes"),
        (f"{coex_summary['stable_mechanism_module']}", "stable mechanism-module genes"),
        ("6/6", "GSE210248 leave-one-patient-out\ndirection concordance in SMC1/SMC2"),
    ]
    for i, (n, lab) in enumerate(counts):
        yy = 0.82 - i * 0.21
        axs[2].text(0.08, yy, n, color=ORANGE, fontweight="bold", fontsize=18, va="center")
        axs[2].text(0.32, yy, lab, color=NAVY, fontsize=8.5, va="center")
    axs[2].text(0.05, 0.04, "Selection claim: convergent candidate evidence\n—not a causal target-validation result.",
                color=RED, fontsize=9, fontweight="bold")
    save(fig, "Figure_3_PDE8B_candidate_prioritization")


def fig4_localization():
    effects = pd.read_csv(OUTPUTS / "scrna-patient-pseudobulk-audit" /
                          "pde8b_patient_level_effect_summary.csv")
    bounds = pd.read_csv(OUTPUTS / "gse293580-reanalysis" / "sensitivity" /
                         "gse293580_annotation_boundary_tests.csv")
    q = bounds[(bounds.contrast == "PAH_combined_vs_Donor") &
               (bounds.metric == "PDE8B_log1p_CPM") &
               bounds.boundary.isin(["core_SMC", "strict_SMC_with_like", "contractile_SMC"])]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    title(fig, 4, "PDE8B localizes to SMC and replicates at the patient level")
    e = effects.iloc[:2].copy(); y = np.arange(2)[::-1]
    axs[0].axvline(0, color=GREY, lw=0.8); axs[0].scatter(e.hedges_g, y, s=75, c=[BLUE, ORANGE])
    axs[0].set_yticks(y, ["GSE210248 SMC1", "GSE210248 SMC2"])
    axs[0].set_xlim(-0.04, 1.12)
    axs[0].set_xlabel("Patient-level Hedges g"); axs[0].set_title("SMC localization", fontweight="bold")
    for yy, (_, row) in zip(y, e.iterrows()):
        axs[0].text(row.hedges_g + .06, yy, f"g={row.hedges_g:.2f}; LOO {int(row.loo_direction_concordant_n)}/{int(row.loo_iterations_n)}", va="center", fontsize=7)
    panel(axs[0], "A")
    labels = {"core_SMC": "Core SMC", "strict_SMC_with_like": "Strict SMC + like",
              "contractile_SMC": "Contractile SMC"}
    q = q.assign(label=q.boundary.map(labels)).sort_values("hedges_g")
    yy = np.arange(len(q)); axs[1].axvline(0, color=GREY, lw=0.8)
    axs[1].scatter(q.hedges_g, yy, s=65, c=[BLUE, ORANGE, TEAL])
    axs[1].set_yticks(yy, q.label); axs[1].set_xlabel("Hedges g (PAH vs donor)")
    axs[1].set_title("GSE293580 boundary sensitivity", fontweight="bold"); panel(axs[1], "B")
    for y0, (_, row) in zip(yy, q.iterrows()):
        axs[1].text(row.hedges_g + .12, y0, f"P={row.exact_permutation_p:.3f}", va="center", fontsize=7)
    axs[2].axis("off"); panel(axs[2], "C")
    box(axs[2], (0.10, 0.62), 0.75, 0.18, "Localization", "PDE8B signal is enriched\nin SMC compartments", edge=BLUE, fill=LIGHT_BLUE)
    box(axs[2], (0.10, 0.35), 0.75, 0.18, "Replication", "Positive patient-level effects\nin two datasets", edge=TEAL, fill="#EAF7F4")
    box(axs[2], (0.10, 0.08), 0.75, 0.18, "Boundary", "Small n; exact P values shown;\nno cell-level pseudoreplication", edge=RED, fill="#FBEDED")
    save(fig, "Figure_4_single_cell_localization_patient_validation")


def fig5_state_programs():
    g210 = pd.read_csv(OUTPUTS / "smc-composition-within-state" /
                       "gse210248_smc_composition_within_state_effects.csv")
    g293 = pd.read_csv(OUTPUTS / "gse293580-reanalysis" /
                       "gse293580_smc_program_group_tests.csv")
    get210 = lambda ct, feature: float(g210[(g210.layer == "Within-state") & (g210.cell_type == ct) & (g210.feature == feature)].hedges_g.iloc[0])
    combined = g293[g293.contrast == "PAH_combined_vs_Donor"].set_index("module")
    mat = np.array([
        [get210("SMC 1", "PDE8B"), get210("SMC 2", "PDE8B"), 2.137],
        [get210("SMC 1", "Contractile"), get210("SMC 2", "Contractile"), combined.loc["contractile", "hedges_g"]],
        [get210("SMC 1", "Synthetic/ECM"), get210("SMC 2", "Synthetic/ECM"), combined.loc["matrix_remodeling", "hedges_g"]],
    ])
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.55), gridspec_kw={"width_ratios": [1.12, 1.00, 1.42]})
    title(fig, 5, "PDE8B elevation accompanies contractile loss and remodeling-state gain")
    norm = TwoSlopeNorm(vmin=-2.2, vcenter=0, vmax=2.2)
    im = axs[0].imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")
    axs[0].set_yticks(range(3), ["PDE8B", "Contractile", "ECM/remodeling"])
    axs[0].set_xticks(range(3), ["GSE210\nSMC1", "GSE210\nSMC2", "GSE293580\nstrict SMC"])
    for i in range(3):
        for j in range(3): axs[0].text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", color="white" if abs(mat[i,j]) > 1.1 else NAVY)
    axs[0].set_title("Cross-dataset direction matrix", fontweight="bold"); panel(axs[0], "A")
    comp = g210[(g210.layer == "Composition") & g210.feature.isin(["SMC 1", "SMC 2", "SMC2 fraction within SMC"])]
    axs[1].barh(np.arange(len(comp)), comp.hedges_g, color=[BLUE, ORANGE, TEAL])
    short_labels = {"SMC 1": "SMC1", "SMC 2": "SMC2", "SMC2 fraction within SMC": "SMC2 / total SMC"}
    axs[1].axvline(0, color=GREY, lw=.8); axs[1].set_yticks(np.arange(len(comp)), [short_labels.get(x, x) for x in comp.feature])
    axs[1].set_xlabel("Hedges g"); axs[1].set_title("Composition is analyzed separately", fontweight="bold"); panel(axs[1], "B")
    axs[2].axis("off"); panel(axs[2], "C")
    box(axs[2], (0.05, .57), .36, .25, "Contractile\nSMC", "", edge=BLUE, fill=LIGHT_BLUE)
    box(axs[2], (.59, .57), .36, .25, "Remodeling\nSMC", "", edge=ORANGE, fill="#FFF3E7")
    axs[2].add_patch(FancyArrowPatch((.42, .695), (.58, .695), arrowstyle="-|>", mutation_scale=16, color=ORANGE, lw=2))
    axs[2].text(.5, .47, "PAH-associated co-occurrence", ha="center", color=NAVY, fontweight="bold")
    axs[2].text(.5, .32, "PDE8B ↑  ·  contractile ↓  ·  ECM/remodeling ↑", ha="center", color=RED, fontsize=7.2)
    axs[2].text(.5, .12, "Cross-sectional association; lineage trajectory\nand PDE8B-driven causality remain untested.", ha="center", color=GREY, fontsize=6.5)
    fig.subplots_adjust(top=.80, wspace=.62)
    save_figure(fig, OUT, "Figure_5_SMC_state_program_remodeling")


def fig6_interactions():
    lr = pd.read_csv(OUTPUTS / "gse293580-reanalysis" / "ligand-receptor-proxy" /
                     "gse293580_mural_lr_patient_tests.csv")
    q = lr[(lr.smc_boundary == "core_SMC") & (lr.contrast == "SSc-PAH_vs_Donor")]
    axes = ["ECM_integrin", "TGFb", "PDGF", "JAG_NOTCH3"]
    q = q[q.axis.isin(axes)].copy()
    fig, axs = plt.subplots(1, 3, figsize=(13, 5), gridspec_kw={"width_ratios": [1.2, 1, 1.05]})
    title(fig, 6, "Exploratory SMC–pericyte communication places remodeling in a mural-cell niche")
    pivot = q.pivot(index="axis", columns="direction", values="hedges_g").reindex(axes)
    x = np.arange(len(pivot)); width=.34
    direction_labels = {"SMC_to_pericyte": "SMC → pericyte", "pericyte_to_SMC": "Pericyte → SMC"}
    for i, (direction, color) in enumerate([("SMC_to_pericyte", ORANGE), ("pericyte_to_SMC", TEAL)]):
        axs[0].bar(x + (i-.5)*width, pivot[direction], width, color=color, label=direction_labels[direction])
    axs[0].axhline(0, color=GREY, lw=.8); axs[0].set_xticks(x, [a.replace("_", "–") for a in pivot.index], rotation=25, ha="right")
    axs[0].set_ylabel("Patient-level Hedges g"); axs[0].legend(frameon=False, fontsize=7)
    axs[0].set_title("SSc-PAH vs donor", fontweight="bold"); panel(axs[0], "A")
    ecm = q[q.axis == "ECM_integrin"].set_index("direction")
    vals = [ecm.loc["SMC_to_pericyte", "hedges_g"], ecm.loc["pericyte_to_SMC", "hedges_g"]]
    axs[1].barh([1, 0], vals, color=[ORANGE, TEAL])
    axs[1].set_yticks([1, 0], ["SMC → pericyte", "Pericyte → SMC"])
    axs[1].set_xlabel("Hedges g"); axs[1].set_title("ECM–integrin axis", fontweight="bold"); panel(axs[1], "B")
    for yy, direction in zip([1,0], ["SMC_to_pericyte", "pericyte_to_SMC"]):
        row=ecm.loc[direction]; axs[1].text(row.hedges_g+.12, yy, f"P={row.exact_permutation_p:.3f}\nq={row.BH_q_within_boundary_direction_contrast:.3f}", va="center", fontsize=7)
    axs[2].axis("off"); panel(axs[2], "C")
    box(axs[2], (.04,.63), .34,.20,"SMC","PDE8B ↑\nremodeling program",edge=ORANGE,fill="#FFF3E7")
    box(axs[2], (.62,.63), .34,.20,"Pericyte","mural-cell\nresponse compartment",edge=TEAL,fill="#EAF7F4")
    axs[2].add_patch(FancyArrowPatch((.39,.75),(.61,.75),arrowstyle="-|>",mutation_scale=15,color=ORANGE,lw=2,linestyle="--"))
    axs[2].add_patch(FancyArrowPatch((.61,.67),(.39,.67),arrowstyle="-|>",mutation_scale=15,color=TEAL,lw=2,linestyle="--"))
    axs[2].text(.5,.50,"ECM–integrin / TGFβ / PDGF",ha="center",fontweight="bold",color=NAVY)
    axs[2].text(.5,.34,"Exploratory patient-level association",ha="center",color=RED,fontweight="bold")
    axs[2].text(.5,.17,"Small sample; family-wise BH correction is not significant.\nDashed arrows denote a testable model, not proven signaling.",ha="center",color=GREY,fontsize=7.5)
    save(fig, "Figure_6_SMC_pericyte_interaction_model")


def main():
    fig1_design(); fig2_core_signature(); fig3_prioritization()
    fig4_localization(); fig5_state_programs(); fig6_interactions()
    manifest = pd.DataFrame({"figure": sorted(p.stem for p in OUT.glob("*.png")),
                             "status": "editorial draft v1"})
    manifest.to_csv(OUT / "six_figure_manifest.csv", index=False)
    print("\n".join(str(p) for p in sorted(OUT.glob("*.png"))))


if __name__ == "__main__":
    main()
