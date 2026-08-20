"""Separate SMC composition shifts from within-state expression changes in GSE210248."""

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
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


STATE = PROJECT_DATA_ROOT / "outputs" / "pde8b-coexpression" / "state-deconvolution"
OUT = PROJECT_DATA_ROOT / "outputs" / "smc-composition-within-state"
OUT.mkdir(parents=True, exist_ok=True)
PANELS = {
    "PDE8B": ["PDE8B"],
    "Contractile": ["ACTA2", "TAGLN", "MYH11", "CNN1", "MYL9", "TPM2", "CALD1"],
    "Synthetic/ECM": ["KLF4", "COL1A1", "COL3A1", "FN1", "VIM"],
    "Mural remodeling": ["RGS5", "CSPG4", "MCAM", "NOTCH3", "PDGFRB"],
}

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
    pooled = math.sqrt(((len(case) - 1) * np.var(case, ddof=1) + (len(control) - 1) * np.var(control, ddof=1)) / (len(case) + len(control) - 2))
    if not np.isfinite(pooled) or pooled == 0:
        return 0.0
    return (np.mean(case) - np.mean(control)) / pooled * (1 - 3 / (4 * (len(case) + len(control)) - 9))


def exact_p(case: np.ndarray, control: np.ndarray) -> float:
    values = np.r_[case, control]
    observed = abs(np.mean(case) - np.mean(control))
    stats = []
    for ids in itertools.combinations(range(len(values)), len(case)):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(ids)] = True
        stats.append(abs(values[mask].mean() - values[~mask].mean()))
    return float(np.mean(np.asarray(stats) >= observed - 1e-12))


def effect(values: pd.DataFrame, value_col: str, feature: str, layer: str, cell_type: str = "Composition") -> dict:
    case = values.loc[values.disease.eq("PAH"), value_col].to_numpy()
    control = values.loc[values.disease.eq("Donor"), value_col].to_numpy()
    return {
        "layer": layer, "cell_type": cell_type, "feature": feature,
        "pah_mean": float(case.mean()), "donor_mean": float(control.mean()),
        "delta": float(case.mean() - control.mean()), "hedges_g": hedges_g(case, control),
        "exact_permutation_p": exact_p(case, control), "n_PAH": len(case), "n_Donor": len(control),
    }


def main() -> None:
    prop = pd.read_csv(STATE / "gse210248_celltype_proportions.csv")
    smc = prop[prop.Cell_annotation.isin(["SMC 1", "SMC 2"])].copy()
    wide = smc.pivot_table(index=["sample", "disease"], columns="Cell_annotation", values="proportion", fill_value=0).reset_index()
    wide["Total SMC"] = wide["SMC 1"] + wide["SMC 2"]
    wide["SMC2 fraction within SMC"] = np.where(wide["Total SMC"] > 0, wide["SMC 2"] / wide["Total SMC"], np.nan)
    composition_rows = [effect(wide, feature, feature, "Composition") for feature in ["SMC 1", "SMC 2", "Total SMC", "SMC2 fraction within SMC"]]

    genes = pd.read_csv(STATE / "gse210248_smc_marker_pseudobulk.csv")
    genes["z"] = genes.groupby(["cell_type", "gene"])["log_cp10k"].transform(
        lambda v: (v - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else 0
    )
    score_records = []
    effect_rows = []
    for panel, members in PANELS.items():
        scores = (genes[genes.gene.isin(members)]
                  .groupby(["sample", "disease", "cell_type"], as_index=False)
                  .agg(score=("z", "mean"), genes_observed=("gene", "nunique")))
        scores["module"] = panel
        score_records.append(scores)
        for cell_type, block in scores.groupby("cell_type"):
            effect_rows.append(effect(block, "score", panel, "Within-state", cell_type))
    scores = pd.concat(score_records, ignore_index=True)
    effects = pd.DataFrame(composition_rows + effect_rows)

    wide.to_csv(OUT / "gse210248_smc_composition_patient_values.csv", index=False, encoding="utf-8-sig")
    scores.to_csv(OUT / "gse210248_smc_within_state_module_scores.csv", index=False, encoding="utf-8-sig")
    effects.to_csv(OUT / "gse210248_smc_composition_within_state_effects.csv", index=False, encoding="utf-8-sig")
    interpretation = {
        "primary": "PAH shows lower SMC1 abundance, modestly higher SMC2 abundance, and within-subtype loss of contractile programs.",
        "pde8b": "PDE8B is directionally higher within both SMC1 and SMC2, arguing against a purely compositional explanation.",
        "boundary": "All comparisons are descriptive (3 PAH vs 3 donors) and exact permutation P values are not multiplicity-adjusted.",
    }
    (OUT / "smc_composition_within_state_summary.json").write_text(json.dumps(interpretation, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), gridspec_kw={"width_ratios": [1.05, 1.25, 1.15]})
    comp = effects[(effects.layer == "Composition") & effects.feature.isin(["SMC 1", "SMC 2", "Total SMC", "SMC2 fraction within SMC"])].copy()
    y = np.arange(len(comp))
    axes[0].hlines(y, 0, comp.hedges_g, color="#CCD4DA")
    axes[0].scatter(comp.hedges_g, y, c=np.where(comp.hedges_g >= 0, "#D97732", "#477AA8"), s=42)
    axes[0].axvline(0, color="#68747E", lw=0.8)
    axes[0].set_yticks(y, comp.feature)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Hedges g (PAH vs donor)")
    axes[0].set_title("A  SMC composition", loc="left", fontweight="bold")

    within = effects[(effects.layer == "Within-state") & effects.feature.isin(["PDE8B", "Contractile", "Synthetic/ECM"])].copy()
    modules = ["PDE8B", "Contractile", "Synthetic/ECM"]
    x = np.arange(len(modules))
    for offset, cell_type, color in [(-0.16, "SMC 1", "#477AA8"), (0.16, "SMC 2", "#D97732")]:
        q = within[within.cell_type.eq(cell_type)].set_index("feature").loc[modules]
        axes[1].bar(x + offset, q.hedges_g, width=0.30, color=color, label=cell_type)
    axes[1].axhline(0, color="#68747E", lw=0.8)
    axes[1].set_xticks(x, modules, rotation=20, ha="right")
    axes[1].set_ylabel("Hedges g (PAH vs donor)")
    axes[1].set_title("B  Changes within each SMC subtype", loc="left", fontweight="bold")
    axes[1].legend(frameon=False)

    ax = axes[2]
    for xpos, disease in enumerate(["Donor", "PAH"]):
        q = wide[wide.disease.eq(disease)]
        jitter = np.linspace(-0.05, 0.05, len(q))
        ax.scatter(xpos + jitter, q["SMC2 fraction within SMC"], s=48,
                   color={"Donor": "#477AA8", "PAH": "#B84343"}[disease], edgecolor="white", lw=0.5)
        for xx, (_, row) in zip(xpos + jitter, q.iterrows()):
            ax.annotate(row["sample"], (xx, row["SMC2 fraction within SMC"]), xytext=(2, 3), textcoords="offset points", fontsize=6)
        ax.hlines(q["SMC2 fraction within SMC"].mean(), xpos - 0.18, xpos + 0.18, color="#263238", lw=1.5)
    state_row = comp[comp.feature.eq("SMC2 fraction within SMC")].iloc[0]
    ax.set_xticks([0, 1], ["Donor", "PAH"])
    ax.set_ylabel("SMC2 / (SMC1 + SMC2)")
    ax.set_title("C  Relative shift toward SMC2", loc="left", fontweight="bold")
    ax.text(0.03, 0.97, f"g={state_row.hedges_g:.2f}\nexact P={state_row.exact_permutation_p:.2f}", transform=ax.transAxes, va="top", fontsize=7)

    fig.suptitle("PDE8B elevation accompanies SMC state remodeling rather than a simple abundance shift", x=0.04, ha="left", fontsize=12.2, fontweight="bold")
    fig.text(0.04, 0.015, "Each estimate uses patient-level summaries (3 PAH and 3 donors). State-module values are standardized within gene and SMC subtype before averaging.", fontsize=7, color="#56616B")
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    fig.savefig(OUT / "Figure_SMC_composition_within_state_draft.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "Figure_SMC_composition_within_state_draft.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Figure_SMC_composition_within_state_draft.svg", bbox_inches="tight")
    plt.close(fig)
    print(effects.to_string(index=False))


if __name__ == "__main__":
    main()
