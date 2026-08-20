from pathlib import Path
import argparse
import gzip
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figure_style import COLORS, apply_publication_style, panel_label, save_figure

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "legend.frameon": False,
})

apply_publication_style()
COL = {"SMC 1": COLORS["smc1"], "SMC 2": COLORS["smc2"],
       "Donor": COLORS["muted"], "PAH": COLORS["signal"]}
MODULE_ORDER = ["Contractile", "Synthetic/ECM", "Mural remodeling", "PDE8B"]
MODULE_LABEL = {"Contractile":"Contractile", "Synthetic/ECM":"Synthetic / ECM",
                "Mural remodeling":"Mural remodeling", "PDE8B":"PDE8B"}

def save_all(fig, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "Figure_SMC1_SMC2_state_continuum_v1"
    save_figure(fig, outdir, stem.name)

def main(data_root: Path, outdir: Path):
    umap_file = data_root / "processed/gse210248_harmony/gse210248_raw_vs_harmony_umap.csv.gz"
    score_file = data_root / "outputs/smc-composition-within-state/gse210248_smc_within_state_module_scores.csv"
    comp_file = data_root / "outputs/smc-composition-within-state/gse210248_smc_composition_patient_values.csv"
    effect_file = data_root / "outputs/smc-composition-within-state/gse210248_smc_composition_within_state_effects.csv"
    um = pd.read_csv(umap_file)
    um = um[um["cell_type"].isin(["SMC 1", "SMC 2"])].copy()
    sc = pd.read_csv(score_file)
    comp = pd.read_csv(comp_file)
    effects = pd.read_csv(effect_file)

    # Patient-level phenotype index: higher = less contractile and more ECM/remodeling.
    wide = sc.pivot_table(index=["sample","disease","cell_type"], columns="module", values="score").reset_index()
    wide["remodeling_index"] = (-wide["Contractile"] + wide["Synthetic/ECM"] + wide["Mural remodeling"]) / 3
    wide.to_csv(outdir / "Figure_SMC1_SMC2_state_continuum_source_data.csv", index=False)

    fig = plt.figure(figsize=(7.20, 5.75), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1.05, 1],
                          left=.075, right=.98, top=.90, bottom=.12, wspace=.30, hspace=.38)
    axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1])
    axC = fig.add_subplot(gs[1,0]); axD = fig.add_subplot(gs[1,1])

    # A: cell-level localization only.
    for ct in ["SMC 1", "SMC 2"]:
        z = um[um.cell_type.eq(ct)]
        axA.scatter(z.harmony_UMAP1, z.harmony_UMAP2, s=2.2, alpha=.42,
                    c=COL[ct], linewidths=0, rasterized=True, label=f"{ct} (cells={len(z):,})")
        axA.text(z.harmony_UMAP1.median(), z.harmony_UMAP2.median(), ct,
                 color=COL[ct], fontsize=8, fontweight="bold", ha="center",
                 bbox=dict(boxstyle="round,pad=.18", fc="white", ec="none", alpha=.82))
    axA.set(xlabel="Harmony UMAP 1", ylabel="Harmony UMAP 2", title="SMC states occupy adjacent transcriptional space")
    axA.legend(loc="upper right", fontsize=6, markerscale=2.4)
    panel_label(axA, "a")

    # B: paired patient-level state scores, preserving biological replicate.
    modules = MODULE_ORDER[:3]
    x = np.arange(len(modules))
    offsets = {"SMC 1": -.11, "SMC 2": .11}
    for i, mod in enumerate(modules):
        q = wide[["sample","disease","cell_type",mod]].pivot_table(index=["sample","disease"], columns="cell_type", values=mod).reset_index()
        for _, r in q.iterrows():
            axB.plot([i+offsets["SMC 1"], i+offsets["SMC 2"]], [r["SMC 1"], r["SMC 2"]],
                     color=COL[r.disease], alpha=.45, lw=.8, zorder=1)
            for ct in ["SMC 1","SMC 2"]:
                axB.scatter(i+offsets[ct], r[ct], s=22, facecolor=COL[ct],
                            edgecolor=COL[r.disease], linewidth=1.0, zorder=2)
    axB.axhline(0, color="#B7C1C8", lw=.7, ls="--")
    axB.set_xticks(x, [MODULE_LABEL[m] for m in modules], rotation=18, ha="right")
    axB.set_ylabel("Patient pseudobulk z-score")
    axB.set_title("Patient-level programs distinguish SMC1 and SMC2")
    axB.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COL[k],label=k,markersize=5) for k in ["SMC 1","SMC 2"]], loc="upper right", fontsize=6)
    panel_label(axB, "b")

    # C: disease effects for each state; biological replicate is the patient.
    y = np.arange(len(MODULE_ORDER))[::-1]
    for j, mod in enumerate(MODULE_ORDER):
        for ct in ["SMC 1","SMC 2"]:
            row = effects[(effects.layer.eq("Within-state")) & effects.cell_type.eq(ct) & effects.feature.eq(mod)].iloc[0]
            mean = row.hedges_g
            yy = y[j] + (-.10 if ct=="SMC 1" else .10)
            axC.scatter(mean, yy, s=34, color=COL[ct], zorder=3)
            axC.plot([0,mean],[yy,yy],color=COL[ct],lw=1.3,alpha=.8)
    axC.axvline(0,color="#9AA7B0",lw=.8,ls="--")
    axC.set_yticks(y,[MODULE_LABEL[m] for m in MODULE_ORDER])
    axC.set_xlabel("Hedges g: PAH versus donor\n(3 PAH, 3 donors)")
    axC.set_title("Disease effects reveal contractile loss and PDE8B gain")
    axC.legend(handles=[Line2D([0],[0],marker='o',color=COL[k],label=k,markersize=5) for k in ["SMC 1","SMC 2"]], loc="lower right", fontsize=6)
    panel_label(axC, "c")

    # D: patient-level PDE8B vs remodeling index; correlations are exploratory.
    for ct, marker in [("SMC 1","o"),("SMC 2","s")]:
        z = wide[wide.cell_type.eq(ct)]
        for disease in ["Donor","PAH"]:
            q = z[z.disease.eq(disease)]
            axD.scatter(q.remodeling_index, q.PDE8B, s=38, marker=marker,
                        facecolor=COL[disease], edgecolor="white", linewidth=.7,
                        alpha=.9, label=f"{ct}, {disease}")
    rx = wide.remodeling_index.rank(method="average").to_numpy()
    ry = wide.PDE8B.rank(method="average").to_numpy()
    rho = np.corrcoef(rx, ry)[0, 1]
    axD.axhline(0,color="#C0C8CE",lw=.7); axD.axvline(0,color="#C0C8CE",lw=.7)
    axD.text(.04,.96,f"Spearman ρ = {rho:.2f}\nn = 12 state–patient pseudobulks\n(exploratory; no P value)",
             transform=axD.transAxes,ha="left",va="top",fontsize=6.4)
    axD.set(xlabel="Remodeling index\n(−contractile + ECM + remodeling)/3", ylabel="PDE8B standardized score",
            title="PDE8B is not tightly coupled to remodeling state")
    axD.legend(fontsize=5.6, ncol=2, loc="lower right", handletextpad=.3, columnspacing=.7)
    panel_label(axD, "d")

    fig.suptitle("SMC1 and SMC2 define a patient-resolved phenotypic continuum in PAH lung",
                 x=.075, ha="left", fontsize=11, fontweight="bold")
    fig.text(.075,.018,"Cells are used for visualization only; quantitative panels use patient-level pseudobulk values (3 donors, 3 PAH).",
             fontsize=6.5, color="#4A5560")
    save_all(fig, outdir)
    plt.close(fig)

if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args=ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True); main(args.data_root,args.outdir)
