"""Build manuscript Figure 2: transparent candidate convergence and validation.

The main figure compares all three prespecified candidates. Cohort-dependence
(LOCO) is exported separately as a supplementary figure. No statistical model
is rerun here; this script visualizes validated tabular outputs.
"""

from __future__ import annotations

import sys
from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402

AUDIT = OUTPUT_ROOT / "pah_audit"
OUTPUTS = PROJECT_DATA_ROOT / "outputs"
OUT = OUTPUTS / "manuscript-figures-v3"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]
COHORT_LABELS = ["15197", "113439", "254617", "208592"]
COHORT_DETAIL = ["Array\n18 / 13", "Array\n14 / 11", "RNA-seq\n82 / 52", "RNA-seq\n15 / 18"]
GENES = ["PDE8B", "PIEZO2", "SLC16A12"]
COLORS = {"PDE8B": "#B84343", "PIEZO2": "#477AA8", "SLC16A12": "#187B7B"}
RED, BLUE, TEAL, GREY, DARK = "#B84343", "#477AA8", "#187B7B", "#AEB8C2", "#263238"

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def panel(ax, label: str, x: float) -> None:
    """Place a panel letter at the author-adjusted F2 panel boundary."""
    ax.text(x, 1.12, label.upper(), transform=ax.transAxes, fontsize=10,
            fontweight=700, va="bottom", ha="left", clip_on=False)


def save_all(fig, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor="white")


def load_data():
    meta_all = pd.read_csv(AUDIT / "lung_mechanism_random_effects_meta.csv")
    robust = meta_all[(meta_all.same_direction.astype(bool)) & (meta_all.fdr < .05) &
                      (meta_all.I2 < 50) & (meta_all.meta_g.abs() >= .50)].copy()
    meta = robust.set_index("symbol").loc[GENES]
    comp = pd.read_csv(OUTPUTS / "candidate-comparison" / "three_candidate_evidence_matrix.csv").set_index("gene").loc[GENES]
    raw = pd.read_csv(AUDIT / "gse53408_raw_cel_candidate_validation.csv").set_index("symbol").loc[GENES]
    raw_all = pd.read_csv(PROJECT_DATA_ROOT / "interim" / "gse53408_rma" /
                          "gse53408_full_differential_expression.csv.gz")
    expr = pd.read_csv(PROJECT_DATA_ROOT / "interim" / "gse53408_rma" / "gse53408_rma_gene_expression.csv.gz")
    loo = pd.read_csv(OUTPUTS / "gse53408-sensitivity" / "gse53408_candidate_leave_one_sample_out.csv")
    loco = pd.read_csv(OUTPUTS / "bulk-meta-robustness" / "lung_meta_loco_summary.csv")
    return robust, meta, comp, raw, raw_all, expr, loo, loco


def draw_main(robust, meta, comp, raw, raw_all, expr, loo):
    fig = plt.figure(figsize=(7.2, 7.3))
    gs = fig.add_gridspec(3, 2, left=.09, right=.98, top=.91, bottom=.10,
                          height_ratios=[1.35, 1.0, 1.05], width_ratios=[1.08, 1],
                          hspace=.62, wspace=.58)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, :]); d = fig.add_subplot(gs[2, :])

    # a: all 188 robust genes, sorted by pooled effect for reproducible display.
    z = robust.sort_values("meta_g", ascending=False).reset_index(drop=True)
    gene_mat = z[[f"g_{x}" for x in COHORTS]].to_numpy(float)
    im = a.imshow(gene_mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2),
                  aspect="auto", interpolation="nearest", rasterized=True)
    a.set_xticks(range(4), [f"{g}\n{detail}" for g, detail in zip(COHORT_LABELS, COHORT_DETAIL)])
    a.set_yticks([0, len(z)-1], ["Higher in PAH", "Lower in PAH"])
    a.set_ylabel("188 robust genes")
    a.yaxis.set_label_coords(-.020, .50)
    cb=fig.colorbar(im,ax=a,fraction=.038,pad=.025); cb.set_label("Hedges g")
    a.set_title("The 188-gene robust set is consistent across four cohorts",
                loc="left", fontweight="bold", fontsize=8.4)
    panel(a, "a", -.316)

    # b: complete validation landscape for the robust set, not a funnel.
    val = robust[["symbol","meta_g"]].merge(
        raw_all[["symbol","hedges_g","adj.P.Val"]], on="symbol", how="left")
    # Two of 188 discovery genes lack a mapped raw-CEL gene-level estimate;
    # the displayed count is therefore 186 and is reported on-panel.
    shown = val[val["hedges_g"].notna()].copy()
    sig = shown["adj.P.Val"].lt(.05)
    b.scatter(shown.loc[~sig,"meta_g"],shown.loc[~sig,"hedges_g"],s=9,color="#C8D0D5",alpha=.65)
    b.scatter(shown.loc[sig,"meta_g"],shown.loc[sig,"hedges_g"],s=14,color="#79A99F",alpha=.78)
    lim=max(abs(shown[["meta_g","hedges_g"]].to_numpy()).max(),1.5)
    b.plot([-lim,lim],[-lim,lim],ls="--",lw=.8,color="#89949B")
    b.axhline(0,color="#B4BCC1",lw=.7); b.axvline(0,color="#B4BCC1",lw=.7)
    label_offsets={"PDE8B":(6,8),"PIEZO2":(6,-12),"SLC16A12":(6,5)}
    for gene in GENES:
        q=shown[shown.symbol.eq(gene)].iloc[0]
        b.scatter(q.meta_g,q.hedges_g,s=46,color=COLORS[gene],edgecolor="white",lw=.6,zorder=4)
        b.annotate(gene,(q.meta_g,q.hedges_g),xytext=label_offsets[gene],textcoords="offset points",
                   fontsize=6.5,fontweight="bold",color=COLORS[gene])
    b.set(xlabel="Four-cohort meta Hedges g",ylabel="GSE53408 raw-CEL Hedges g")
    b.set_title("Independent-cohort validation of the robust set",
                loc="left",fontweight="bold",fontsize=8.4)
    panel(b,"b", -.178)

    # c: return the three candidates to the four discovery cohorts.
    mat = meta[[f"g_{x}" for x in COHORTS]].to_numpy(float)
    im = c.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5), aspect="auto")
    c.set_yticks(range(3), GENES, fontweight="bold")
    c.set_xticks(range(4), [f"{g}\n{detail}" for g, detail in zip(COHORT_LABELS, COHORT_DETAIL)])
    for i in range(3):
        for j in range(4):
            c.text(j,i,f"{mat[i,j]:+.2f}",ha="center",va="center",fontsize=7,fontweight="bold",
                   color="white" if abs(mat[i,j])>.72 else DARK)
    cb=fig.colorbar(im,ax=c,fraction=.024,pad=.02); cb.set_label("Hedges g")
    c.set_title("All three validated candidates retain consistent effects in each discovery cohort",
                loc="left", fontweight="bold"); panel(c,"c", -.126)

    # d: retain sample-level RNA expression for the three focal candidates.
    sub=gs[2,:].subgridspec(1,3,wspace=.42)
    fig.delaxes(d); axes=[fig.add_subplot(sub[0,i]) for i in range(3)]
    gsm_cols=[x for x in expr.columns if x!="symbol"]
    groups=np.array(["PAH" if int(x.replace("GSM",""))<=1290999 else "Control" for x in gsm_cols])
    for k,(ax,gene) in enumerate(zip(axes,GENES)):
        vals=expr.loc[expr.symbol.eq(gene),gsm_cols].iloc[0].astype(float).to_numpy()
        for xpos,grp,col in [(0,"Control",BLUE),(1,"PAH",RED)]:
            q=vals[groups==grp]; jitter=np.linspace(-.08,.08,len(q))
            ax.scatter(xpos+jitter,q,s=20,color=col,alpha=.78,edgecolor="white",lw=.35)
            ax.hlines(q.mean(),xpos-.18,xpos+.18,color=DARK,lw=1.3)
        ax.set_xticks([0,1],["Control\n(n=11)","PAH\n(n=12)"])
        ax.set_title(gene,color=COLORS[gene],fontweight="bold")
        annotation_y = .905 if gene == "SLC16A12" else .96
        ax.text(.04,annotation_y,f"g={raw.loc[gene,'hedges_g']:+.2f}\nFDR={raw.loc[gene,'adj.P.Val']:.3f}",
                transform=ax.transAxes,va="top",fontsize=6.2)
        if k==0: ax.set_ylabel("RMA expression"); panel(ax,"d", -.465)
        else: ax.set_yticklabels([])

    save_all(fig, "Figure_2_gene_set_and_independent_validation_v6_draft")
    plt.close(fig)


def draw_loco_supplement(meta, loco):
    pde = meta.loc["PDE8B"]
    fig, ax = plt.subplots(figsize=(4.8, 3.1))
    fig.subplots_adjust(left=.22, right=.97, top=.84, bottom=.25)
    y = np.arange(len(loco))[::-1]
    ax.errorbar(loco.PDE8B_meta_g, y,
                xerr=[loco.PDE8B_meta_g-loco.PDE8B_ci_low,
                      loco.PDE8B_ci_high-loco.PDE8B_meta_g],
                fmt="o", color=TEAL, ecolor="#7D8992", capsize=2.5, markersize=5)
    ax.axvline(0, color="#7D8992", lw=.8)
    ax.axvline(pde.meta_g, color=RED, lw=1, ls="--", label="Full four-cohort meta")
    ax.set_yticks(y, [f"Omit {x}" for x in loco.omitted_cohort])
    ax.set_xlabel("Three-cohort random-effects Hedges g")
    ax.set_title("Supplementary Figure | PDE8B leave-one-cohort-out analysis",
                 loc="left", fontweight="bold")
    ax.legend(loc="lower right")
    save_all(fig, "Supplementary_Figure_PDE8B_LOCO_v4_draft")
    plt.close(fig)


def main():
    robust, meta, comp, raw, raw_all, expr, loo, loco = load_data()
    draw_main(robust, meta, comp, raw, raw_all, expr, loo)
    draw_loco_supplement(meta, loco)


if __name__ == "__main__":
    main()
