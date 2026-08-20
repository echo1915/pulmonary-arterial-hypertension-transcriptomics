"""Figure 1 v2: bulk discovery linked to single-cell localization and validation."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402

OUT = PROJECT_DATA_ROOT / "outputs" / "figure1-integrative-v2"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#17324D"; BLUE = "#3977A8"; TEAL = "#2A9D8F"
ORANGE = "#E58B2B"; RED = "#C84D4D"; GREY = "#707B84"
PALE_BLUE = "#EAF2F8"; PALE_TEAL = "#EAF7F4"; PALE_ORANGE = "#FFF2E4"
PALE_RED = "#FBEDED"; PALE_GREY = "#F1F3F4"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"
})


def rounded(ax, x, y, w, h, title, subtitle="", ec=BLUE, fc="white",
            title_size=9, subtitle_size=6.8, dashed=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                       ec=ec, fc=fc, lw=1.5, linestyle="--" if dashed else "-")
    ax.add_patch(p)
    ax.text(x+w/2, y+h*0.64, title, ha="center", va="center",
            fontsize=title_size, fontweight="bold", color=ec, linespacing=1.05)
    if subtitle:
        ax.text(x+w/2, y+h*0.28, subtitle, ha="center", va="center",
                fontsize=subtitle_size, color=GREY, linespacing=1.12)


def arrow(ax, start, end, color=NAVY, dashed=False, width=1.2):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10,
                                 lw=width, color=color,
                                 linestyle="--" if dashed else "-"))


def panel_a(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.00, 1.02, "A", fontsize=13, fontweight="bold", color=NAVY)
    ax.text(0.08, 1.01, "Cohorts and modality-specific roles", fontsize=10.5,
            fontweight="bold", color=NAVY)
    ax.text(0.02, .88, "BULK RNA-seq", fontsize=9.5, fontweight="bold", color=BLUE)
    cohorts = [("GSE15197", "18 PAH / 13 donor"), ("GSE113439", "14 / 11"),
               ("GSE254617", "82 / 52"), ("GSE208592", "15 / 18")]
    for i, (name, n) in enumerate(cohorts):
        rounded(ax, .03+i*.235, .69, .205, .13, name, n, BLUE, PALE_BLUE, 8.5, 6.5)
    arrow(ax, (.50, .68), (.50, .58), BLUE)
    rounded(ax, .25, .45, .50, .12, "Cross-cohort disease discovery",
            "subjects/cohorts are biological replicates", BLUE, PALE_BLUE)
    ax.text(.02, .34, "SINGLE-CELL RNA-seq", fontsize=9.5, fontweight="bold", color=ORANGE)
    rounded(ax, .05, .13, .37, .15, "GSE210248 pulmonary artery",
            "SMC1/SMC2 localization\n3 PAH + 3 donor", ORANGE, PALE_ORANGE)
    rounded(ax, .56, .13, .37, .15, "GSE293580 lung",
            "SMC boundary, state programs, mural interaction\nIPAH / SSc-PAH / donor", ORANGE, PALE_ORANGE)
    arrow(ax, (.42, .44), (.24, .29), ORANGE)
    arrow(ax, (.58, .44), (.74, .29), ORANGE)
    ax.text(.50, .035, "Bulk defines reproducibility; single cell resolves cellular origin and state context.",
            ha="center", fontsize=7.1, color=NAVY, fontweight="bold")


def panel_b(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.00, 1.02, "B", fontsize=13, fontweight="bold", color=NAVY)
    ax.text(0.08, 1.01, "Convergent prioritization of PDE8B", fontsize=10.5,
            fontweight="bold", color=NAVY)
    funnel = [
        (.07, .82, .86, "15,015 common genes", BLUE, PALE_BLUE),
        (.12, .68, .76, "2,413 meta-analysis FDR < 0.05", BLUE, PALE_BLUE),
        (.17, .54, .66, "817 direction- and heterogeneity-robust", TEAL, PALE_TEAL),
        (.22, .40, .56, "543 stable mechanism-module genes", TEAL, PALE_TEAL),
        (.27, .26, .46, "201 high-confidence core genes", ORANGE, PALE_ORANGE),
    ]
    for x, y, w, text, edge, fill in funnel:
        poly = Polygon([(x,y+.10),(x+w,y+.10),(x+w-.035,y),(x+.035,y)],
                       closed=True, ec=edge, fc=fill, lw=1.3)
        ax.add_patch(poly); ax.text(.5, y+.052, text, ha="center", va="center",
                                    fontsize=7.4, color=edge, fontweight="bold")
    arrow(ax, (.50,.25),(.50,.18), ORANGE)
    rounded(ax, .28, .04, .44, .13, "PDE8B",
            "bulk specificity + SMC localization + patient replication",
            RED, PALE_RED, 12, 6.4)
    ax.text(.02, .015, "Candidate selection is prespecified and comparative—not based on one volcano plot.",
            fontsize=6.5, color=GREY)


def panel_c(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.00, 1.02, "C", fontsize=13, fontweight="bold", color=NAVY)
    ax.text(0.08, 1.01, "Evidence ladder and inference boundary", fontsize=10.5,
            fontweight="bold", color=NAVY)
    steps = [
        (.05,.77,.57,.12,"1  Association","PDE8B increased across 4 independent lung cohorts",BLUE,PALE_BLUE,False),
        (.12,.61,.57,.12,"2  Localization","positive patient-level effects in SMC1/SMC2",TEAL,PALE_TEAL,False),
        (.19,.45,.57,.12,"3  State context","contractile ↓; ECM/remodeling ↑",ORANGE,PALE_ORANGE,False),
        (.26,.29,.57,.12,"4  Mural-cell niche","SMC–pericyte ECM–integrin signal; exploratory",RED,PALE_RED,True),
        (.33,.10,.57,.12,"5  Causal validation","perturbation and functional rescue not yet tested",GREY,PALE_GREY,True),
    ]
    for x,y,w,h,t,s,e,f,d in steps:
        rounded(ax,x,y,w,h,t,s,e,f,8.7,6.4,d)
    for i in range(len(steps)-1):
        x1,y1,w1,h1,*_=steps[i]; x2,y2,w2,h2,*_=steps[i+1]
        arrow(ax,(x1+w1*.68,y1-.005),(x2+w2*.45,y2+h2+.005),
              GREY, dashed=i>=2, width=1.0)
    ax.text(.93,.52,"SUPPORTED",rotation=90,ha="center",va="center",
            color=TEAL,fontweight="bold",fontsize=8)
    ax.plot([.91,.91],[.44,.88],color=TEAL,lw=3,solid_capstyle="round")
    ax.text(.93,.34,"EXPLORATORY",rotation=90,ha="center",va="center",
            color=RED,fontweight="bold",fontsize=7)
    ax.plot([.91,.91],[.28,.41],color=RED,lw=3,linestyle="--",solid_capstyle="round")
    ax.text(.93,.15,"PENDING",rotation=90,ha="center",va="center",
            color=GREY,fontweight="bold",fontsize=8)
    ax.plot([.91,.91],[.09,.23],color=GREY,lw=3,linestyle="--",solid_capstyle="round")


def main():
    fig = plt.figure(figsize=(13, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.36, .92, 1.02], wspace=.10,
                          left=.035, right=.985, bottom=.08, top=.83)
    panel_a(fig.add_subplot(gs[0,0])); panel_b(fig.add_subplot(gs[0,1])); panel_c(fig.add_subplot(gs[0,2]))
    fig.text(.035,.94,"Figure 1 | Integrative bulk-to-single-cell strategy identifies an SMC-associated PDE8B program",
             fontsize=15, fontweight="bold", color=NAVY)
    fig.text(.035,.885,
             "Four independent lung cohorts establish reproducibility; two single-cell datasets resolve cellular localization, state remodeling and an exploratory mural-cell niche.",
             fontsize=8.5, color=GREY)
    fig.text(.035,.018,
             "GSE117261 overlaps PHBI/GSE254617 and is not counted independently. GSE53408 remains secondary validation. Dashed elements denote exploratory or untested inference.",
             fontsize=6.8, color=GREY)
    stem = OUT / "Figure_1_integrative_bulk_single_cell_v2"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
