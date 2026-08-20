"""Build Figure 1 study workflow and Table 1 cohort overview.

Only audited cohort-level fields are reported. Unavailable demographic and
subtype fields are not inferred. Figure 2-5 outputs are not touched.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402

OUTPUTS = PROJECT_DATA_ROOT / "outputs"
OUT = OUTPUTS / "figure1-publication-v4"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#20262D"
MUTED = "#66727B"
LINE = "#B7C0C6"
PALE = "#EFF2F4"
BLUE = "#6F98B5"
TEAL = "#438E86"
ORANGE = "#C9773D"
RED = "#B94747"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.2,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

COHORTS = [
    {"name": "GSE15197", "role": "Discovery", "platform": "Microarray (GPL6480)",
     "tissue": "Lung", "pah": 18, "control": 13, "total": 31,
     "notes": "PAH vs control; IPF-PH samples used separately for specificity"},
    {"name": "GSE113439", "role": "Discovery", "platform": "Microarray (GPL6244)",
     "tissue": "Lung", "pah": 14, "control": 11, "total": 25,
     "notes": "One CTEPH sample excluded"},
    {"name": "GSE254617", "role": "Discovery", "platform": "RNA-seq",
     "tissue": "Lung", "pah": 82, "control": 52, "total": 134,
     "notes": "Technical columns aggregated to 134 unique PHBI subjects"},
    {"name": "GSE208592", "role": "Discovery", "platform": "RNA-seq",
     "tissue": "Lung", "pah": 15, "control": 18, "total": 33,
     "notes": "PAH vs control lung samples"},
    {"name": "GSE53408", "role": "Independent validation", "platform": "Microarray raw CEL (GPL6244)",
     "tissue": "Lung", "pah": 12, "control": 11, "total": 23,
     "notes": "Raw CEL files reprocessed with RMA; excluded from primary meta-analysis"},
]


def box(ax, xy, wh, title, body, face=PALE, edge="white", title_color=INK,
        title_size=7.4, body_size=6.0):
    x, y = xy; w, h = wh
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.012",
                                facecolor=face, edgecolor=edge, linewidth=.9))
    ax.text(x+w/2, y+h*.70, title, ha="center", va="center", fontweight="bold",
            fontsize=title_size, color=title_color)
    ax.text(x+w/2, y+h*.33, body, ha="center", va="center", fontsize=body_size,
            color=INK, linespacing=1.24)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9,
                                 linewidth=1.0, color=LINE, shrinkA=3, shrinkB=3))


def save_figure(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor="white")


def draw_workflow():
    fig, ax = plt.subplots(figsize=(7.2, 6.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    # Stage 1: bulk discovery and independent validation.
    ax.add_patch(Rectangle((.035,.89),.93,.038,facecolor="#D9DDE0",edgecolor="none"))
    ax.text(.50,.909,"1  Bulk-lung discovery and validation",ha="center",va="center",
            fontsize=8.2,fontweight="bold")
    box(ax,(.04,.71),(.16,.135),"Discovery cohorts","4 lung cohorts\n223 subjects",face="#E5EBEF")
    box(ax,(.235,.71),(.16,.135),"Meta-analysis","15,165 genes\n251 FDR < 0.05",face="#DDE8ED")
    box(ax,(.43,.71),(.16,.135),"Robust set","188 genes\n4-cohort concordance",face="#A9C8C4")
    box(ax,(.625,.71),(.16,.135),"Raw-CEL validation","GSE53408\n23 subjects",face="#DDE8ED")
    box(ax,(.82,.71),(.145,.135),"Focal candidates","PDE8B\nPIEZO2 | SLC16A12",face="#A9C8C4")
    for x0,x1 in [(.20,.235),(.395,.43),(.59,.625),(.785,.82)]: arrow(ax,(x0,.777),(x1,.777))

    # Stage 2: two analytically distinct branches.
    ax.add_patch(Rectangle((.035,.635),.93,.038,facecolor="#D9DDE0",edgecolor="none"))
    ax.text(.50,.654,"2  Candidate-focused localization and program analysis",ha="center",va="center",
            fontsize=8.2,fontweight="bold")
    box(ax,(.055,.445),(.27,.135),"scRNA-seq discovery",
        "GSE210248\ncell localization\npatient pseudobulk",
        face="#E6F0EE",title_color=TEAL,title_size=7.0,body_size=5.8)
    box(ax,(.365,.445),(.27,.135),"scRNA-seq replication",
        "GSE293580\nindependent cohort\npatient-level evidence",
        face="#D5E8E4",title_color=TEAL,title_size=7.0,body_size=5.8)
    box(ax,(.675,.445),(.27,.135),"PDE8B-centered bulk program",
        "Four-cohort coexpression\nLOCO stability\nPathway / PAH-gene context",
        face="#F0DDCC",title_color=ORANGE,title_size=7.0,body_size=5.8)
    arrow(ax,(.325,.512),(.365,.512))
    arrow(ax,(.635,.512),(.675,.512))

    # Stage 3: synthesis, with a visible non-causal boundary.
    ax.add_patch(Rectangle((.035,.365),.93,.038,facecolor="#D9DDE0",edgecolor="none"))
    ax.text(.50,.384,"3  Evidence synthesis and mechanistic interpretation",ha="center",va="center",
            fontsize=8.2,fontweight="bold")
    box(ax,(.055,.13),(.27,.16),"Integrated evidence",
        "Cross-cohort association\nTwo-cohort SMC localization",face="#E5EBEF")
    box(ax,(.365,.13),(.27,.16),"Supported interpretation",
        "SMC-associated PDE8B signal\nLinked state remodeling",
        face="#D9E5E2",title_color=RED)
    box(ax,(.675,.13),(.27,.16),"Working hypothesis",
        "PDE8B-linked cAMP remodeling\nin pulmonary vascular SMCs",
        face="#F0DDCC",title_color=ORANGE)
    arrow(ax,(.325,.21),(.365,.21))
    arrow(ax,(.635,.21),(.675,.21))
    save_figure(fig,"Figure_1_overall_project_workflow_v5_draft")
    plt.close(fig)


def write_table_csv():
    path = OUT / "Table_1_bulk_lung_cohort_overview.csv"
    fields = ["Cohort", "Analysis role", "Platform", "Tissue", "PAH subjects",
              "Control subjects", "Total subjects", "Key inclusion/exclusion note"]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for c in COHORTS:
            writer.writerow({"Cohort":c["name"],"Analysis role":c["role"],"Platform":c["platform"],
                             "Tissue":c["tissue"],"PAH subjects":c["pah"],
                             "Control subjects":c["control"],"Total subjects":c["total"],
                             "Key inclusion/exclusion note":c["notes"]})
        writer.writerow({"Cohort":"Total","Analysis role":"All bulk cohorts","Platform":"—","Tissue":"Lung",
                         "PAH subjects":141,"Control subjects":105,"Total subjects":246,
                         "Key inclusion/exclusion note":"Discovery: 223 subjects; independent validation: 23 subjects"})


def draw_table():
    columns=["Characteristic"]+[c["name"] for c in COHORTS]+["Total"]
    rows=[
        ["Analysis role"]+[c["role"].replace("Independent validation","Independent\nvalidation") for c in COHORTS]+["—"],
        ["Platform"]+[c["platform"].replace(" (GPL6480)","").replace(" (GPL6244)","").replace("Microarray raw CEL","Microarray\nraw CEL") for c in COHORTS]+["—"],
        ["Tissue"]+[c["tissue"] for c in COHORTS]+["Lung"],
        ["PAH subjects"]+[str(c["pah"]) for c in COHORTS]+["141"],
        ["Control subjects"]+[str(c["control"]) for c in COHORTS]+["105"],
        ["Total subjects"]+[str(c["total"]) for c in COHORTS]+["246"],
        ["Primary meta-analysis"]+["Yes","Yes","Yes","Yes","No","4 cohorts"],
    ]
    fig,ax=plt.subplots(figsize=(7.2,3.15)); ax.axis("off")
    table=ax.table(cellText=rows,colLabels=columns,cellLoc="center",colLoc="center",
                   loc="upper left",bbox=[0,.18,1,.76],colWidths=[.19,.132,.132,.132,.132,.15,.132])
    table.auto_set_font_size(False); table.set_fontsize(6.3)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor("#3B4247"); cell.set_linewidth(.6)
        if r==0:
            cell.set_facecolor("#D9DDE0"); cell.set_text_props(fontweight="bold")
        elif c==0:
            cell.set_facecolor("#F0F2F3"); cell.set_text_props(ha="left")
        elif r%2==0: cell.set_facecolor("#F7F8F8")
    save_figure(fig,"Table_1_bulk_lung_cohort_overview_v1_draft")
    plt.close(fig)


def main():
    draw_workflow(); write_table_csv(); draw_table()


if __name__ == "__main__":
    main()
