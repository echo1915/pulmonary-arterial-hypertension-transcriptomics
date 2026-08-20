from __future__ import annotations

from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


O = PROJECT_DATA_ROOT / "outputs"
P = PROJECT_DATA_ROOT / "processed" / "gse210248_harmony" / "gse210248_raw_vs_harmony_umap.csv.gz"
TREND = O / "smc-gene-state-trends-review-v1" / "Figure_SMC_gene_state_trends_source_data.csv"
MARKERS = O / "pde8b-coexpression" / "state-deconvolution" / "gse210248_smc_marker_pseudobulk.csv"
ENRICH = O / "smc-symmetric-enrichment-v1" / "Figure_SMC1_SMC2_symmetric_enrichment_source_data.csv"
SCORES = O / "smc-state-continuum-v1" / "Figure_SMC1_SMC2_state_continuum_source_data.csv"
OUT = O / "figure4-smc-main-supp-review-v1"

SMC1, SMC2 = "#E86B6B", "#42B883"
BLUE, PURPLE, RED, INK, GRID = "#4169A1", "#7652A6", "#F21A24", "#263746", "#E2E6E9"
DOT_CMAP = LinearSegmentedColormap.from_list("dot", [BLUE, PURPLE, RED])

mpl.rcParams.update({
    "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7, "axes.spines.top":False, "axes.spines.right":False,
    "axes.linewidth":.75, "svg.fonttype":"none", "pdf.fonttype":42,
    "legend.frameon":False,
})


def lab(ax, s, x_fig=None):
    fig=ax.figure
    y_fig=fig.transFigure.inverted().transform(ax.transAxes.transform((0,1.10)))[1]
    if x_fig is None:
        x_fig=.015 if ax.get_position().x0<.50 else .55
    fig.text(x_fig,y_fig,s.upper(),fontsize=10,fontweight=700,va="bottom",ha="left")


def shift_axes(ax, dx_points):
    """Reproduce the author's horizontal panel adjustment in physical points."""
    p=ax.get_position()
    dx=dx_points/(72*ax.figure.get_figwidth())
    ax.set_position([p.x0+dx,p.y0,p.width,p.height])


def smooth(x, y, bins=26):
    edges=np.unique(np.quantile(x,np.linspace(0,1,bins+1))); xx=[]; yy=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(x>=lo)&(x<hi if hi!=edges[-1] else x<=hi)
        if m.sum()>=10: xx.append(np.median(x[m])); yy.append(np.mean(y[m]))
    return np.asarray(xx),np.asarray(yy)


def save(fig, stem):
    fig.savefig(stem.with_suffix(".png"),dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"),dpi=600,bbox_inches="tight",facecolor="white",
                pil_kwargs={"compression":"tiff_lzw"})


def plot_umap(ax, um, focused=False):
    for ct,c in [("SMC 1",SMC1),("SMC 2",SMC2)]:
        q=um[um.cell_type.eq(ct)]
        ax.scatter(q.harmony_UMAP1,q.harmony_UMAP2,s=3.2,alpha=.42,c=c,lw=0,rasterized=True,label=ct)
    if focused:
        xlo,xhi=um.harmony_UMAP1.quantile([.01,.99]); ylo,yhi=um.harmony_UMAP2.quantile([.01,.99])
        ax.set_xlim(xlo,xhi); ax.set_ylim(ylo,yhi)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(loc="best",markerscale=2.2)


def plot_marker_dot(ax, marker):
    genes=["MYH11","TAGLN","ACTA2","CNN1","COL1A1","COL3A1","FN1","VIM"]
    q=marker[marker.gene.isin(genes)].copy()
    agg=q.groupby(["cell_type","gene"]).agg(avg=("log_cp10k","mean"), frac=("log_cp10k",lambda x:(x>0).mean())).reset_index()
    agg["scaled"]=agg.groupby("gene")["avg"].transform(lambda x:(x-x.min())/(x.max()-x.min()) if x.max()>x.min() else .5)
    for yi,ct in enumerate(["SMC 1","SMC 2"]):
        z=agg[agg.cell_type.eq(ct)].set_index("gene").reindex(genes)
        ax.scatter(np.arange(len(genes)),np.full(len(genes),yi),s=30+190*z.frac,
                   c=z.scaled,cmap=DOT_CMAP,vmin=0,vmax=1,edgecolor="white",lw=.6)
    ax.axvspan(-.5,3.5,color="#F3F4F5",zorder=-2); ax.axvspan(3.5,7.5,color="#E8EBED",zorder=-2)
    ax.set_xticks(range(len(genes)),genes,rotation=90); ax.set_yticks([0,1],["SMC 1","SMC 2"])
    ax.set_xlim(-.6,7.6); ax.set_ylim(-.65,1.65)
    ax.text(1.5,1.12,"Contractile",transform=ax.get_xaxis_transform(),ha="center",color=SMC1,fontsize=6.5)
    ax.text(5.5,1.12,"ECM",transform=ax.get_xaxis_transform(),ha="center",color=SMC2,fontsize=6.5)


def plot_trend(ax, trend, gene):
    for ct,c in [("SMC 1",SMC1),("SMC 2",SMC2)]:
        q=trend[trend.cell_type.eq(ct)]
        ax.scatter(q.state_score,q[gene],s=4.2,alpha=.22,c=c,lw=0,rasterized=True,label=ct)
    q=trend.sort_values("state_score"); xx,yy=smooth(q.state_score.to_numpy(),q[gene].to_numpy())
    ax.plot(xx,yy,c="#D64A61",lw=1.9)
    ax.set_title(gene,fontstyle="italic",pad=5); ax.set_xlabel("SMC state score\nECM − contractile")
    ax.set_ylabel("Expression (log CP10K)")


def plot_enrichment(ax,enrich):
    q=enrich.copy(); q["sig"]=-np.log10(q.fdr.clip(lower=1e-30)); q["x"]=np.where(q.side.eq("SMC1"),-q.sig,q.sig)
    q=q.sort_values(["side","sig"],ascending=[True,False])
    y=np.arange(len(q)); cols=np.where(q.side.eq("SMC1"),BLUE,"#8B5FBF")
    ax.barh(y,q.x,color=cols,height=.68,alpha=.94)
    short = {
        "Extracellular matrix organization": "ECM organization",
        "extracellular structure organization": "ECM structure organization",
        "Collagen biosynthesis and modifying enzymes": "Collagen biosynthesis",
        "actin filament-based process": "Actin filament process",
    }
    labels=[short.get(term,term) for term in q.term]
    ax.axvline(0,c=INK,lw=.8); ax.set_yticks(y,labels,fontsize=5.35); ax.invert_yaxis()
    ax.set_xlabel("SMC1-associated  ←  −log10(FDR)  →  SMC2-associated")
    ax.grid(axis="x",color=GRID,lw=.55,zorder=0)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    um=pd.read_csv(P); um=um[um.cell_type.isin(["SMC 1","SMC 2"])].copy()
    trend=pd.read_csv(TREND); marker=pd.read_csv(MARKERS); enrich=pd.read_csv(ENRICH); scores=pd.read_csv(SCORES)

    # Main Figure 4: compact, near-square panels. Marker dot plots are omitted;
    # patient-level program evidence is promoted from the supplement.
    fig=plt.figure(figsize=(7.45,8.2)); gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.15],hspace=.66,wspace=.62,left=.12,right=.97,top=.93,bottom=.08)
    a=fig.add_subplot(gs[0,0]); plot_umap(a,um,True); lab(a,"a")
    b=fig.add_subplot(gs[0,1]); shift_axes(b,-31)
    for i,mod in enumerate(["Contractile","Synthetic/ECM","Mural remodeling"]):
        q=scores[scores.cell_type.isin(["SMC 1","SMC 2"])][["sample","disease","cell_type",mod]].pivot(index=["sample","disease"],columns="cell_type",values=mod).reset_index()
        for _,r in q.iterrows():
            b.plot([i-.13,i+.13],[r["SMC 1"],r["SMC 2"]],c="#AAB2B8",alpha=.55,lw=.7)
            b.scatter([i-.13,i+.13],[r["SMC 1"],r["SMC 2"]],c=[SMC1,SMC2],s=18,zorder=3)
    b.axhline(0,c=INK,lw=.7); b.set_xticks(range(3),["Contractile","ECM","Remodeling"],rotation=20,ha="right"); b.set_ylabel("Patient pseudobulk z score"); lab(b,"b",.475)
    c=fig.add_subplot(gs[1,0]); plot_trend(c,trend,"MYH11"); lab(c,"c")
    d=fig.add_subplot(gs[1,1]); shift_axes(d,-31); plot_trend(d,trend,"FN1"); lab(d,"d",.475)
    e=fig.add_subplot(gs[2,0]); shift_axes(e,7)
    pah=scores[scores.disease.eq("PAH")]; wide=pah.pivot(index="sample",columns="cell_type",values=["Contractile","Synthetic/ECM","Mural remodeling"])
    names=["Mural remodeling","Contractile","Synthetic/ECM"]; vals=np.array([(wide[(n,"SMC 2")]-wide[(n,"SMC 1")]).mean() for n in names])
    e.barh(range(3),vals,color=[BLUE,PURPLE,"#C51B8A"]); e.axvline(0,c=INK,lw=.8); e.set_yticks(range(3),names); e.invert_yaxis(); e.set_xlabel("Mean paired SMC2 − SMC1 score in PAH"); lab(e,"e")
    f=fig.add_subplot(gs[2,1]); shift_axes(f,-7); plot_enrichment(f,enrich); lab(f,"f",.475)
    save(fig,OUT/"Figure_4_SMC_state_continuum_main_review_v1"); plt.close(fig)

    # Supplement: complete UMAP and the PDE8B state-specificity boundary only.
    fig=plt.figure(figsize=(6.7,3.25)); gs=fig.add_gridspec(1,2,wspace=.48,left=.10,right=.96,top=.84,bottom=.18)
    a=fig.add_subplot(gs[0,0]); plot_umap(a,um,False); lab(a,"a")
    b=fig.add_subplot(gs[0,1]); plot_trend(b,trend,"PDE8B"); lab(b,"b",.527)
    save(fig,OUT/"Supplementary_Figure_SMC_state_support_review_v1"); plt.close(fig)

    pd.DataFrame({"display":"main","total_cells":[len(um)],"x_q_low":[um.harmony_UMAP1.quantile(.01)],"x_q_high":[um.harmony_UMAP1.quantile(.99)],"y_q_low":[um.harmony_UMAP2.quantile(.01)],"y_q_high":[um.harmony_UMAP2.quantile(.99)]}).to_csv(OUT/"Figure4_umap_display_audit.csv",index=False)
    print(OUT)


if __name__=="__main__": main()
