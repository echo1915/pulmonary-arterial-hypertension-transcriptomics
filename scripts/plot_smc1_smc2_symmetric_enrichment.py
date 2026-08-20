from pathlib import Path
import argparse, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_style import COLORS, apply_publication_style, save_figure

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans","sans-serif"],
 "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.spines.top":False,"axes.spines.right":False})
apply_publication_style()
COL={"SMC1":COLORS["smc1"],"SMC2":COLORS["smc2"]}

def select_terms(d, side):
    q=d[d.side.eq(side)].copy()
    patterns = (["muscle system process","muscle contraction","actin filament-based process","actin cytoskeleton organization",
                 "Smooth Muscle Contraction","oxidative phosphorylation","mitochondrial ATP synthesis"] if side=="SMC1" else
                ["Extracellular matrix organization","extracellular structure organization","collagen fibril organization",
                 "ECM proteoglycans","Collagen biosynthesis","vasculature development","regulation of locomotion"])
    rows=[]; used=set()
    for p in patterns:
        z=q[q.term.str.contains(re.escape(p),case=False,regex=True)].sort_values("fdr")
        if len(z):
            r=z.iloc[0]
            if r.term.lower() not in used: rows.append(r); used.add(r.term.lower())
    return pd.DataFrame(rows).head(6)

def main(inp,out):
    out.mkdir(parents=True,exist_ok=True); d=pd.read_csv(inp)
    a=select_terms(d,"SMC1"); b=select_terms(d,"SMC2")
    n=max(len(a),len(b)); a=a.reset_index(drop=True); b=b.reset_index(drop=True)
    rows=[]
    for i in range(n):
        if i<len(a): rows.append(a.iloc[i].to_dict())
        if i<len(b): rows.append(b.iloc[i].to_dict())
    pd.DataFrame(rows).to_csv(out/"Figure_SMC1_SMC2_symmetric_enrichment_source_data.csv",index=False)
    fig,ax=plt.subplots(figsize=(7.2,3.45)); y=np.arange(n)[::-1]
    def side_plot(z,sign):
        x=sign*(-np.log10(z.fdr.clip(lower=1e-50)).clip(upper=20).to_numpy())
        s=28+260*z.gene_ratio.to_numpy()
        ax.scatter(x,y[:len(z)],s=s,c=COL[z.side.iloc[0]],alpha=.88,edgecolor="white",linewidth=.7,zorder=3)
        for xx,yy,t in zip(x,y,z.term):
            label=t.replace("mitochondrial ATP synthesis coupled electron transport","Mitochondrial ATP synthesis")
            ax.text(xx+(.35 if sign>0 else -.35),yy,label,ha="left" if sign>0 else "right",va="center",fontsize=6.6)
    side_plot(a,-1); side_plot(b,1)
    ax.axvline(0,color="#687680",lw=.9); ax.set_yticks([])
    lim=max(21,ax.get_xlim()[1],-ax.get_xlim()[0]); ax.set_xlim(-lim,lim)
    ticks=np.arange(-20,21,5); ax.set_xticks(ticks,[str(abs(x)) for x in ticks])
    ax.set_xlabel("−log10(FDR)")
    ax.text(.02,1.03,"SMC1-enriched programs",transform=ax.transAxes,color=COL["SMC1"],fontsize=9,fontweight="bold",ha="left")
    ax.text(.98,1.03,"SMC2-enriched programs",transform=ax.transAxes,color=COL["SMC2"],fontsize=9,fontweight="bold",ha="right")
    ax.set_title("Patient-level enrichment separates contractile and matrix-remodeling SMC programs",pad=17)
    for ratio,x in zip([.06,.12,.18],[.43,.49,.55]):
        ax.scatter(x,.05,s=28+260*ratio,transform=ax.transAxes,c="#8C99A3",edgecolor="white",linewidth=.6)
        ax.text(x,.012,f"{ratio:.0%}",transform=ax.transAxes,ha="center",fontsize=5.8)
    ax.text(.49,.10,"Gene ratio",transform=ax.transAxes,ha="center",fontsize=6.2,color="#5D6870")
    fig.text(.08,.02,"Top 300 genes per side from paired SMC2−SMC1 patient pseudobulk ranking; g:Profiler GO:BP/Reactome FDR. n=6 patients.",fontsize=6.2,color="#56636C")
    fig.subplots_adjust(left=.06,right=.98,top=.82,bottom=.18)
    stem=out/"Figure_SMC1_SMC2_symmetric_enrichment_v1"
    save_figure(fig, out, stem.name)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True); p.add_argument('--out',type=Path,required=True); z=p.parse_args(); main(z.input,z.out)
