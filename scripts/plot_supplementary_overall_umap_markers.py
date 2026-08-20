"""Supplementary overall UMAP plus cell-type marker dot plot for GSE210248."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import COLORS, panel_label, save_figure  # noqa: E402


ORDER = ["Fibro", "Endo 1", "Endo 2", "Mono / Macs", "T cells", "NK cells",
         "B cells", "DC", "Granulocytes", "Epithelial", "Mast cells", "SMC 1", "SMC 2"]

MARKER_BLOCKS = {
    "Fibro": ["COL1A1", "DCN"],
    "Endo 1": ["PECAM1", "VWF"],
    "Endo 2": ["CA4", "RGCC"],
    "Mono / Macs": ["LYZ", "C1QC"],
    "T cells": ["CD3D", "CD3E"],
    "NK cells": ["NKG7", "GNLY"],
    "B cells": ["CD79A", "MS4A1"],
    "DC": ["FCER1A", "CD1C"],
    "Granulocytes": ["S100A8", "S100A9"],
    "Epithelial": ["EPCAM", "KRT19"],
    "Mast cells": ["TPSAB1", "KIT"],
    "SMC 1": ["TAGLN", "MYH11"],
    "SMC 2": ["COL3A1", "FN1"],
}
MARKERS = [g for ct in ORDER for g in MARKER_BLOCKS[ct]]

PALETTE = {
    "Fibro":"#35A9D4", "Endo 1":"#19B76B", "Endo 2":"#B3A52A",
    "Mono / Macs":"#1FB77A", "T cells":"#E58A24", "NK cells":"#D84C8B",
    "B cells":"#E6C43B", "DC":"#B58BC7", "Granulocytes":"#836AAE",
    "Epithelial":"#E05A61", "Mast cells":"#8B624A", "SMC 1":"#9A78B5",
    "SMC 2":"#7465A7",
}

DOT_CMAP = LinearSegmentedColormap.from_list("marker_bpr", ["#3B6FB6", "#A08AC4", "#F06473"])


def aggregate_markers(base: Path, um: pd.DataFrame) -> pd.DataFrame:
    with gzip.open(base / "features.tsv.gz", "rt") as h:
        genes = [x.rstrip().split("\t")[1] for x in h]
    with gzip.open(base / "barcodes.tsv.gz", "rt") as h:
        barcodes = [x.strip() for x in h]
    if barcodes != um.barcode.tolist():
        raise ValueError("Barcode order differs between matrix and UMAP table")
    gene_to_idx = {g:i for i,g in enumerate(genes)}
    present = [g for g in MARKERS if g in gene_to_idx]
    X = mmread(str(base / "matrix.mtx.gz")).tocsr().T[:, [gene_to_idx[g] for g in present]]
    lib = np.asarray(mmread(str(base / "matrix.mtx.gz")).tocsr().sum(axis=0)).ravel()
    rows=[]
    for ct in ORDER:
        idx=np.where(um.cell_type.to_numpy()==ct)[0]
        z=X[idx]
        mean=np.asarray(z.multiply(np.divide(1e4,lib[idx],out=np.zeros(len(idx)),where=lib[idx]>0)[:,None]).mean(axis=0)).ravel()
        det=np.asarray((z>0).mean(axis=0)).ravel()
        for g,m,p in zip(present,np.log1p(mean),det): rows.append({"cell_type":ct,"gene":g,"average_log_cp10k":m,"fraction_expressed":p})
    return pd.DataFrame(rows)


def main(data_root: Path, outdir: Path):
    um=pd.read_csv(data_root/"processed/gse210248_harmony/gse210248_raw_vs_harmony_umap.csv.gz")
    base=data_root/"raw/current-data-snapshot/pah_scrna/gse210248/Human_PA_Integrated"
    dot=aggregate_markers(base,um); outdir.mkdir(parents=True,exist_ok=True)
    # Gene-wise scaling makes each marker's annotation specificity visible.
    dot["scaled_expression"] = dot.groupby("gene")["average_log_cp10k"].transform(
        lambda x: (x-x.min())/(x.max()-x.min()) if x.max()>x.min() else 0.0
    )
    dot.to_csv(outdir/"Supplementary_overall_marker_dotplot_source_data.csv",index=False)
    fig=plt.figure(figsize=(8.2,4.55)); gs=fig.add_gridspec(1,2,left=.05,right=.94,top=.97,bottom=.20,width_ratios=[1.12,1.72],wspace=.24)
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1])
    for ct in ORDER:
        z=um[um.cell_type.eq(ct)]; a.scatter(z.harmony_UMAP1,z.harmony_UMAP2,s=.8,color=PALETTE[ct],alpha=.65,linewidths=0,rasterized=True,label=ct)
    a.set(xticks=[],yticks=[],xlabel="UMAP 1",ylabel="UMAP 2"); panel_label(a,"a")
    a.legend(loc="upper center",bbox_to_anchor=(.50,-.10),ncol=4,fontsize=5.2,markerscale=3,handletextpad=.25,columnspacing=.65)
    genes=[g for g in MARKERS if g in dot.gene.unique()]
    starts={}; cursor=0
    for ct in ORDER:
        block=[g for g in MARKER_BLOCKS[ct] if g in genes]; starts[ct]=(cursor,cursor+len(block)); cursor+=len(block)
    for k,ct in enumerate(ORDER):
        lo,hi=starts[ct]
        if k%2==1: b.axvspan(lo-.5,hi-.5,color="#E9ECEE",zorder=0)
    for i,ct in enumerate(ORDER):
        z=dot[dot.cell_type.eq(ct)].set_index("gene").reindex(genes)
        x=np.arange(len(z)); b.scatter(x,np.full(len(z),i),s=5+105*z.fraction_expressed,c=z.scaled_expression,cmap=DOT_CMAP,vmin=0,vmax=1,edgecolor="white",linewidth=.20,zorder=2)
    b.set_xticks(range(len(genes)),genes,rotation=90,ha="center",fontsize=5.2); b.set_yticks(range(len(ORDER)),ORDER,fontsize=6.2); b.invert_yaxis(); panel_label(b,"b")
    for ct,(lo,hi) in starts.items():
        b.text((lo+hi-1)/2,-.85,ct.replace("Mono / Macs","Mono/Macs"),ha="center",va="bottom",fontsize=4.8,rotation=35,color=PALETTE[ct])
    # Compact size legend and color scale, matching the reference layout.
    handles=[plt.scatter([],[],s=5+105*x,color="#8D7DB9",edgecolor="none") for x in [.10,.30,.50]]
    leg=b.legend(handles,["10%","30%","50%"],title="Expressed",loc="upper left",bbox_to_anchor=(1.01,1.0),fontsize=5.2,title_fontsize=5.5,frameon=False)
    b.add_artist(leg)
    sm=plt.cm.ScalarMappable(cmap=DOT_CMAP,norm=plt.Normalize(0,1)); cb=fig.colorbar(sm,ax=b,fraction=.025,pad=.02); cb.set_label("Relative expression",fontsize=5.5); cb.set_ticks([0,1]); cb.set_ticklabels(["Low","High"]); cb.ax.tick_params(labelsize=5)
    save_figure(fig,outdir,"Supplementary_Figure_overall_UMAP_and_markers")


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--data-root",type=Path,required=True); p.add_argument("--outdir",type=Path,required=True); z=p.parse_args(); main(z.data_root,z.outdir)
