"""Build the data-dense five-figure manuscript architecture.

This script reorganizes existing validated outputs. It does not rerun statistical
analyses or change biological replicate definitions.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import COLORS, figure_title, panel_label, save_figure  # noqa: E402
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


ROOT = PROJECT_DATA_ROOT
O = ROOT / "outputs"
OUT = O / "manuscript-figures-dense-v2"
INK, MUTED = COLORS["ink"], COLORS["muted"]
BLUE, ORANGE = COLORS["smc1"], COLORS["smc2"]
TEAL, RED, GREY = COLORS["teal"], COLORS["signal"], COLORS["neutral"]
COHORTS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]

# Keep the submission-critical settings explicit in this entry-point as well as
# in figure_style.py so standalone source preflight can verify them.
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def finish(fig, stem):
    save_figure(fig, OUT, stem)


def effect_lollipop(ax, labels, values, colors=None, xlabel="Hedges g"):
    y = np.arange(len(labels))[::-1]
    colors = colors or [BLUE] * len(labels)
    ax.axvline(0, color=GREY, lw=.8)
    for yy, val, col in zip(y, values, colors):
        ax.plot([0, val], [yy, yy], color=col, lw=1.2, alpha=.75)
        ax.scatter(val, yy, s=28, color=col, edgecolor="white", linewidth=.5, zorder=3)
    ax.set_yticks(y, labels)
    ax.set_xlabel(xlabel)


def fig1_bulk():
    meta = pd.read_csv(O / "current-results/pah_audit/lung_mechanism_random_effects_meta.csv")
    loco = pd.read_csv(O / "bulk-meta-robustness/lung_meta_loco_summary.csv")
    pde = meta.loc[meta.symbol.eq("PDE8B")].iloc[0]
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, left=.09, right=.98, top=.91, bottom=.10,
                          width_ratios=[1.05, 1.35], height_ratios=[1, 1],
                          hspace=.42, wspace=.36)
    a, b, c, d = [fig.add_subplot(gs[i, j]) for i, j in [(0,0),(0,1),(1,0),(1,1)]]
    figure_title(fig, "Cross-cohort lung transcriptomics identifies a reproducible PAH program")

    names = COHORTS + ["Random effects"]
    vals = [pde[f"g_{x}"] for x in COHORTS] + [pde.meta_g]
    y = np.arange(5)[::-1]
    a.axvline(0, color=GREY, lw=.8)
    a.scatter(vals[:4], y[:4], s=35, color=BLUE)
    a.errorbar(pde.meta_g, y[-1], xerr=[[pde.meta_g-pde.ci_low], [pde.ci_high-pde.meta_g]],
               fmt="D", color=RED, capsize=2.5, ms=5)
    a.set_yticks(y, names); a.set_xlabel("Hedges g (PAH vs donor)")
    a.set_title("PDE8B is concordantly increased")
    a.text(.02, -.25, f"Meta g={pde.meta_g:.2f}; FDR={pde.fdr:.1e}; I²={pde.I2:.0f}%",
           transform=a.transAxes, color=RED, fontsize=6.5)
    panel_label(a, "a")

    m = meta.copy(); m["mlog"] = -np.log10(m.fdr.clip(lower=1e-300))
    robust = m.same_direction & (m.fdr < .05) & (m.I2 < 50)
    b.scatter(m.loc[~robust,"meta_g"], m.loc[~robust,"mlog"], s=3, color="#D5DBDE", alpha=.45, rasterized=True)
    b.scatter(m.loc[robust,"meta_g"], m.loc[robust,"mlog"], s=7, color=TEAL, alpha=.72, rasterized=True)
    b.scatter(pde.meta_g, -np.log10(pde.fdr), s=45, color=RED, zorder=4)
    b.annotate("PDE8B", (pde.meta_g, -np.log10(pde.fdr)), xytext=(5,5), textcoords="offset points", color=RED, fontweight="bold")
    b.axvline(0, color=GREY, lw=.7); b.set(xlabel="Random-effects Hedges g", ylabel="−log10(FDR)")
    b.set_title("Transcriptome-wide random-effects meta-analysis"); panel_label(b,"b")

    top = pd.concat([m.nsmallest(8, "meta_g"), m.nlargest(8, "meta_g")]).drop_duplicates("symbol")
    top = top.sort_values("meta_g")
    mat = top[[f"g_{x}" for x in COHORTS]].to_numpy()
    im = c.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2), aspect="auto")
    c.set_yticks(range(len(top)), top.symbol); c.set_xticks(range(4), [x.replace("GSE","") for x in COHORTS], rotation=35, ha="right")
    c.set_title("Leading concordant effects across cohorts"); panel_label(c,"c")
    fig.colorbar(im, ax=c, shrink=.65, pad=.02, label="Hedges g")

    d.errorbar(loco.PDE8B_meta_g, np.arange(len(loco)),
               xerr=[loco.PDE8B_meta_g-loco.PDE8B_ci_low, loco.PDE8B_ci_high-loco.PDE8B_meta_g],
               fmt="o", color=BLUE, ecolor=BLUE, capsize=2.5, ms=4)
    d.axvline(0, color=GREY, lw=.8); d.axvline(pde.meta_g, color=RED, lw=.9, ls="--")
    d.set_yticks(np.arange(len(loco)), [f"omit {x}" for x in loco.omitted_cohort])
    d.set_xlabel("PDE8B leave-one-cohort-out g"); d.set_title("Signal persists after omitting each cohort")
    panel_label(d,"d")
    finish(fig, "Figure_1_bulk_discovery_and_robustness")


def fig2_candidate():
    comp = pd.read_csv(O / "candidate-comparison/three_candidate_evidence_matrix.csv")
    loo = pd.read_csv(O / "gse53408-sensitivity/gse53408_candidate_leave_one_sample_out.csv")
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.7), gridspec_kw={"height_ratios":[1,1.05]})
    fig.subplots_adjust(left=.10, right=.98, top=.90, bottom=.10, wspace=.34, hspace=.44)
    figure_title(fig, "Multi-layer validation prioritizes PDE8B over alternative lung candidates")
    a,b,c,d = axs.ravel()

    metrics = ["main_meta_g","PAH_vs_IPF_PH_g","GSE53408_raw_g","GSE210248_best_SMC_g","GSE293580_IPAH_g","GSE293580_SSc_PAH_g"]
    labs = ["Primary meta","PAH vs IPF-PH","GSE53408","GSE210248 SMC","GSE293580 IPAH","GSE293580 SSc-PAH"]
    x=np.arange(len(labs)); width=.23
    gene_cols={"PDE8B":RED,"PIEZO2":BLUE,"SLC16A12":TEAL}
    for j,(_,r) in enumerate(comp.iterrows()):
        a.scatter(x+(j-1)*width, [r[z] for z in metrics], s=28, color=gene_cols[r.gene], label=r.gene)
    a.axhline(0,color=GREY,lw=.8); a.set_xticks(x,labs,rotation=35,ha="right"); a.set_ylabel("Hedges g")
    a.legend(ncol=3, loc="lower left", bbox_to_anchor=(0,1.01), handletextpad=.3, columnspacing=.8)
    a.set_title("Parallel candidate performance", pad=28); panel_label(a,"a")

    p=comp.set_index("gene").loc["PDE8B"]
    vals=[p.main_meta_g,p.PAH_vs_IPF_PH_g,p.GSE53408_raw_g,p.GSE210248_best_SMC_g,p.GSE293580_IPAH_g,p.GSE293580_SSc_PAH_g]
    effect_lollipop(b,labs,vals,[RED,ORANGE,TEAL,BLUE,BLUE,ORANGE])
    b.set_title("PDE8B evidence is directionally coherent"); panel_label(b,"b")

    if "symbol" in loo.columns: q=loo[loo.symbol.eq("PDE8B")]
    elif "gene" in loo.columns: q=loo[loo.gene.eq("PDE8B")]
    else: q=loo
    eff_col = "hedges_g" if "hedges_g" in q.columns else next(z for z in q.columns if z.endswith("_g") or "effect" in z.lower())
    c.scatter(np.arange(len(q)), q[eff_col], s=25, color=TEAL)
    c.axhline(0,color=GREY,lw=.8); c.set_xlabel("Leave-one-sample-out iteration"); c.set_ylabel("PDE8B Hedges g")
    c.set_title("GSE53408 secondary validation is not sample-driven"); panel_label(c,"c")

    evidence=np.array([[1,1,1,1,1],[1,1,1,0,1],[1,1,1,1,0]],float)
    d.imshow(evidence,cmap=matplotlib.colors.ListedColormap(["#E5E9EB",TEAL]),vmin=0,vmax=1,aspect="auto")
    d.set_yticks(range(3),comp.gene); d.set_xticks(range(5),["4-cohort\ndirection","Disease\nspecificity","Secondary\nbulk","SMC\nlocalization","Etiologic\nreplication"],rotation=25)
    for i in range(3):
        for j in range(5): d.text(j,i,"●" if evidence[i,j] else "–",ha="center",va="center",color="white" if evidence[i,j] else MUTED)
    d.set_title("Convergent evidence dimensions"); panel_label(d,"d")
    finish(fig,"Figure_2_candidate_prioritization_and_validation")


def fig3_single_cell():
    u=pd.read_csv(ROOT/"processed/gse210248_umap/gse210248_recomputed_umap.csv.gz")
    h=pd.read_csv(ROOT/"processed/gse210248_harmony/gse210248_raw_vs_harmony_umap.csv.gz",usecols=["barcode","harmony_UMAP1","harmony_UMAP2"])
    u=u.merge(h,on="barcode",how="left"); vals=pd.read_csv(O/"scrna-patient-pseudobulk-audit/gse210248_pde8b_smc_patient_values.csv")
    prop=pd.read_csv(O/"pde8b-coexpression/state-deconvolution/gse210248_celltype_proportions.csv")
    fig=plt.figure(figsize=(7.2,6.1)); gs=fig.add_gridspec(2,3,left=.07,right=.98,top=.91,bottom=.10,hspace=.42,wspace=.36,height_ratios=[1.1,1])
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); c=fig.add_subplot(gs[0,2]); d=fig.add_subplot(gs[1,:2]); e=fig.add_subplot(gs[1,2])
    figure_title(fig,"Single-cell data localize the PDE8B signal to PAH-associated SMC states")
    cats=sorted(u.cell_type.unique()); cmap=plt.cm.tab20(np.linspace(0,1,len(cats))); cc=dict(zip(cats,cmap))
    for ct in cats:
        q=u[u.cell_type.eq(ct)]; a.scatter(q.harmony_UMAP1,q.harmony_UMAP2,s=.8,color=cc[ct],alpha=.55,rasterized=True)
        if len(q)>100: a.text(q.harmony_UMAP1.median(),q.harmony_UMAP2.median(),ct,fontsize=5.0,ha="center")
    a.set(xticks=[],yticks=[],xlabel="Harmony UMAP 1",ylabel="UMAP 2"); a.set_title("Pulmonary-artery cell landscape"); panel_label(a,"a")
    zero=u.PDE8B_log_cp10k<=0; b.scatter(u.loc[zero,"harmony_UMAP1"],u.loc[zero,"harmony_UMAP2"],s=.7,color="#D7DCDF",alpha=.35,rasterized=True)
    pos=u[~zero].sort_values("PDE8B_log_cp10k"); sc=b.scatter(pos.harmony_UMAP1,pos.harmony_UMAP2,c=pos.PDE8B_log_cp10k,s=2.2,cmap="magma",rasterized=True)
    b.set(xticks=[],yticks=[],xlabel="Harmony UMAP 1"); b.set_title("PDE8B expression"); panel_label(b,"b"); fig.colorbar(sc,ax=b,shrink=.55,pad=.02)
    focus=prop[prop.Cell_annotation.isin(["SMC 1","SMC 2","Endo 1","Endo 2","Fibro"])]
    order=["SMC 1","SMC 2","Endo 1","Endo 2","Fibro"]; x=np.arange(len(order)); w=.34
    for i,(grp,col) in enumerate([("Donor",MUTED),("PAH",RED)]):
        z=focus[focus.disease.eq(grp)].groupby("Cell_annotation").proportion.mean().reindex(order)
        c.bar(x+(i-.5)*w,z,w,color=col,label=grp)
    c.set_xticks(x,[z.replace(" ","\n") for z in order]); c.set_ylabel("Cell fraction"); c.legend(); c.set_title("Structural-cell composition"); panel_label(c,"c")
    for i,ct in enumerate(["SMC 1","SMC 2"]):
        q=vals[vals.cell_type.eq(ct)]
        for grp,col,off in [("Donor",MUTED,-.12),("PAH",RED,.12)]:
            z=q[q.disease.eq(grp)]; d.scatter(np.full(len(z),i+off),z.log_cp10k,s=34,color=col,edgecolor="white",linewidth=.5)
    d.set_xticks([0,1],["SMC1","SMC2"]); d.set_ylabel("PDE8B log(CP10K+1)"); d.set_title("Patient-level SMC pseudobulk expression"); panel_label(d,"d")
    det=vals.groupby(["cell_type","disease"]).detection_fraction.mean().unstack()
    e.imshow(det[["Donor","PAH"]],cmap="Blues",aspect="auto")
    e.set_yticks(range(len(det)),det.index); e.set_xticks([0,1],["Donor","PAH"])
    for i in range(len(det)):
        for j in range(2): e.text(j,i,f"{det.iloc[i,j]:.1%}",ha="center",va="center",fontsize=6)
    e.set_title("Detection fraction"); panel_label(e,"e")
    fig.text(.07,.025,"Cells are visualization units; expression comparisons use patient-level pseudobulk (3 PAH, 3 donors).",fontsize=6.3,color=MUTED)
    finish(fig,"Figure_3_single_cell_landscape_and_patient_replication")


def fig4_smc():
    um=pd.read_csv(ROOT/"processed/gse210248_harmony/gse210248_raw_vs_harmony_umap.csv.gz")
    um=um[um.cell_type.isin(["SMC 1","SMC 2"])]
    scores=pd.read_csv(O/"smc-composition-within-state/gse210248_smc_within_state_module_scores.csv")
    effects=pd.read_csv(O/"smc-composition-within-state/gse210248_smc_composition_within_state_effects.csv")
    enr=pd.read_csv(O/"smc-symmetric-enrichment-v1/Figure_SMC1_SMC2_symmetric_enrichment_source_data.csv")
    fig=plt.figure(figsize=(7.2,6.2)); gs=fig.add_gridspec(2,3,left=.08,right=.98,top=.91,bottom=.11,hspace=.43,wspace=.40,width_ratios=[1.15,1,1])
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1:]); c=fig.add_subplot(gs[1,0]); d=fig.add_subplot(gs[1,1:])
    figure_title(fig,"SMC1 and SMC2 form a contractile-to-matrix-remodeling state continuum")
    for ct,col in [("SMC 1",BLUE),("SMC 2",ORANGE)]:
        q=um[um.cell_type.eq(ct)]; a.scatter(q.harmony_UMAP1,q.harmony_UMAP2,s=1.7,color=col,alpha=.45,rasterized=True,label=ct)
    a.set(xticks=[],yticks=[],xlabel="Harmony UMAP 1",ylabel="UMAP 2"); a.legend(); a.set_title("Adjacent SMC transcriptional states"); panel_label(a,"a")
    modules=["Contractile","Synthetic/ECM","Mural remodeling"]
    wide=scores.pivot_table(index=["sample","disease","cell_type"],columns="module",values="score").reset_index()
    for i,mod in enumerate(modules):
        q=wide.pivot_table(index=["sample","disease"],columns="cell_type",values=mod).reset_index()
        for _,r in q.iterrows():
            b.plot([i-.10,i+.10],[r["SMC 1"],r["SMC 2"]],color=MUTED,alpha=.35,lw=.7)
            b.scatter(i-.10,r["SMC 1"],s=18,color=BLUE); b.scatter(i+.10,r["SMC 2"],s=18,color=ORANGE)
    b.axhline(0,color=GREY,lw=.7,ls="--"); b.set_xticks(range(3),modules,rotation=15,ha="right"); b.set_ylabel("Patient pseudobulk z-score"); b.set_title("Patient-resolved state programs"); panel_label(b,"b")
    q=effects[(effects.layer.eq("Within-state")) & effects.feature.isin(modules+["PDE8B"])]
    y=np.arange(4)[::-1]
    for ct,col,off in [("SMC 1",BLUE,-.08),("SMC 2",ORANGE,.08)]:
        z=q[q.cell_type.eq(ct)].set_index("feature").reindex(modules+["PDE8B"])
        for yy,val in zip(y+off,z.hedges_g): c.plot([0,val],[yy,yy],color=col,lw=1); c.scatter(val,yy,s=22,color=col)
    c.axvline(0,color=GREY,lw=.8,ls="--"); c.set_yticks(y,modules+["PDE8B"]); c.set_xlabel("Hedges g: PAH vs donor"); c.set_title("Disease effects within SMC states"); panel_label(c,"c")
    n=max((enr.side=="SMC1").sum(),(enr.side=="SMC2").sum()); yy=np.arange(n)[::-1]
    for side,sign,col in [("SMC1",-1,BLUE),("SMC2",1,ORANGE)]:
        z=enr[enr.side.eq(side)].reset_index(drop=True); xx=sign*(-np.log10(z.fdr.clip(lower=1e-50)).clip(upper=20))
        d.scatter(xx,yy[:len(z)],s=22+180*z.gene_ratio,color=col,edgecolor="white",linewidth=.5)
        for x0,y0,t in zip(xx,yy,z.term): d.text(x0+(.25 if sign>0 else -.25),y0,t,ha="left" if sign>0 else "right",va="center",fontsize=5.6)
    d.axvline(0,color=GREY,lw=.8); d.set_yticks([]); lim=21; d.set_xlim(-lim,lim); ticks=np.arange(-20,21,10); d.set_xticks(ticks,[str(abs(x)) for x in ticks]); d.set_xlabel("−log10(FDR)")
    d.set_title("Paired-patient pathway enrichment"); panel_label(d,"d")
    fig.text(.08,.025,"State differences do not establish lineage conversion from SMC1 to SMC2.",fontsize=6.3,color=MUTED)
    finish(fig,"Figure_4_SMC_state_continuum_and_enrichment")


def fig5_program():
    co=pd.read_csv(O/"pde8b-coexpression/pde8b_cross_cohort_coexpression_meta.csv")
    path=pd.read_csv(O/"pde8b-coexpression/enrichment/pde8b_core_pathway_enrichment.csv")
    classic=pd.read_csv(O/"pde8b-coexpression/state-deconvolution/pde8b_classic_pah_gene_connections.csv")
    lr=pd.read_csv(O/"gse293580-reanalysis/ligand-receptor-proxy/gse293580_mural_lr_patient_tests.csv")
    fig=plt.figure(figsize=(7.2,6.0)); gs=fig.add_gridspec(2,2,left=.10,right=.98,top=.91,bottom=.10,hspace=.42,wspace=.40)
    a,b,c,d=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    figure_title(fig,"PDE8B covaries with a vascular remodeling and cyclic-nucleotide program")
    co["mlog"]=-np.log10(co.adjusted_fdr.clip(lower=1e-300)); core=co.core_module.astype(bool)
    a.scatter(co.loc[~core,"adjusted_meta_r"],co.loc[~core,"mlog"],s=3,color="#D5DBDE",alpha=.45,rasterized=True)
    a.scatter(co.loc[core,"adjusted_meta_r"],co.loc[core,"mlog"],s=7,color=TEAL,alpha=.75,rasterized=True)
    a.axvline(0,color=GREY,lw=.7); a.set(xlabel="Meta correlation with PDE8B",ylabel="−log10(FDR)"); a.set_title("Cross-cohort coexpression landscape"); panel_label(a,"a")
    for _,r in co[core].nlargest(3,"adjusted_meta_r").iterrows():
        a.annotate(r.symbol,(r.adjusted_meta_r,r.mlog),xytext=(3,2),textcoords="offset points",fontsize=5.5)
    z=path.nsmallest(9,"fdr").sort_values("fdr",ascending=False); cols=[BLUE if x=="core_negative" else ORANGE for x in z.query_set]
    short=[]
    for term in z.term_name:
        term=term.replace("plasma membrane bounded cell projection","cell projection")
        term=term.replace("apoptotic process involved in luteolysis","apoptotic process")
        term=term.replace("cellular component organization or biogenesis","cellular organization")
        short.append(term[:30])
    b.barh(np.arange(len(z)),-np.log10(z.fdr),color=cols); b.set_yticks(np.arange(len(z)),short); b.set_xlabel("−log10(FDR)"); b.set_title("Core-module pathway enrichment"); panel_label(b,"b")
    top=classic.reindex(classic.adjusted_meta_r.abs().sort_values(ascending=False).index).head(14).sort_values("adjusted_meta_r")
    effect_lollipop(c,top.symbol,top.adjusted_meta_r,[ORANGE if x>0 else BLUE for x in top.adjusted_meta_r],xlabel="Meta correlation with PDE8B")
    c.set_title("Connections to established PAH genes"); panel_label(c,"c")
    q=lr[(lr.smc_boundary.eq("core_SMC")) & (lr.contrast.eq("SSc-PAH_vs_Donor"))]
    q=q[q.axis.isin(["ECM_integrin","TGFb","PDGF","JAG_NOTCH3"])]
    piv=q.pivot(index="axis",columns="direction",values="hedges_g"); x=np.arange(len(piv)); w=.34
    for i,(direction,col) in enumerate([("SMC_to_pericyte",ORANGE),("pericyte_to_SMC",TEAL)]): d.bar(x+(i-.5)*w,piv[direction],w,color=col,label=direction.replace("_"," "))
    d.axhline(0,color=GREY,lw=.8); d.set_xticks(x,[z.replace("_","–") for z in piv.index],rotation=25,ha="right"); d.set_ylabel("Patient-level Hedges g"); d.legend(fontsize=5.8); d.set_title("Exploratory mural-cell communication"); panel_label(d,"d")
    fig.text(.10,.025,"Coexpression and ligand–receptor proxy results are associative; interaction FDRs are not significant.",fontsize=6.3,color=MUTED)
    finish(fig,"Figure_5_PDE8B_vascular_program_and_interactions")


def main():
    fig1_bulk(); fig2_candidate(); fig3_single_cell(); fig4_smc(); fig5_program()
    print("\n".join(str(x) for x in sorted(OUT.glob("*.png"))))


if __name__ == "__main__":
    main()
