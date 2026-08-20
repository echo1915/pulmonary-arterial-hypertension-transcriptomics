import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from figure_style import COLORS, apply_publication_style, panel_label, save_figure

PROJECT='project-001-pulmonary-arterial-hypertension-transcriptomics'
ROOT=Path(r'G:\workdata\projects')/PROJECT
AUDIT=ROOT/'outputs'/'current-results'/'pah_audit'; COEX=ROOT/'outputs'/'pde8b-coexpression'
OUT=ROOT/'outputs'/'supplementary-figures-sci-v1'; OUT.mkdir(parents=True,exist_ok=True)
apply_publication_style()
BLUE=COLORS['smc1']; TEAL=COLORS['teal']; ORANGE=COLORS['smc2']; RED=COLORS['signal']; NAVY=COLORS['ink']; GREY=COLORS['muted']; LIGHT=COLORS['pale']
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
plt.rcParams['font.size']=8
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.spines.right']=False
plt.rcParams['pdf.fonttype']=42
plt.rcParams['svg.fonttype']='none'

def panel(ax,label): panel_label(ax,label)
def save(fig,name):
    save_figure(fig,OUT,name)

def figure1():
    fig,ax=plt.subplots(figsize=(11.5,6.2)); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.text(.02,.95,'Figure 1. Study design and multi-stage prioritization of PDE8B',fontsize=15,fontweight='bold',color=NAVY)
    cohorts=[('GSE15197','18 PAH / 13 donor'),('GSE113439','14 PAH / 11 donor'),('GSE254617','82 PAH / 52 donor'),('GSE208592','15 PAH / 18 donor')]
    for i,(a,b) in enumerate(cohorts):
        x=.04+i*.235; p=FancyBboxPatch((x,.72),.19,.12,boxstyle='round,pad=.012',fc='#EDF4F8',ec=BLUE,lw=1.2); ax.add_patch(p); ax.text(x+.095,.79,a,ha='center',fontweight='bold'); ax.text(x+.095,.75,b,ha='center',color=GREY)
    stages=[('15,165','genes shared'),('251','meta FDR < 0.05'),('188','robust genes'),('specificity','PAH vs IPF-PH'),('single cell','localization'),('PDE8B','lead signal')]
    for i,(a,b) in enumerate(stages):
        x=.035+i*.16; col=RED if i==5 else TEAL
        ax.add_patch(FancyBboxPatch((x,.39),.13,.13,boxstyle='round,pad=.012',fc='white',ec=col,lw=1.4)); ax.text(x+.065,.465,a,ha='center',fontweight='bold',color=col,fontsize=10); ax.text(x+.065,.425,b,ha='center',color=GREY,fontsize=7)
        if i<5: ax.add_patch(FancyArrowPatch((x+.132,.455),(x+.157,.455),arrowstyle='-|>',mutation_scale=10,color=NAVY))
    ax.text(.05,.25,'Independent support',fontweight='bold',color=NAVY,fontsize=10)
    notes=['GSE53408 lung validation: g = 1.08','GSE210248: SMC1/SMC2 localization','GSE293580: IPAH g = 1.33;\\nSSc-PAH g = 0.72']
    for i,t in enumerate(notes): ax.text(.07,.20-i*.060,u'\u2022 '+t,fontsize=8.5,va='top',linespacing=1.15)
    ax.text(.68,.22,'Core interpretation',fontweight='bold',color=NAVY,fontsize=10)
    ax.text(.70,.16,'PDE8B marks an SMC-enriched\ncyclic-nucleotide/remodeling state;\ncausality remains untested.',fontsize=10,color=RED)
    save(fig,'Figure_1_study_design_and_selection')

def figure2():
    meta=pd.read_csv(AUDIT/'lung_mechanism_random_effects_meta.csv'); p=meta[meta.symbol.eq('PDE8B')].iloc[0]
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1,1.35]}); panel(axs[0],'A'); panel(axs[1],'B')
    names=['GSE15197','GSE113439','GSE254617','GSE208592','Random-effects meta']; vals=[p.g_GSE15197,p.g_GSE113439,p.g_GSE254617,p.g_GSE208592,p.meta_g]
    y=np.arange(5)[::-1]; axs[0].axvline(0,color=GREY,lw=.8); axs[0].scatter(vals[:4],y[:4],s=45,c=BLUE,zorder=3); axs[0].errorbar(p.meta_g,y[4],xerr=[[p.meta_g-p.ci_low],[p.ci_high-p.meta_g]],fmt='D',c=RED,capsize=3)
    axs[0].set_yticks(y,names); axs[0].set_xlabel('Hedges g (PAH vs donor)'); axs[0].set_title('PDE8B across lung cohorts',fontweight='bold'); axs[0].text(.03,-.22,f'Meta g={p.meta_g:.2f} (95% CI {p.ci_low:.2f} to {p.ci_high:.2f}); I²={p.I2:.0f}%',transform=axs[0].transAxes,color=RED)
    m=meta.copy(); m['mlog']=-np.log10(m.fdr.clip(lower=1e-300)); robust=m.same_direction & (m.fdr<.05) & (m.I2<50)
    axs[1].scatter(m.loc[~robust,'meta_g'],m.loc[~robust,'mlog'],s=3,c='#CCD5DB',alpha=.45,rasterized=True); axs[1].scatter(m.loc[robust,'meta_g'],m.loc[robust,'mlog'],s=6,c=TEAL,alpha=.7,rasterized=True)
    q=m[m.symbol.eq('PDE8B')].iloc[0]; axs[1].scatter(q.meta_g,q.mlog,s=55,c=RED,zorder=4); axs[1].annotate('PDE8B',(q.meta_g,q.mlog),xytext=(8,8),textcoords='offset points',fontweight='bold',color=RED)
    axs[1].axvline(0,color=GREY,lw=.7); axs[1].set(xlabel='Random-effects Hedges g',ylabel='−log10(FDR)',title='Cross-cohort lung meta-analysis'); axs[1].set_title('Cross-cohort lung meta-analysis',fontweight='bold')
    fig.suptitle('Figure 2. Reproducible upregulation of PDE8B in PAH lung',fontsize=14,fontweight='bold',color=NAVY); fig.tight_layout(); save(fig,'Figure_2_bulk_meta_analysis')

def figure3():
    u=pd.read_csv(ROOT/'processed'/'gse210248_umap'/'gse210248_recomputed_umap.csv.gz')
    harmony=pd.read_csv(ROOT/'processed'/'gse210248_harmony'/'gse210248_raw_vs_harmony_umap.csv.gz',usecols=['barcode','harmony_UMAP1','harmony_UMAP2'])
    u=u.merge(harmony,on='barcode',how='left',validate='one_to_one')
    if u[['harmony_UMAP1','harmony_UMAP2']].isna().any().any(): raise ValueError('Harmony coordinate merge failed')
    u['UMAP1']=u.harmony_UMAP1; u['UMAP2']=u.harmony_UMAP2
    prop=pd.read_csv(COEX/'state-deconvolution'/'gse210248_celltype_proportions.csv')
    cats=sorted(u.cell_type.unique()); cmap=plt.cm.tab20(np.linspace(0,1,len(cats))); colors=dict(zip(cats,cmap))
    fig,axs=plt.subplots(1,3,figsize=(13,4.4));
    for ct in cats:
        q=u[u.cell_type.eq(ct)]; axs[0].scatter(q.UMAP1,q.UMAP2,s=1.2,color=colors[ct],alpha=.65,rasterized=True,label=ct)
        axs[0].text(q.UMAP1.median(),q.UMAP2.median(),ct,fontsize=5.5)
    axs[0].set_title('Harmony-integrated cell-type landscape'); axs[0].set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2'); panel(axs[0],'A')
    zero=u.PDE8B_log_cp10k<=0; axs[1].scatter(u.loc[zero,'UMAP1'],u.loc[zero,'UMAP2'],s=1,c='#D9DEE2',alpha=.35,rasterized=True)
    pos=u[~zero].sort_values('PDE8B_log_cp10k'); sc=axs[1].scatter(pos.UMAP1,pos.UMAP2,c=pos.PDE8B_log_cp10k,s=4,cmap='magma',rasterized=True); fig.colorbar(sc,ax=axs[1],shrink=.65,label='log(CP10K+1)')
    axs[1].set_title('PDE8B on Harmony UMAP'); axs[1].set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2'); panel(axs[1],'B')
    focus=prop[prop.Cell_annotation.isin(['SMC 1','SMC 2','Endo 1','Endo 2','Fibro'])].copy(); order=['SMC 1','SMC 2','Endo 1','Endo 2','Fibro']; x=np.arange(len(order)); width=.32
    for i,(d,c) in enumerate([('Donor',BLUE),('PAH',RED)]):
        means=focus[focus.disease.eq(d)].groupby('Cell_annotation').proportion.mean().reindex(order); sem=focus[focus.disease.eq(d)].groupby('Cell_annotation').proportion.sem().reindex(order)
        axs[2].bar(x+(i-.5)*width,means,width,yerr=sem,color=c,alpha=.85,label=d,capsize=2)
    axs[2].set_xticks(x,order,rotation=35,ha='right'); axs[2].set_ylabel('Fraction of cells'); axs[2].set_title('Structural-cell composition'); axs[2].legend(frameon=False); panel(axs[2],'C')
    fig.text(.01,.01,'Harmony batch variable: sample; disease status retained. Cell annotations and all downstream pseudobulk statistics are unchanged.',fontsize=7,color=GREY)
    fig.suptitle('Figure 3. Single-cell localization of PDE8B in PAH pulmonary artery',fontsize=14,fontweight='bold',color=NAVY); fig.tight_layout(); save(fig,'Figure_3_scRNA_UMAP_and_localization')

def figure4():
    e=pd.read_csv(COEX/'state-deconvolution'/'gse210248_smc_state_marker_effects.csv'); modules=e[e.panel.eq('module')]; genes=e[e.panel.eq('gene')]
    fig=plt.figure(figsize=(12,7)); gs=fig.add_gridspec(2,3,width_ratios=[1,1.25,1.25],height_ratios=[1,1],hspace=.48,wspace=.42)
    ax=fig.add_subplot(gs[:,0]); panel(ax,'A'); names=['PDE8B','contractile','synthetic_ecm','mural_remodeling']; lab=['PDE8B gene','Contractile module','Synthetic/ECM module','Mural-remodeling module']; y=np.arange(4)
    for j,(ct,c) in enumerate([('SMC 1',BLUE),('SMC 2',RED)]):
        q=modules[modules.cell_type.eq(ct)].set_index('feature').reindex(names); ax.barh(y+(j-.5)*.32,q.hedges_g,.32,color=c,label=ct)
    ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(y,lab); ax.invert_yaxis(); ax.set_xlabel('Hedges g (PAH vs donor)')
    ax.set_ylabel('Feature or predefined module'); ax.set_title('Gene and state-program effects',fontweight='bold'); ax.legend(frameon=False)
    axh=fig.add_subplot(gs[:,1]); panel(axh,'B'); mark=['PDE8B','ACTA2','TAGLN','MYH11','CNN1','MYL9','TPM2','COL1A1','COL3A1','FN1','VIM','RGS5','CSPG4','NOTCH3']
    mat=genes.pivot(index='feature',columns='cell_type',values='hedges_g').reindex(mark)[['SMC 1','SMC 2']]; im=axh.imshow(mat,cmap='RdBu_r',vmin=-2.5,vmax=2.5,aspect='auto'); axh.set_xticks([0,1],['SMC 1','SMC 2']); axh.set_yticks(np.arange(len(mark)),mark); axh.set_title('Marker-level effect sizes',fontweight='bold'); fig.colorbar(im,ax=axh,shrink=.55,label='Hedges g')
    ax2=fig.add_subplot(gs[0,2]); panel(ax2,'C'); ax2.axis('off'); ax2.set_xlim(0,1); ax2.set_ylim(0,1)
    for x0,title,sub,col in [(.03,'Contractile-like SMC','ACTA2 TAGLN MYH11\nCNN1 MYL9 TPM2',BLUE),(.63,'Remodeling-associated SMC','PDE8B COL3A1 FN1\nVIM / ECM program',RED)]:
        ax2.add_patch(FancyBboxPatch((x0,.30),.34,.42,boxstyle='round,pad=.02',fc='white',ec=col,lw=2)); ax2.text(x0+.17,.60,title,ha='center',fontweight='bold',color=col); ax2.text(x0+.17,.43,sub,ha='center',va='center')
    ax2.add_patch(FancyArrowPatch((.39,.51),(.61,.51),arrowstyle='simple',mutation_scale=18,color=ORANGE,alpha=.85)); ax2.text(.50,.64,'PAH-associated state shift',ha='center',fontweight='bold',color=ORANGE)
    ax3=fig.add_subplot(gs[1,2]); panel(ax3,'D'); p=pd.read_csv(COEX/'state-deconvolution'/'gse210248_celltype_proportion_effects.csv'); q=p[p.cell_type.isin(['SMC 1','SMC 2'])]
    ax3.bar(q.cell_type,q.hedges_g,color=[BLUE,RED]); ax3.axhline(0,color=GREY,lw=.8); ax3.set_ylabel('Hedges g for cell fraction'); ax3.set_title('SMC abundance difference',fontweight='bold'); ax3.text(.02,-.28,'n=3 PAH and 3 donors; descriptive, not a lineage trajectory',transform=ax3.transAxes,color=GREY)
    fig.suptitle('Figure 4. PDE8B upregulation accompanies a PAH-associated SMC state shift',fontsize=14,fontweight='bold',color=NAVY); save(fig,'Figure_4_SMC_state_transition')

def figure5():
    m=pd.read_csv(COEX/'pde8b_cross_cohort_coexpression_meta.csv'); enr=pd.read_csv(COEX/'enrichment'/'pde8b_core_pathway_enrichment.csv'); rep=pd.read_csv(COEX/'single-cell-module'/'gse293580_pde8b_module_effects.csv')
    fig,axs=plt.subplots(1,3,figsize=(13,4.5)); panel(axs[0],'A'); panel(axs[1],'B'); panel(axs[2],'C')
    m['mlog']=-np.log10(m.adjusted_fdr.clip(lower=1e-300)); core=m.core_module.astype(bool); axs[0].scatter(m.loc[~core,'adjusted_meta_r'],m.loc[~core,'mlog'],s=2,c='#D5DDE2',alpha=.4,rasterized=True); axs[0].scatter(m.loc[core,'adjusted_meta_r'],m.loc[core,'mlog'],s=6,c=ORANGE,alpha=.7,rasterized=True)
    for g0 in ['TRPC3','PRKG1','PDE5A','PLCB4','BMPR2','KCNK3']:
        q=m[m.symbol.eq(g0)].iloc[0]; axs[0].annotate(g0,(q.adjusted_meta_r,q.mlog),fontsize=6,xytext=(3,3),textcoords='offset points')
    axs[0].set(xlabel='Meta-correlation with PDE8B',ylabel='−log10(FDR)',title='Stable cross-cohort co-expression')
    wanted=['cilium organization','actin filament-based process','muscle contraction','cGMP-PKG signaling pathway','smooth muscle contraction','cellular response to calcium ion']; q=enr[enr.term_name.isin(wanted)].drop_duplicates('term_name').sort_values('fdr'); yy=np.arange(len(q)); axs[1].scatter(-np.log10(q.fdr),yy,s=15+q.intersection_size*4,c=q.intersection_size,cmap='viridis'); axs[1].set_yticks(yy,q.term_name); axs[1].invert_yaxis(); axs[1].set_xlabel('−log10(FDR)'); axs[1].set_title('Functional enrichment')
    keep=rep[rep.module.isin(['PDE8B_gene','positive_core','contractile_cyclic','signed_full'])].copy(); order=['PDE8B_gene','positive_core','contractile_cyclic','signed_full']; xx=np.arange(4); width=.34
    for i,(ct,c) in enumerate([('IPAH_vs_Donor',RED),('SSc-PAH_vs_Donor',TEAL)]):
        z=keep[keep.contrast.eq(ct)].set_index('module').reindex(order); axs[2].bar(xx+(i-.5)*width,z.hedges_g,width,color=c,label=ct.replace('_vs_Donor',''))
    axs[2].axhline(0,color=GREY,lw=.8); axs[2].set_xticks(xx,['PDE8B\ngene','Positive-core\nmodule','Contractile/cyclic\nmodule','Signed-full\nmodule'],rotation=25); axs[2].set_ylabel('Hedges g (PAH vs donor)'); axs[2].set_title('Independent replication by PAH etiology')
    axs[2].text(.02,.98,'Positive values indicate higher expression/module score in PAH',transform=axs[2].transAxes,va='top',fontsize=6.5,color=GREY); axs[2].legend(frameon=False)
    fig.suptitle('Figure 5. PDE8B-centered pathways and etiologic heterogeneity',fontsize=14,fontweight='bold',color=NAVY); fig.tight_layout(); save(fig,'Figure_5_coexpression_pathways_replication')

def supplements():
    loco=pd.read_csv(COEX/'pde8b_module_leave_one_cohort_out.csv'); genes=['TRPC3','PRKG1','PDE5A','PLCB4','PRKAA2','PRKAB2','PPP1R12B','CALD1','BMPR2','KCNK3']; q=loco[loco.symbol.isin(genes)].pivot(index='symbol',columns='omitted_cohort',values='meta_r').reindex(genes)
    fig,ax=plt.subplots(figsize=(6.5,5)); im=ax.imshow(q,cmap='RdBu_r',vmin=-.65,vmax=.65); ax.set_xticks(range(q.shape[1]),q.columns,rotation=35,ha='right'); ax.set_yticks(range(q.shape[0]),q.index); plt.colorbar(im,ax=ax,label='Leave-one-cohort-out meta-r'); ax.set_title('Supplementary Figure 1. Sensitivity of key PDE8B connections',fontweight='bold'); fig.tight_layout(); save(fig,'Supplementary_Figure_1_LOCO')
    c=pd.read_csv(COEX/'state-deconvolution'/'pde8b_classic_pah_gene_connections.csv').sort_values('adjusted_meta_r'); fig,ax=plt.subplots(figsize=(6.5,6)); y=np.arange(len(c)); ax.hlines(y,0,c.adjusted_meta_r,color=LIGHT); ax.scatter(c.adjusted_meta_r,y,c=np.where(c.adjusted_fdr<.05,TEAL,GREY),s=30); ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(y,c.symbol); ax.set_xlabel('Meta-correlation with PDE8B'); ax.set_title('Supplementary Figure 2. Connections to established PAH genes',fontweight='bold'); fig.tight_layout(); save(fig,'Supplementary_Figure_2_classic_PAH_connections')

def main():
    figure1(); figure2(); figure3(); figure4(); figure5(); supplements(); print('\n'.join(str(x) for x in sorted(OUT.glob('*.png'))))
if __name__=='__main__': main()






