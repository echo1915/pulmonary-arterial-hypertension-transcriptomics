import gzip, itertools, math, json
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT="project-001-pulmonary-arterial-hypertension-transcriptomics"
ROOT=Path(r"G:\workdata\projects")/PROJECT
RAW=ROOT/"raw"/"current-data-snapshot"/"pah_scrna"
COEX=ROOT/"outputs"/"pde8b-coexpression"; OUT=COEX/"state-deconvolution"; OUT.mkdir(parents=True,exist_ok=True)

def g(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float); sp=math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return 0 if sp==0 else ((a.mean()-b.mean())/sp)*(1-3/(4*(len(a)+len(b))-9))
def perm(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float); v=np.r_[a,b]; obs=abs(a.mean()-b.mean()); n=len(a); z=[]
    for idx in itertools.combinations(range(len(v)),n):
        q=np.zeros(len(v),bool); q[list(idx)]=1; z.append(abs(v[q].mean()-v[~q].mean()))
    return sum(x>=obs-1e-12 for x in z)/len(z)

def main():
    meta=pd.read_csv(RAW/"GSE210248_Human_PA_Metadata.csv.gz",sep=";",decimal=",",index_col=0)
    meta["sample"]=meta["new.ident"].astype(str); meta["disease"]=np.where(meta["sample"].str.startswith("PAH"),"PAH","Donor")
    tab=meta.groupby(["sample","disease","Cell_annotation"]).size().rename("cells").reset_index()
    totals=meta.groupby("sample").size(); tab["proportion"]=tab.cells/tab['sample'].map(totals)
    tab.to_csv(OUT/"gse210248_celltype_proportions.csv",index=False)
    effects=[]
    for ct,b in tab.groupby("Cell_annotation"):
        a=b.loc[b.disease.eq('PAH'),'proportion']; d=b.loc[b.disease.eq('Donor'),'proportion']
        if len(a)==3 and len(d)==3: effects.append({'cell_type':ct,'pah_mean':a.mean(),'donor_mean':d.mean(),'delta':a.mean()-d.mean(),'hedges_g':g(a,d),'permutation_p':perm(a,d)})
    eff=pd.DataFrame(effects).sort_values('delta',ascending=False); eff.to_csv(OUT/"gse210248_celltype_proportion_effects.csv",index=False)

    scores=pd.read_csv(COEX/"single-cell-module"/"gse210248_pde8b_module_score_variants.csv")
    prop=tab.rename(columns={'Cell_annotation':'cell_type'})[['sample','cell_type','proportion']]
    assoc=scores.merge(prop,on=['sample','cell_type'])
    rows=[]
    for (module,ct),b in assoc.groupby(['module','cell_type']):
        if len(b)>=5:
            rows.append({'module':module,'cell_type':ct,'n_samples':len(b),'spearman_module_vs_proportion':b.score.rank().corr(b.proportion.rank())})
    pd.DataFrame(rows).to_csv(OUT/"gse210248_module_proportion_correlations.csv",index=False)

    co=pd.read_csv(COEX/"pde8b_cross_cohort_coexpression_meta.csv")
    classic=['BMPR2','KCNK3','SOX17','TBX4','CAV1','ENG','ACVRL1','SMAD9','EIF2AK4','ATP13A3','AQP1','GDF2','KDR','ABCC8','KLK1','EDNRA','ROCK2','NOS3','PDE5A','PRKG1','TRPC3','PLCB4']
    c=co[co.symbol.isin(classic)].copy().sort_values('adjusted_meta_r',ascending=False)
    c.to_csv(OUT/"pde8b_classic_pah_gene_connections.csv",index=False)
    summary={'SMC_proportions':eff[eff.cell_type.isin(['SMC 1','SMC 2'])].to_dict('records'),'classic_connections':c[['symbol','adjusted_meta_r','adjusted_fdr','adjusted_I2','pah_only_meta_r','pah_only_fdr']].to_dict('records')}
    (OUT/"state_deconvolution_summary.json").write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('Proportion effects:'); print(eff.to_string(index=False)); print('\nClassic PAH connections:'); print(c[['symbol','adjusted_meta_r','adjusted_fdr','adjusted_I2','pah_only_meta_r','pah_only_fdr']].to_string(index=False))
if __name__=='__main__': main()
