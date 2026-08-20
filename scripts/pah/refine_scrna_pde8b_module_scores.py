import json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(r"G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\outputs\pde8b-coexpression")
IN=ROOT/"single-cell-module"/"gse210248_pde8b_core_gene_pseudobulk.csv"
OUT=ROOT/"single-cell-module"
CONTRACTILE={"PRKG1","PDE5A","TRPC3","PLCB4","PRKAA2","PRKAB2","PPP1R12B","CALD1","RYR3","ROCK2","EDNRA","ATP1A2","PDE3B","PDE4C"}

def g(a,b):
    a,b=np.asarray(a),np.asarray(b); sp=math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return 0 if sp==0 else ((a.mean()-b.mean())/sp)*(1-3/(4*(len(a)+len(b))-9))

def main():
    x=pd.read_csv(IN)
    x["z_within_celltype"]=x.groupby(["cell_type","gene"])["log_cp10k"].transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else 0)
    x["z_across_celltypes"]=x.groupby("gene")["log_cp10k"].transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else 0)
    definitions={"positive_core":x.module_sign>0,"negative_core":x.module_sign<0,"signed_full":pd.Series(True,index=x.index),"contractile_cyclic":x.gene.isin(CONTRACTILE)}
    all_scores=[]; all_effects=[]; localization=[]
    for name,mask in definitions.items():
        b=x[mask].copy(); sign=b.module_sign if name=="signed_full" else 1
        b["scorepart"]=b.z_within_celltype*sign
        s=b.groupby(["sample","disease","cell_type"],as_index=False).agg(score=("scorepart","mean"),genes=("gene","nunique"))
        s["module"]=name; all_scores.append(s)
        for ct,q in s.groupby("cell_type"):
            a=q.loc[q.disease.eq("PAH"),"score"]; c=q.loc[q.disease.eq("Donor"),"score"]
            all_effects.append({"module":name,"cell_type":ct,"delta":a.mean()-c.mean(),"hedges_g":g(a,c),"genes":int(q.genes.max())})
        loc=b.assign(locpart=b.z_across_celltypes*sign).groupby("cell_type",as_index=False).agg(localization_score=("locpart","mean"),genes=("gene","nunique"))
        loc["module"]=name; localization.append(loc)
    scores=pd.concat(all_scores); effects=pd.DataFrame(all_effects); loc=pd.concat(localization)
    scores.to_csv(OUT/"gse210248_pde8b_module_score_variants.csv",index=False)
    effects.to_csv(OUT/"gse210248_pde8b_module_variant_effects.csv",index=False)
    loc.to_csv(OUT/"gse210248_pde8b_module_localization.csv",index=False)
    focus=effects[effects.cell_type.isin(["SMC 1","SMC 2","Endo 1","Endo 2","Fibro"])].sort_values(["module","cell_type"])
    print(focus.to_string(index=False)); print("\nLocalization top:"); print(loc.sort_values(["module","localization_score"],ascending=[True,False]).groupby("module").head(4).to_string(index=False))

if __name__=="__main__": main()
