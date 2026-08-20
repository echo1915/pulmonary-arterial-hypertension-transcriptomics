import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = "project-001-pulmonary-arterial-hypertension-transcriptomics"
DATA_ROOT = Path(r"G:\workdata\projects") / PROJECT
RAW = DATA_ROOT / "raw" / "current-data-snapshot" / "pah_scrna"
COEX = DATA_ROOT / "outputs" / "pde8b-coexpression"
OUT = COEX / "single-cell-module"

def hedges_g(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2: return np.nan
    sp = math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    if sp == 0: return 0.0
    return ((a.mean()-b.mean())/sp)*(1-3/(4*(len(a)+len(b))-9))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = RAW / "gse210248" / "Human_PA_Integrated"
    meta = pd.read_csv(RAW / "GSE210248_Human_PA_Metadata.csv.gz", sep=";", decimal=",", index_col=0)
    module = pd.read_csv(COEX / "pde8b_high_confidence_core_module.csv")
    signs = dict(zip(module.symbol.astype(str), np.sign(module.adjusted_meta_r)))
    with gzip.open(base / "features.tsv.gz", "rt") as h: features=[x.rstrip().split("\t")[1] for x in h]
    with gzip.open(base / "barcodes.tsv.gz", "rt") as h: barcodes=[x.strip() for x in h]
    if barcodes != list(meta.index): raise ValueError("Barcode/metadata mismatch")
    wanted={i+1:g for i,g in enumerate(features) if g in signs}
    groups=list(zip(meta["new.ident"].astype(str),meta["Cell_annotation"].astype(str)))
    ncells=defaultdict(int)
    for g in groups: ncells[g]+=1
    lib=defaultdict(float); counts=defaultdict(float); detected=defaultdict(int)
    with gzip.open(base/"matrix.mtx.gz","rt") as h:
        for line in h:
            if line.startswith("%"): continue
            if len(line.split())==3: break
        for line in h:
            row,col,val=line.split(); row,col,val=int(row),int(col),float(val); grp=groups[col-1]
            lib[grp]+=val
            if row in wanted:
                gene=wanted[row]; counts[(gene,)+grp]+=val; detected[(gene,)+grp]+=1
    records=[]
    for gene in sorted(wanted.values()):
        for (sample,celltype),n in sorted(ncells.items()):
            c=counts[(gene,sample,celltype)]; total=lib[(sample,celltype)]
            records.append({"gene":gene,"sample":sample,"disease":"PAH" if sample.startswith("PAH") else "Donor",
                            "cell_type":celltype,"n_cells":n,"log_cp10k":math.log1p(10000*c/total) if total else 0,
                            "detection_fraction":detected[(gene,sample,celltype)]/n,"module_sign":int(signs[gene])})
    pb=pd.DataFrame(records); pb.to_csv(OUT/"gse210248_pde8b_core_gene_pseudobulk.csv",index=False)
    pb["z"] = pb.groupby(["cell_type","gene"])["log_cp10k"].transform(lambda x:(x-x.mean())/x.std(ddof=1) if x.std(ddof=1)>0 else 0)
    pb["signed_z"] = pb.z*pb.module_sign
    scores=(pb.groupby(["sample","disease","cell_type"],as_index=False)
              .agg(module_score=("signed_z","mean"),genes_detected=("detection_fraction",lambda x:int((x>0).sum())),genes_total=("gene","size"),n_cells=("n_cells","first")))
    scores.to_csv(OUT/"gse210248_pde8b_module_scores.csv",index=False)
    effects=[]
    for ct,b in scores.groupby("cell_type"):
        a=b.loc[b.disease.eq("PAH"),"module_score"]; c=b.loc[b.disease.eq("Donor"),"module_score"]
        if len(a)==3 and len(c)==3:
            effects.append({"cell_type":ct,"pah_mean":a.mean(),"donor_mean":c.mean(),"delta":a.mean()-c.mean(),"hedges_g":hedges_g(a,c),
                            "mean_genes_detected":b.genes_detected.mean(),"mean_cells":b.n_cells.mean()})
    effects=pd.DataFrame(effects).sort_values("delta",ascending=False)
    effects.to_csv(OUT/"gse210248_pde8b_module_celltype_effects.csv",index=False)
    summary={"requested_core_genes":len(signs),"present_genes":len(wanted),"cells":len(meta),"samples":meta['new.ident'].nunique(),
             "cell_types_tested":len(effects),"top_increased":effects.head(10).to_dict('records'),"top_decreased":effects.tail(10).to_dict('records'),
             "caution":"Sample-level descriptive effect sizes (3 PAH vs 3 donor); no cell-level pseudo-replication."}
    (OUT/"gse210248_pde8b_module_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if not k.startswith('top_')},indent=2)); print(effects.head(15).to_string(index=False))

if __name__ == "__main__": main()
