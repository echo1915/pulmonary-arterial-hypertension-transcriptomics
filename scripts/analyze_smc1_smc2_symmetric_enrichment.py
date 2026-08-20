from pathlib import Path
import argparse, gzip, json, ssl
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import ttest_rel
from urllib.request import Request, urlopen
import certifi

def main(raw: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    de_file=out/"gse210248_smc2_vs_smc1_patient_pseudobulk_de.csv"
    if de_file.exists():
        de=pd.read_csv(de_file).sort_values("paired_t")
        left=de.head(300).gene.tolist(); right=de.tail(300).gene.tolist()
        return run_enrichment(left,right,out)
    meta = pd.read_csv(raw / "GSE210248_Human_PA_Metadata.csv.gz", sep=";", index_col=0)
    base = raw / "gse210248/Human_PA_Integrated"
    with gzip.open(base / "barcodes.tsv.gz", "rt") as f: barcodes=[x.strip() for x in f]
    feat = pd.read_csv(base / "features.tsv.gz", sep="\t", header=None, compression="gzip")
    genes = feat.iloc[:,1].astype(str).to_numpy()
    X = mmread(base / "matrix.mtx.gz").tocsr()  # genes x cells
    if list(meta.index) != barcodes:
        meta = meta.loc[barcodes]
    keep = meta["Cell_annotation"].isin(["SMC 1","SMC 2"]).to_numpy()
    sm = meta.loc[keep, ["new.ident","Cell_annotation"]].copy()
    X = X[:,keep]
    groups=[]
    for sample in sorted(sm["new.ident"].unique()):
        for state in ["SMC 1","SMC 2"]:
            idx=((sm["new.ident"]==sample)&(sm["Cell_annotation"]==state)).to_numpy()
            counts=np.asarray(X[:,idx].sum(axis=1)).ravel()
            log=np.log1p(counts/counts.sum()*1e6)
            groups.append(pd.DataFrame({"gene":genes,"sample":sample,"state":state,"log_cpm":log}))
    pb=pd.concat(groups,ignore_index=True)
    w=pb.pivot_table(index="gene",columns=["sample","state"],values="log_cpm")
    samples=sorted(sm["new.ident"].unique())
    d=np.column_stack([w[(s,"SMC 2")]-w[(s,"SMC 1")] for s in samples])
    stat,p=ttest_rel(np.column_stack([w[(s,"SMC 2")] for s in samples]),np.column_stack([w[(s,"SMC 1")] for s in samples]),axis=1)
    de=pd.DataFrame({"gene":w.index,"mean_delta_smc2_minus_smc1":d.mean(1),"paired_t":stat,"p_value":p})
    de=de.replace([np.inf,-np.inf],np.nan).dropna().sort_values("paired_t")
    de.to_csv(de_file,index=False)
    # Balanced, expression-aware top lists; enrichment is descriptive due n=6.
    left=de.head(300).gene.tolist()       # SMC1 higher
    right=de.tail(300).gene.tolist()      # SMC2 higher
    run_enrichment(left,right,out)

def run_enrichment(left,right,out):
    def gp(query,label):
        payload={"organism":"hsapiens","query":query,"sources":["GO:BP","REAC"],"user_threshold":0.05,
                 "significance_threshold_method":"fdr","no_evidences":False}
        req=Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=json.dumps(payload).encode(),headers={"Content-Type":"application/json"})
        with urlopen(req,timeout=90,context=ssl.create_default_context(cafile=certifi.where())) as res: js=json.loads(res.read().decode())
        rows=[]
        for x in js.get("result",[]):
            rows.append({"side":label,"source":x["source"],"term_id":x["native"],"term":x["name"],
                         "fdr":x["p_value"],"intersection_size":x["intersection_size"],
                         "term_size":x["term_size"],"query_size":x["query_size"],
                         "gene_ratio":x["intersection_size"]/max(x["query_size"],1)})
        return rows
    enr=pd.DataFrame(gp(left,"SMC1")+gp(right,"SMC2"))
    enr.to_csv(out/"gse210248_smc1_smc2_gprofiler_enrichment.csv",index=False)
    (out/"gse210248_smc1_top300.txt").write_text("\n".join(left),encoding="utf-8")
    (out/"gse210248_smc2_top300.txt").write_text("\n".join(right),encoding="utf-8")

if __name__=="__main__":
    a=argparse.ArgumentParser(); a.add_argument("--raw",type=Path,required=True); a.add_argument("--out",type=Path,required=True)
    z=a.parse_args(); main(z.raw,z.out)
