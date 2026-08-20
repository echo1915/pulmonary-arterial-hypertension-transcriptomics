import gzip, math, itertools
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT="project-001-pulmonary-arterial-hypertension-transcriptomics"
ROOT=Path(r"G:\workdata\projects")/PROJECT
DATA=ROOT/"raw"/"current-data-snapshot"/"pah_scrna"/"gse293580"
COEX=ROOT/"outputs"/"pde8b-coexpression"; OUT=COEX/"single-cell-module"
GROUPS={"GSM8885793":"Donor","GSM8885794":"Donor","GSM8885795":"Donor","GSM8885796":"Donor","GSM8885797":"SSc-PAH","GSM8885798":"SSc-PAH","GSM8885799":"SSc-PAH","GSM8885800":"SSc-PAH","GSM8885801":"IPAH","GSM8885802":"IPAH","GSM8885803":"IPAH","GSM8885804":"IPAH"}
CONTRACTILE={"PRKG1","PDE5A","TRPC3","PLCB4","PRKAA2","PRKAB2","PPP1R12B","CALD1","RYR3","ROCK2","EDNRA","ATP1A2","PDE3B","PDE4C"}
def effect(a,b):
    a,b=np.asarray(a),np.asarray(b); sp=math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return 0 if sp==0 else ((a.mean()-b.mean())/sp)*(1-3/(4*(len(a)+len(b))-9))
def perm_p(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float); values=np.r_[a,b]; observed=abs(a.mean()-b.mean()); n=len(a)
    stats=[]
    for idx in itertools.combinations(range(len(values)),n):
        take=np.zeros(len(values),dtype=bool); take[list(idx)]=True
        stats.append(abs(values[take].mean()-values[~take].mean()))
    return sum(v>=observed-1e-12 for v in stats)/len(stats)
def main():
    OUT.mkdir(parents=True,exist_ok=True); mod=pd.read_csv(COEX/"pde8b_high_confidence_core_module.csv"); signs=dict(zip(mod.symbol,np.sign(mod.adjusted_meta_r)))
    rec=[]
    for matrix in sorted(DATA.glob("GSM*_matrix.mtx.gz")):
        acc=matrix.name.split('_')[0]; feat=Path(str(matrix).replace('matrix.mtx.gz','features.tsv.gz'))
        with gzip.open(feat,'rt') as h: genes=[x.rstrip().split('\t')[1] for x in h]
        wanted={i+1:g for i,g in enumerate(genes) if g in signs or g=='PDE8B'}; counts=defaultdict(float); lib=0
        with gzip.open(matrix,'rt') as h:
            for line in h:
                if line.startswith('%'): continue
                break
            for line in h:
                row,col,val=line.split(); row,val=int(row),float(val); lib+=val
                if row in wanted: counts[wanted[row]]+=val
        for gene in wanted.values(): rec.append({'accession':acc,'condition':GROUPS[acc],'gene':gene,'log_cp10k':math.log1p(10000*counts[gene]/lib)})
    x=pd.DataFrame(rec); x['z']=x.groupby('gene').log_cp10k.transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else 0)
    defs={'positive_core':lambda q:q.gene.map(signs).eq(1),'negative_core':lambda q:q.gene.map(signs).eq(-1),'signed_full':lambda q:q.gene.isin(signs),'contractile_cyclic':lambda q:q.gene.isin(CONTRACTILE)}
    scores=[]
    for name,fn in defs.items():
        b=x[fn(x)].copy(); b['part']=b.z*(b.gene.map(signs) if name=='signed_full' else 1)
        s=b.groupby(['accession','condition'],as_index=False).agg(score=('part','mean'),genes=('gene','nunique')); s['module']=name; scores.append(s)
    scores=pd.concat(scores); scores.to_csv(OUT/'gse293580_pde8b_module_scores.csv',index=False)
    rows=[]
    for name,b in scores.groupby('module'):
        donor=b.loc[b.condition.eq('Donor'),'score']
        for cond in ['IPAH','SSc-PAH']:
            a=b.loc[b.condition.eq(cond),'score']; rows.append({'module':name,'contrast':cond+'_vs_Donor','delta':a.mean()-donor.mean(),'hedges_g':effect(a,donor),'permutation_p':perm_p(a,donor),'genes':int(b.genes.max())})
    p=x[x.gene.eq('PDE8B')]
    donor=p.loc[p.condition.eq('Donor'),'log_cp10k']
    for cond in ['IPAH','SSc-PAH']:
        a=p.loc[p.condition.eq(cond),'log_cp10k']; rows.append({'module':'PDE8B_gene','contrast':cond+'_vs_Donor','delta':a.mean()-donor.mean(),'hedges_g':effect(a,donor),'permutation_p':perm_p(a,donor),'genes':1})
    out=pd.DataFrame(rows); out.to_csv(OUT/'gse293580_pde8b_module_effects.csv',index=False); print(out.to_string(index=False))
if __name__=='__main__': main()
