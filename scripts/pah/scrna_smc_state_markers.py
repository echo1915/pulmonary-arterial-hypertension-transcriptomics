import gzip, math, itertools
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
PROJECT="project-001-pulmonary-arterial-hypertension-transcriptomics"; ROOT=Path(r"G:\workdata\projects")/PROJECT
RAW=ROOT/"raw"/"current-data-snapshot"/"pah_scrna"; OUT=ROOT/"outputs"/"pde8b-coexpression"/"state-deconvolution"
PANELS={'PDE8B':['PDE8B'],'contractile':['ACTA2','TAGLN','MYH11','CNN1','MYL9','TPM2','CALD1'],'synthetic_ecm':['KLF4','COL1A1','COL3A1','FN1','VIM'],'mural_remodeling':['RGS5','CSPG4','MCAM','NOTCH3','PDGFRB']}
def g(a,b):
 a,b=np.asarray(a),np.asarray(b); sp=math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)); return 0 if sp==0 else ((a.mean()-b.mean())/sp)*(1-3/(4*(len(a)+len(b))-9))
def main():
 base=RAW/'gse210248'/'Human_PA_Integrated'; meta=pd.read_csv(RAW/'GSE210248_Human_PA_Metadata.csv.gz',sep=';',decimal=',',index_col=0)
 with gzip.open(base/'features.tsv.gz','rt') as h: features=[x.rstrip().split('\t')[1] for x in h]
 with gzip.open(base/'barcodes.tsv.gz','rt') as h: barcodes=[x.strip() for x in h]
 if barcodes!=list(meta.index): raise ValueError('alignment')
 wanted_genes={g for v in PANELS.values() for g in v}; wanted={i+1:g for i,g in enumerate(features) if g in wanted_genes}; groups=list(zip(meta['new.ident'].astype(str),meta['Cell_annotation'].astype(str)))
 lib=defaultdict(float); counts=defaultdict(float)
 with gzip.open(base/'matrix.mtx.gz','rt') as h:
  for line in h:
   if line.startswith('%'): continue
   break
  for line in h:
   row,col,val=line.split(); row,col,val=int(row),int(col),float(val); grp=groups[col-1]; lib[grp]+=val
   if row in wanted: counts[(wanted[row],)+grp]+=val
 rec=[]
 for gene in wanted.values():
  for sample,ct in sorted(set(groups)):
   rec.append({'gene':gene,'sample':sample,'disease':'PAH' if sample.startswith('PAH') else 'Donor','cell_type':ct,'log_cp10k':math.log1p(10000*counts[(gene,sample,ct)]/lib[(sample,ct)])})
 x=pd.DataFrame(rec); x=x[x.cell_type.isin(['SMC 1','SMC 2'])]; x.to_csv(OUT/'gse210248_smc_marker_pseudobulk.csv',index=False)
 rows=[]
 for (gene,ct),b in x.groupby(['gene','cell_type']):
  a=b.loc[b.disease.eq('PAH'),'log_cp10k']; d=b.loc[b.disease.eq('Donor'),'log_cp10k']; rows.append({'feature':gene,'panel':'gene','cell_type':ct,'delta':a.mean()-d.mean(),'hedges_g':g(a,d)})
 x['z']=x.groupby(['cell_type','gene']).log_cp10k.transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else 0)
 for panel,genes in PANELS.items():
  s=x[x.gene.isin(genes)].groupby(['sample','disease','cell_type'],as_index=False).agg(score=('z','mean'),genes=('gene','nunique'))
  for ct,b in s.groupby('cell_type'):
   a=b.loc[b.disease.eq('PAH'),'score']; d=b.loc[b.disease.eq('Donor'),'score']; rows.append({'feature':panel,'panel':'module','cell_type':ct,'delta':a.mean()-d.mean(),'hedges_g':g(a,d),'genes':int(b.genes.max())})
 out=pd.DataFrame(rows); out.to_csv(OUT/'gse210248_smc_state_marker_effects.csv',index=False); print(out[out.panel.eq('module')].to_string(index=False)); print('\nGenes:'); print(out[out.panel.eq('gene')].sort_values(['cell_type','hedges_g']).to_string(index=False))
if __name__=='__main__': main()
