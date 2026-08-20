import gzip, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
import umap

PROJECT="project-001-pulmonary-arterial-hypertension-transcriptomics"
ROOT=Path(r"G:\workdata\projects")/PROJECT
RAW=ROOT/"raw"/"current-data-snapshot"/"pah_scrna"; BASE=RAW/"gse210248"/"Human_PA_Integrated"
OUT=ROOT/"processed"/"gse210248_umap"; OUT.mkdir(parents=True,exist_ok=True)

def main():
    meta=pd.read_csv(RAW/"GSE210248_Human_PA_Metadata.csv.gz",sep=";",decimal=",",index_col=0)
    with gzip.open(BASE/"barcodes.tsv.gz","rt") as h: barcodes=[x.strip() for x in h]
    with gzip.open(BASE/"features.tsv.gz","rt") as h: genes=[x.rstrip().split('\t')[1] for x in h]
    if barcodes!=list(meta.index): raise ValueError("barcode mismatch")
    X=mmread(str(BASE/"matrix.mtx.gz")).tocsr().T.astype(np.float32)
    lib=np.asarray(X.sum(axis=1)).ravel(); scale=np.divide(1e4,lib,out=np.zeros_like(lib),where=lib>0)
    X=X.multiply(scale[:,None]).tocsr(); X.data=np.log1p(X.data)
    mean=np.asarray(X.mean(axis=0)).ravel(); mean2=np.asarray(X.power(2).mean(axis=0)).ravel(); var=mean2-mean**2
    valid=np.array([not (g.startswith('MT-') or g.startswith('RPL') or g.startswith('RPS')) for g in genes])
    idx=np.where(valid)[0][np.argsort(var[valid])[-2000:]]; Xh=X[:,idx]
    pcs=TruncatedSVD(n_components=30,random_state=42).fit_transform(Xh)
    emb=umap.UMAP(n_neighbors=20,min_dist=0.30,n_components=2,metric="euclidean",random_state=42,n_jobs=1).fit_transform(pcs)
    pde_idx=genes.index('PDE8B'); pde=np.asarray(X[:,pde_idx].todense()).ravel()
    out=pd.DataFrame({'barcode':barcodes,'UMAP1':emb[:,0],'UMAP2':emb[:,1],'PDE8B_log_cp10k':pde,
                      'cell_type':meta.Cell_annotation.astype(str).values,'compartment':meta.Compartments.astype(str).values,
                      'sample':meta['new.ident'].astype(str).values,'disease':meta.disease.astype(str).values})
    out.to_csv(OUT/"gse210248_recomputed_umap.csv.gz",index=False,compression='gzip')
    params={'cells':X.shape[0],'genes':X.shape[1],'hvg':2000,'pca_components':30,'n_neighbors':20,'min_dist':0.30,'random_state':42,'note':'UMAP recomputed from public raw counts; not the authors original embedding.'}
    (OUT/"umap_parameters.json").write_text(json.dumps(params,indent=2),encoding='utf-8'); print(json.dumps(params,indent=2))
if __name__=='__main__': main()
