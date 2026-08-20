"""Fetch a deduplicated PubMed reference pool for the supervisor manuscript."""
from __future__ import annotations
import json, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

QUERIES = [
    '(pulmonary arterial hypertension[Title/Abstract]) AND (guideline OR pathobiology OR vascular remodeling)',
    '(pulmonary arterial hypertension[Title/Abstract]) AND (smooth muscle cell OR phenotypic switching OR single-cell RNA sequencing)',
    '(pulmonary hypertension[Title/Abstract]) AND (cAMP OR cGMP OR phosphodiesterase OR PDE3 OR PDE5)',
    '(PDE8B[Title/Abstract]) OR (phosphodiesterase 8B[Title/Abstract])',
    '(pulmonary arterial hypertension[Title/Abstract]) AND (transcriptomics OR gene expression OR meta-analysis)',
]
BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def get(url):
    req=urllib.request.Request(url,headers={'User-Agent':'project-001-manuscript/1.0'})
    return urllib.request.urlopen(req,timeout=40).read()

ids=[]
for query in QUERIES:
    url=BASE+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':query,'retmax':25,'sort':'relevance','retmode':'json'})
    data=json.loads(get(url))
    ids.extend(data['esearchresult']['idlist'])
ids=list(dict.fromkeys(ids))

url=BASE+'efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'})
root=ET.fromstring(get(url))
records=[]
for art in root.findall('.//PubmedArticle'):
    cit=art.find('./MedlineCitation'); a=cit.find('./Article')
    pmid=cit.findtext('./PMID','')
    title=''.join(a.find('./ArticleTitle').itertext()) if a.find('./ArticleTitle') is not None else ''
    journal=a.findtext('./Journal/Title','')
    year=a.findtext('./Journal/JournalIssue/PubDate/Year','') or a.findtext('./Journal/JournalIssue/PubDate/MedlineDate','')[:4]
    volume=a.findtext('./Journal/JournalIssue/Volume',''); issue=a.findtext('./Journal/JournalIssue/Issue','')
    pages=a.findtext('./Pagination/MedlinePgn','')
    authors=[]
    for x in a.findall('./AuthorList/Author'):
        coll=x.findtext('./CollectiveName')
        if coll: authors.append(coll); continue
        last=x.findtext('./LastName',''); init=x.findtext('./Initials','')
        if last: authors.append((last+' '+init).strip())
    doi=''
    for eid in art.findall('./PubmedData/ArticleIdList/ArticleId'):
        if eid.attrib.get('IdType')=='doi': doi=eid.text or ''
    records.append({'pmid':pmid,'title':title,'journal':journal,'year':year,'volume':volume,'issue':issue,'pages':pages,'authors':authors,'doi':doi})

out=Path(__file__).resolve().parents[1]/'manuscript'/'reference_pool_pubmed.json'
out.write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
print(out, len(records))
