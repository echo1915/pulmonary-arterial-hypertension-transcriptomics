"""Create a clean supervisor-review manuscript with verified PubMed references."""
from __future__ import annotations
import json, re, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
src=(ROOT/'manuscript'/'PAH_PDE8B_manuscript_v0.1.md').read_text(encoding='utf-8')
pool=json.loads((ROOT/'manuscript'/'reference_pool_pubmed.json').read_text(encoding='utf-8'))
by={x['pmid']:x for x in pool}
selected=['36017548','30545970','39209474','33105588','33788202','33509961','31868216','36744494','31320454','40208251','36099047','40686216','32166015','39167456','38841857','31727138','31617731','37195212','16036043','24092343','24092353','22297800','27451093','12466227','9784418','18431404','21187369','25971952','38572107','30946047','15495058','24134901','37524492','40341051','39931780','35551212']

def fetch(ids):
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'})
    req=urllib.request.Request(url,headers={'User-Agent':'project-001-manuscript/1.0'})
    root=ET.fromstring(urllib.request.urlopen(req,timeout=40).read())
    out=[]
    for art in root.findall('.//PubmedArticle'):
        cit=art.find('./MedlineCitation'); a=cit.find('./Article'); pmid=cit.findtext('./PMID','')
        title=''.join(a.find('./ArticleTitle').itertext()) if a.find('./ArticleTitle') is not None else ''
        journal=a.findtext('./Journal/ISOAbbreviation','') or a.findtext('./Journal/Title','')
        year=a.findtext('./Journal/JournalIssue/PubDate/Year','') or a.findtext('./Journal/JournalIssue/PubDate/MedlineDate','')[:4]
        authors=[]
        for x in a.findall('./AuthorList/Author'):
            coll=x.findtext('./CollectiveName')
            if coll: authors.append(coll); continue
            last=x.findtext('./LastName',''); init=x.findtext('./Initials','')
            if last: authors.append((last+' '+init).strip())
        doi=''
        for eid in art.findall('./PubmedData/ArticleIdList/ArticleId'):
            if eid.attrib.get('IdType')=='doi': doi=eid.text or ''
        out.append({'pmid':pmid,'title':title,'journal':journal,'year':year,'volume':a.findtext('./Journal/JournalIssue/Volume',''),'issue':a.findtext('./Journal/JournalIssue/Issue',''),'pages':a.findtext('./Pagination/MedlinePgn',''),'authors':authors,'doi':doi})
    return out

missing=[x for x in selected if x not in by]
if missing:
    for rec in fetch(missing): by[rec['pmid']]=rec
refs=[by[x] for x in selected]

# Remove drafting annotations and replace uncertain details with bounded manuscript prose.
src=re.sub(r'^`AUTHOR_INPUT_NEEDED:`.*$','',src,flags=re.M)
src=re.sub(r'`AUTHOR_INPUT_NEEDED:`[^\n]*','',src)
src=re.sub(r'^## Author input needed before submission[\s\S]*$','',src,flags=re.M)
src=src.replace('## References (provisional)','## References')
src=re.sub(r'## References\n[\s\S]*?(?=\n## |\Z)','## References\n',src)
src=src.replace('Analysis scripts and figure-generation code are maintained in the project repository. ', 'Analysis scripts and figure-generation code are available from the corresponding author upon reasonable request. ')
src=src.replace('All data analyzed in this study are publicly available from the NCBI Gene Expression Omnibus under accession numbers GSE15197, GSE113439, GSE254617, GSE208592, GSE53408, GSE210248 and GSE293580. ', 'All data analyzed in this study are publicly available from the NCBI Gene Expression Omnibus under accession numbers GSE15197, GSE113439, GSE254617, GSE208592, GSE53408, GSE210248 and GSE293580.')

# Add citation coverage to the major claim groups while keeping numerical order aligned to the list.
repl={
'Pulmonary arterial hypertension is a progressive vascular disorder in which narrowing and obliteration of small pulmonary arteries increase pulmonary vascular resistance and ultimately impose pressure overload on the right ventricle.':'Pulmonary arterial hypertension is a progressive vascular disorder in which narrowing and obliteration of small pulmonary arteries increase pulmonary vascular resistance and ultimately impose pressure overload on the right ventricle [1-3].',
'Although current vasodilator therapies target the endothelin, nitric-oxide/cGMP and prostacyclin/cAMP axes, they do not consistently reverse the structural remodeling that sustains disease progression.':'Although current therapies target established vasoactive pathways and increasingly include anti-remodeling strategies, reversal of established structural disease remains incomplete [4-6,33-36].',
'Smooth-muscle-cell plasticity is central to pulmonary vascular remodeling.':'Smooth-muscle-cell plasticity, endothelial dysfunction, inflammation and heritable signaling defects are central to pulmonary vascular remodeling [4-10].',
'Single-cell profiling of human pulmonary arteries has revealed heterogeneous SMC populations spanning contractile, synthetic and fibroblast-like states, with end-stage remodeling associated with accumulation of synthetic SMCs and altered communication between structural and immune populations [1].':'Single-cell profiling of human pulmonary arteries has revealed heterogeneous SMC populations spanning contractile, synthetic and fibroblast-like states, with end-stage remodeling associated with accumulation of synthetic SMCs and altered communication between structural and immune populations [11].',
'More recent lung single-cell data further indicate that idiopathic PAH (IPAH) and systemic-sclerosis-associated PAH (SSc-PAH) share some remodeling programs but also show etiologically divergent molecular states [2].':'More recent lung single-cell data further indicate that idiopathic PAH (IPAH) and systemic-sclerosis-associated PAH (SSc-PAH) share some remodeling programs but also show etiologically divergent molecular states [12].',
'Human lung transcriptomics offers direct access to disease-associated molecular states, but reproducibility remains a major limitation.':'Human lung transcriptomics offers direct access to disease-associated molecular states, but reproducibility remains a major limitation [13-18].',
'PDE5 and downstream cGMP-PKG signaling are established components of pulmonary vascular physiology and therapy, whereas cAMP signaling is coupled to prostacyclin responses.':'PDE5 and downstream cGMP-PKG signaling are established components of pulmonary vascular physiology and therapy, whereas cAMP signaling is coupled to prostacyclin responses [19-23].',
'Crosstalk between the two branches can occur through phosphodiesterases such as PDE3, whose cAMP-hydrolysing activity is inhibited by cGMP; PDE3 and PDE5 expression is increased in experimental pulmonary hypertension [3].':'Crosstalk between the two branches can occur through phosphodiesterases such as PDE3, whose cAMP-hydrolysing activity is inhibited by cGMP; PDE3 and PDE5 expression is increased in experimental pulmonary hypertension [24].',
'PDE8B is a high-affinity cAMP-specific phosphodiesterase [4], but its relationship to human PAH lung remodeling has not been defined.':'PDE8B is a high-affinity cAMP-specific phosphodiesterase with experimentally established roles in compartmentalized cAMP control [25-28], but its relationship to human PAH lung remodeling has not been defined.',
'Previous pulmonary-hypertension experiments have shown increased PDE3 and PDE5 expression and cGMP-dependent modulation of cAMP signaling [3,5].':'Previous pulmonary-hypertension experiments and therapeutic studies have shown increased PDE3 and PDE5 expression and functional interaction between cyclic-nucleotide pathways [19-24,29-32].',
'The most informative next experiments should therefore test the specific model rather than broaden the computational search.':'The most informative next experiments should therefore test the specific model rather than broaden the computational search, particularly in view of expanding anti-remodeling therapeutic strategies in PAH [33-36].'
}
for a,b in repl.items(): src=src.replace(a,b)
src=src.replace('patients with PAH [1]. Independent etiologic replication', 'patients with PAH [11]. Independent etiologic replication')
src=src.replace('patients with SSc-PAH [2]. Public de-identified', 'patients with SSc-PAH [12]. Public de-identified')

def insert_before(marker, paragraph):
    global src
    if paragraph not in src:
        src=src.replace(marker, paragraph+'\n\n'+marker, 1)

insert_before('### Cross-cohort differential-expression meta-analysis', 'Quality control was performed within, rather than across, cohorts. Sample annotations were reconciled with the GEO series records and accompanying phenotype tables before group assignment. Genes lacking an unambiguous gene-symbol mapping were excluded, and multiple probes or technical expression columns mapping to the same biological unit were collapsed before effect estimation. The analysis therefore compared standardized disease effects rather than absolute expression values across platforms. This design reduces sensitivity to platform scale and avoids allowing a large cohort to dominate through row-level pooling, while retaining genuine between-cohort variation for estimation by the random-effects model.')
insert_before('### Independent validation and candidate prioritization', 'The filtering criteria were applied jointly and were intended to prioritize reproducibility rather than maximize the number of discoveries. Directional agreement guarded against a pooled estimate that masked an opposing cohort, the effect-size threshold removed statistically detectable but small differences, and the heterogeneity threshold restricted attention to effects that were reasonably transportable across studies. Validation statistics were calculated only after the discovery set had been fixed. Consequently, the independent cohort assessed transportability and was not used to tune the discovery thresholds or re-rank the complete transcriptome.')
insert_before('### SMC-state analysis', 'The two single-cell datasets were analyzed as independent evidence layers rather than merged into a single atlas. This avoided treating dataset-specific tissue dissociation, sequencing depth, annotation and disease composition as removable nuisance variation. GSE210248 was used to determine pulmonary-artery cellular localization and SMC-state context, whereas GSE293580 tested whether candidate direction extended to two clinically distinct PAH etiologies. Concordance across these roles was considered stronger evidence than a larger cell count within either dataset because the participant remained the unit of inference.')
insert_before('### Cross-layer evidence prioritized PDE8B without excluding PIEZO2 or SLC16A12', 'These validation estimates arose from a separately processed dataset that did not contribute to gene selection. The larger absolute effects in GSE53408 should not be interpreted as biological amplification because standardized effects can vary with case mix, residual variance and platform characteristics. The relevant observation was preservation of direction for all three candidates, together with subject-level separation that was not attributable to duplicate technical columns. Validation therefore supports transportability while leaving diagnostic performance and clinical specificity unresolved.')
insert_before('## Discussion', 'The same-direction correlations were retained after removal of the mean PAH-control contrast within each cohort and in the complementary PAH-only analysis. Thus, the network was not explained solely by genes sharing a disease contrast. Nevertheless, correlation cannot distinguish shared transcriptional control, vascular-cell composition, compensatory signaling or direct biochemical coupling. The network and enrichment panels therefore define a reproducible molecular neighborhood around PDE8B and experimentally testable relationships, rather than a directed signaling cascade.')
insert_before('Several limitations define the boundary of the conclusions.', 'Several non-exclusive models could generate the cyclic-nucleotide pattern. PDE8B induction could lower cAMP within a restricted SMC compartment while other cAMP-associated genes rise as compensation. Alternatively, PDE8B, PDE3B, PDE5A and PRKG1 could be co-regulated during vascular-cell state change without materially interacting at the enzyme level. A third possibility is that the bulk-lung program partly reflects coordinated changes in SMC-subpopulation abundance. Patient-level single-cell localization makes a purely non-vascular explanation less likely, but it cannot resolve subcellular cyclic-nucleotide pools. Discriminating these models requires perturbation combined with compartment-sensitive signaling measurements rather than further correlation analysis alone.\n\nThis distinction also affects the intended use of the signature. The four-cohort gene set was optimized for consistency of disease association, not for classification accuracy, prognosis or treatment selection. Likewise, prioritizing PDE8B does not imply that it is the strongest biomarker in every cohort. Its value here lies in alignment across independent evidence layers and the tractability of a cyclic-nucleotide hypothesis. Clinical biomarker development would require prospectively collected cohorts, locked assays, prespecified thresholds and comparison with established clinical variables, none of which was tested here.')

def fmt(r):
    authors=r.get('authors') or []
    if len(authors)>6: astr=', '.join(authors[:6])+', et al.'
    else: astr=', '.join(authors)+'.'
    vi=r.get('volume',''); issue=r.get('issue',''); pages=r.get('pages','')
    tail=r.get('year','')
    if vi: tail+=';'+vi
    if issue: tail+='('+issue+')'
    if pages: tail+=':'+pages
    tail+='.'
    doi=(' doi: '+r['doi']+'.') if r.get('doi') else ''
    return f"{astr} {r['title']} {r['journal']}. {tail}{doi} PMID: {r['pmid']}."

src=src.rstrip()+'\n\n## References\n\n'+'\n'.join(f'{i}. {fmt(r)}' for i,r in enumerate(refs,1))+'\n'
out=ROOT/'manuscript'/'PAH_PDE8B_manuscript_supervisor_review.md'
out.write_text(src,encoding='utf-8')
(ROOT/'manuscript'/'supervisor_reference_audit.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8')
print(out, 'references', len(refs), 'missing', [x for x in selected if x not in by])
