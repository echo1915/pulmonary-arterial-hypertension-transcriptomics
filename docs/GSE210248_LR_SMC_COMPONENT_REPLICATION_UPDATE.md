# GSE210248 ligand-receptor SMC-component replication update

Date: 2026-08-10

## Project identity

- Project number: `project-001`
- Canonical project name: `project-001-pulmonary-arterial-hypertension-transcriptomics`
- Code directory: `C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- External data directory: `G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`

## Feasibility boundary

GSE210248 contains SMC1 and SMC2 annotations but no independent pericyte
annotation. Fibroblasts were not relabeled as pericytes. Consequently, this
dataset cannot independently replicate SMC-pericyte communication.

The analysis was restricted to the SMC-side molecular components of the seven
prespecified GSE293580 ligand-receptor axes. SMC ligand modules represent the
SMC-to-pericyte model side, and SMC receptor modules represent the
pericyte-to-SMC model side.

Raw counts were streamed from the original sparse matrix and aggregated by
patient and SMC subtype. Statistical comparisons used three PAH and three donor
patients; cells were not treated as independent replicates.

## Results

Ligand-module directions were more reproducible than receptor-module
directions:

- 12 of 14 ligand axis-by-SMC-subtype comparisons were directionally
  concordant between GSE210248 and GSE293580.
- 6 of 14 receptor comparisons were concordant.
- Overall concordance was 18 of 28 comparisons (64.3%).

PDGF, TGF-beta, JAG, ECM, and THBS1 ligand modules retained the GSE293580
direction in both SMC1 and SMC2. EDN and CXCL12 ligand directions differed in
SMC2.

TGF-beta, EDN, and ECM-integrin receptor directions were concordant in both SMC
subtypes. PDGFRB, NOTCH3, CXCR4, and CD47 receptor directions were discordant.
This argues against treating the full receptor-side communication model as a
stable cross-cohort program.

The strongest GSE210248 component was increased CXCR4 in SMC1 (g = 2.49,
exact P = 0.10), but this was not directionally concordant with GSE293580 and
did not survive multiplicity adjustment. No GSE210248 component family passed
BH correction.

## Interpretation

The most defensible cross-cohort conclusion is that several SMC ligand and
matrix-remodeling programs show directional replication. The complete
SMC-pericyte communication architecture does not replicate because GSE210248
lacks a pericyte population and receptor-side directions are heterogeneous.

These results support prioritizing ECM-integrin and TGF-beta for follow-up,
while treating NOTCH3, PDGFRB, and CXCR4 receptor behavior as state- or
cohort-dependent.

## Outputs

Outputs are under `outputs/gse210248-lr-smc-component-replication/` relative to
the external data directory:

- `gse210248_lr_gene_smc_patient_pseudobulk.csv`
- `gse210248_lr_smc_patient_module_scores.csv`
- `gse210248_lr_smc_component_effects.csv`
- `gse210248_gse293580_lr_smc_component_concordance.csv`
- `gse210248_lr_smc_component_summary.json`
