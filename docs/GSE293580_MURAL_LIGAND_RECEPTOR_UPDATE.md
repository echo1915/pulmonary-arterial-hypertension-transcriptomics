# GSE293580 mural ligand-receptor proxy update

Date: 2026-08-10

## Project identity

- Project number: `project-001`
- Canonical project name: `project-001-pulmonary-arterial-hypertension-transcriptomics`
- Code directory: `C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- External data directory: `G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`

## Analysis definition

Seven prespecified axes were audited: PDGF, TGF-beta, JAG-NOTCH3, EDN,
CXCL12-CXCR4, ECM-integrin, and THBS1-CD47. All component genes were present in
the GSE293580 matrix.

For each patient, ligand and receptor genes were aggregated separately in SMCs
and pericytes. A directional compatibility score was defined as the geometric
mean of the sender ligand-module and receiver receptor-module log1p CPM. Both
SMC-to-pericyte and pericyte-to-SMC directions were tested. Each compartment
required at least five cells.

The primary SMC boundary excluded `SMC_like` and combined
`contractile_SMC` with `modulated_SMC`. A secondary boundary included
`SMC_like`. Patients/samples, not cells, were the statistical units.

## Combined PAH results

For core SMCs, the largest positive compatibility effects were:

- SMC to pericyte TGF-beta: g = 1.22, exact P = 0.179.
- SMC to pericyte ECM-integrin: g = 0.99, exact P = 0.143.
- SMC to pericyte PDGF: g = 0.86, exact P = 0.464.
- Pericyte to SMC CXCL12-CXCR4: g = 1.22, exact P = 0.321.
- Pericyte to SMC JAG-NOTCH3: g = 0.97, exact P = 0.143.
- Pericyte to SMC ECM-integrin: g = 0.72, exact P = 0.321.

SMC-to-pericyte JAG-NOTCH3 compatibility decreased (g = -1.73, exact
P = 0.071). All proxy-score BH-adjusted q values were at least 0.417.

All fourteen direction-by-axis effect signs were unchanged when `SMC_like`
was added to the SMC boundary.

## Component audit

The proxy scores should not be interpreted without their components. Notable
component effects included increased SMC-side TGF-beta ligand expression,
increased pericyte-side CXCL12 ligand expression, increased SMC-side integrin
receptor expression, and decreased SMC-side NOTCH3 receptor expression. No
component family survived multiplicity adjustment.

The positive pericyte-to-SMC JAG-NOTCH3 proxy therefore reflects increased
pericyte ligand expression despite lower SMC NOTCH3 receptor expression. It is
not evidence of uniformly activated Notch signaling.

## Disease-subtype heterogeneity

ECM-integrin compatibility was concentrated in SSc-PAH:

- SMC to pericyte: g = 5.08, exact P = 0.067, leave-one-out direction 6/6.
- Pericyte to SMC: g = 2.31, exact P = 0.067, leave-one-out direction 6/6.

The corresponding IPAH effects were close to zero. This supports an
SSc-PAH-specific matrix-remodeling hypothesis rather than a universal PAH
communication claim.

## Inference boundary

These are expression-compatibility proxies. They do not establish spatial
proximity, ligand secretion, receptor protein abundance, pathway activation,
or causal cell-cell communication. Only two donors and six eligible PAH
patients contributed to the combined primary comparisons.

## Outputs

Outputs are under `outputs/gse293580-reanalysis/ligand-receptor-proxy/`
relative to the external data directory:

- `gse293580_mural_lr_gene_patient_pseudobulk.csv`
- `gse293580_mural_lr_patient_proxy_scores.csv`
- `gse293580_mural_lr_patient_tests.csv`
- `gse293580_mural_lr_component_tests.csv`
- `gse293580_mural_lr_summary.json`
