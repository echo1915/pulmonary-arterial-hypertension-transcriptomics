# GSE293580 Phase 4 patient-level interaction update

Date: 2026-08-10

## Project identity

- Project number: `project-001`
- Canonical project name: `project-001-pulmonary-arterial-hypertension-transcriptomics`
- Code directory: `C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- External data directory: `G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`

## Continuation point

This analysis continued from the locked Phase 3 mural-cell annotations and the
existing patient-level SMC program checkpoint. It did not repeat clustering,
annotation, or the completed patient pseudobulk analysis.

Statistical units are patients/samples. A compartment requires at least five
cells. All tests are exploratory because only two donor samples have eligible
paired SMC and pericyte observations.

## Added analyses

1. Disease-by-cell-state interaction: exact label-permutation comparisons of
   each patient's SMC-minus-pericyte PDE8B detection contrast against donors.
2. Patient-level SMC program comparisons for contractile, matrix-remodeling,
   proliferation, and cAMP-response programs.
3. Leave-one-patient-out ranges for PDE8B-program Spearman correlations.

## Results and interpretation

- The combined PAH SMC-minus-pericyte contrast was 0.079 higher than donors
  (Hedges g = 0.455, exact P = 0.536). IPAH showed essentially no interaction
  shift (difference-in-differences = 0.004, exact P = 1.000).
- In combined PAH, the contractile program was lower (g = -0.842) and the
  matrix-remodeling program was higher (g = 0.713), but neither survived exact
  testing or within-contrast BH adjustment (all q >= 0.611).
- PAH-only PDE8B-program correlations were weak. Their leave-one-patient-out
  ranges crossed zero for every module, so stable within-PAH coupling was not
  demonstrated.

These results are directionally compatible with SMC phenotypic remodeling but
do not establish a disease-specific SMC-pericyte interaction or a stable PDE8B
state-program relationship in GSE293580. They support association/localization
language only, not causality.

## Reproducible outputs

- `outputs/gse293580-reanalysis/gse293580_pde8b_disease_by_cellstate_interaction.csv`
- `outputs/gse293580-reanalysis/gse293580_smc_program_group_tests.csv`
- `outputs/gse293580-reanalysis/gse293580_pde8b_smc_program_correlations.csv`
- `outputs/gse293580-reanalysis/gse293580_phase4_state_interaction_summary.json`

The output paths above are relative to the external data directory.
