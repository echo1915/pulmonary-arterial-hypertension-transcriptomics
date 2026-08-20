# GSE293580 sensitivity and cross-cohort update

Date: 2026-08-10

## Project identity

- Project number: `project-001`
- Canonical project name: `project-001-pulmonary-arterial-hypertension-transcriptomics`
- Code directory: `C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- External data directory: `G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`

## Continuation strategy

The analysis reused the locked Phase 3 mural annotations, patient PDE8B
pseudobulk table, and patient SMC program table. It did not repeat clustering or
annotation. Patients/samples remained the statistical units.

## Threshold sensitivity

PDE8B expression and SMC-pericyte localization were audited at minimum cell
thresholds of 3, 5, 10, and 20. The patient-program checkpoint was originally
created at a five-SMC minimum, so program sensitivity is identifiable only at
thresholds 5, 10, and 20.

For combined PAH versus donor strict SMC pseudobulk, PDE8B remained higher at
thresholds 3, 5, and 10:

| Minimum SMCs | PAH n | Donor n | Hedges g | Exact P | LOO direction |
|---:|---:|---:|---:|---:|---:|
| 3 | 8 | 2 | 1.55 | 0.133 | 10/10 |
| 5 | 7 | 2 | 2.14 | 0.083 | 9/9 |
| 10 | 6 | 2 | 1.91 | 0.107 | 8/8 |

At threshold 20 no donor retained sufficient strict SMCs, so the contrast is
not estimable. The SMC-minus-pericyte interaction weakened as the threshold
increased (g 0.57, 0.45, and 0.16 at thresholds 3, 5, and 10) and remained
non-significant.

The combined-PAH contractile decrease and matrix-remodeling increase retained
their directions at thresholds 5 and 10. Formal evidence remained limited by
only two eligible donors.

## Locked annotation-boundary sensitivity

At the primary five-cell threshold, excluding `SMC_like` did not remove the
PDE8B signal:

- Core SMC (`contractile_SMC` + `modulated_SMC`): g = 2.05, exact P = 0.083,
  leave-one-out direction 9/9.
- Strict SMC including `SMC_like`: g = 2.14, exact P = 0.083, leave-one-out
  direction 9/9.

Individual-state inference is not reliable: no donor had at least five
`modulated_SMC` cells and only one donor had at least five `SMC_like` cells.

## GSE210248 directional replication

The two cohorts agreed in direction for all prespecified harmonized features:

- PDE8B: higher in GSE210248 SMC1 and SMC2 and in GSE293580 strict SMC.
- Contractile program: lower in both cohorts.
- ECM/matrix-remodeling program: higher in both cohorts.

This is directional replication only. Cell-state definitions and module-score
scales differ between cohorts, and all patient groups are small. The result
supports an association between PDE8B elevation and SMC remodeling but does not
establish causal regulation.

## Outputs

Outputs are under
`outputs/gse293580-reanalysis/sensitivity/` relative to the external data root.
The key files are:

- `gse293580_pde8b_threshold_sensitivity.csv`
- `gse293580_smc_pericyte_interaction_sensitivity.csv`
- `gse293580_smc_program_threshold_sensitivity.csv`
- `gse293580_annotation_boundary_tests.csv`
- `gse210248_gse293580_directional_replication.csv`
- `gse293580_sensitivity_summary.json`
- `gse293580_annotation_boundary_summary.json`
