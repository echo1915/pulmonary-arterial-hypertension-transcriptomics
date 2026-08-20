# Data

Large and generated project data are external to this repository:

`G:\workdata\projects\project-001-project-001-pulmonary-arterial-hypertension-transcriptomics`

- `raw/current-data-snapshot/`: preserved snapshot of all inputs used before migration.
- `interim/`: temporary and transformed working data.
- `processed/`: standardized analysis-ready data.
- `outputs/current-results/`: preserved result set at the migration checkpoint.
- `manifests/`: source and integrity records.

Override the location with `PAH_DATA_ROOT`. Do not place GEO matrices or scRNA-seq
objects in this repository.

