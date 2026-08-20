# Pulmonary Arterial Hypertension Transcriptomics

This project continues the existing PAH multi-cohort lung transcriptomics and
single-cell analysis. The current article direction is a robust pulmonary
vascular-remodeling signature centered on an SMC-associated PDE8B/cAMP axis,
with an endothelial ABCG2/DCBLD1/SELP module as supporting context.

## Current status

- Four independent human lung bulk cohorts completed for random-effects meta-analysis.
- Disease-specific comparisons against IPF-PH/PHPF-PH completed.
- Pulmonary-artery scRNA-seq localization completed in GSE210248.
- Independent whole-lung replication completed in GSE293580.
- PDE8B is the current primary article candidate; causality is not yet established.

See [docs/PROJECT_STATUS.md](docs/PROJECT_STATUS.md) for the migration checkpoint and
`docs/PAH_SCHEME2_FINAL_ARTICLE_ASSESSMENT.md` for the current scientific assessment.

## Layout

- `scripts/`: reproducible audit and analysis scripts retained from the completed work.
- `docs/`: protocols, decisions, and article assessments.
- `config/`: project configuration documentation.
- `data/`: documentation and optional small fixtures only.
- Large data and generated outputs are kept outside the repository and resolved through `scripts/project_paths.py`.

## Path configuration

Defaults are defined in `scripts/project_paths.py`.

```powershell
$env:WORKDATA_ROOT = 'G:\workdata'
# Optional full override:
$env:PAH_DATA_ROOT = 'G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics'
# Use a separate output location for validation runs:
$env:PAH_OUTPUT_ROOT = 'G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\outputs\validation'
```

## Continue the analysis

Run commands from the repository root. The bundled Codex Python runtime used for
the completed analyses includes pandas and NumPy.

```powershell
python scripts\pah\lung_mechanism_meta.py
python scripts\pah\scrna_candidate_localization.py
python scripts\pah\scrna_ph_specificity.py
python scripts\pah\scrna_gse293580_replication.py
```

The PDE8B cross-cohort coexpression, pathway, subject-level sensitivity and
manuscript-figure analyses are complete. Before rerunning the full workflow,
validate the configured external paths:

```powershell
python tests\test_project_paths.py
```

The Respiratory Research submission revision is documented in
[`docs/RESPIRATORY_RESEARCH_REVISION_AND_RELEASE_PLAN.md`](docs/RESPIRATORY_RESEARCH_REVISION_AND_RELEASE_PLAN.md).
Public release requirements are tracked in [`CODE_RELEASE_CHECKLIST.md`](CODE_RELEASE_CHECKLIST.md).

## Project R package

[`packages/pahscrna`](packages/pahscrna) contains guardrailed single-cell
utilities adapted from a published GEO/Seurat workflow. It validates curated
subject metadata, builds raw-count pseudobulk matrices, and manages external
checkpoints. Optional CellChat, pseudotime, and hdWGCNA modules are disabled by
default and do not replace the completed subject-level PAH analyses.


