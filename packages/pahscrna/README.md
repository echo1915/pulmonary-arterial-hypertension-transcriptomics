# pahscrna

`pahscrna` is a conservative extraction of reusable ideas from the WeChat
article “一键跑完GEO单细胞分析：从Seurat到CellChat、拟时序、hdWGCNA全流程自动化”.
It is not a verbatim repackaging of the one-click script.

The package deliberately makes subject IDs, curated phenotype groups, raw-count
pseudobulk aggregation, external output paths, and checkpoint fingerprints
explicit. It never installs packages or guesses case/control labels.

```r
devtools::load_all("packages/pahscrna")
cfg <- pah_default_config()
paths <- pah_project_paths(cfg)
pah_validate_metadata(cell_metadata)
pb <- pah_pseudobulk_counts(raw_counts, cell_metadata)
```

On this workstation, run R through the project wrapper so the external package
library and ASCII-only temporary directory are configured consistently:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_project_r.ps1 -Expression "library(pahscrna); pah_optional_modules()"
```

CellChat, pseudotime, and hdWGCNA are optional exploratory modules. Their
availability can be inspected with `pah_optional_modules()`; they are disabled
in the default configuration and must not replace subject-level inference.
