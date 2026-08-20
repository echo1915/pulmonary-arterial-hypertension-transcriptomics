\name{pahscrna-package}
\alias{pahscrna-package}
\alias{pah_default_config}
\alias{pah_project_paths}
\alias{pah_validate_metadata}
\alias{pah_pseudobulk_counts}
\alias{pah_checkpoint}
\alias{pah_optional_modules}
\title{Guardrailed single-cell utilities for the PAH project}
\description{
Utilities for explicit external paths, subject-level metadata validation,
raw-count pseudobulk aggregation, resumable checkpoints, and optional-module
dependency reporting.
}
\usage{
pah_default_config(data_root = Sys.getenv("PAH_DATA_ROOT", unset =
  file.path("G:/workdata/projects",
  "project-001-pulmonary-arterial-hypertension-transcriptomics")),
  output_root = Sys.getenv("PAH_OUTPUT_ROOT", unset =
  file.path(data_root, "outputs", "scrna-rpkg-validation")), seed = 123L)

pah_project_paths(config = pah_default_config(), create = FALSE)

pah_validate_metadata(metadata, subject_col = "subject_id",
  sample_col = "sample_id", group_col = "group",
  celltype_col = "cell_type", min_subjects_per_group = 3L)

pah_pseudobulk_counts(counts, metadata, subject_col = "subject_id",
  group_col = "group", celltype_col = "cell_type", min_cells = 20L,
  min_subjects_per_group = 3L)

pah_checkpoint(path, code, fingerprint, force = FALSE)

pah_optional_modules()
}
\arguments{
  \item{data_root}{External PAH project data directory.}
  \item{output_root}{External directory for newly generated outputs.}
  \item{seed}{Integer random seed.}
  \item{config}{Configuration returned by \code{pah_default_config()}.}
  \item{create}{Whether to create writable external directories.}
  \item{metadata}{Cell-level metadata data frame.}
  \item{subject_col}{Name of the independent subject identifier column.}
  \item{sample_col}{Name of the sample identifier column.}
  \item{group_col}{Name of the curated biological group column.}
  \item{celltype_col}{Name of the curated cell-type column.}
  \item{min_subjects_per_group}{Minimum independent subjects required per group.}
  \item{counts}{Gene-by-cell matrix of non-negative raw integer counts.}
  \item{min_cells}{Minimum cells retained in each pseudobulk library.}
  \item{path}{Destination RDS checkpoint path.}
  \item{code}{Expression evaluated when a valid checkpoint is unavailable.}
  \item{fingerprint}{Non-empty input and parameter identity string.}
  \item{force}{Whether to ignore an existing checkpoint.}
}
\details{
Phenotype groups are never guessed. Confirmatory comparisons must use subjects,
not cells, as independent biological replicates. Optional packages are reported
but never installed automatically.
}
\value{
Configuration and path functions return named lists or vectors. Metadata
validation invisibly returns a group audit. Pseudobulk aggregation returns a
list with counts, library metadata, and exclusions. The checkpoint function
returns the cached or evaluated value. Optional-module reporting returns a data
frame.
}
\examples{
cfg <- pah_default_config()
pah_project_paths(cfg)
pah_optional_modules()
}
