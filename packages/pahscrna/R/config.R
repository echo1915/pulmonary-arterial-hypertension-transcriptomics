#' Create a conservative PAH single-cell configuration
#'
#' Defaults are intentionally analysis-safe: automatic phenotype guessing is
#' disabled and subject-level replication is required.
#' @param data_root External project data directory.
#' @param output_root Directory for newly generated outputs.
#' @param seed Random seed.
#' @return A named list.
#' @export
pah_default_config <- function(
    data_root = Sys.getenv("PAH_DATA_ROOT", unset = file.path(
      "G:/workdata/projects",
      "project-001-pulmonary-arterial-hypertension-transcriptomics"
    )),
    output_root = Sys.getenv("PAH_OUTPUT_ROOT", unset = file.path(
      data_root, "outputs", "scrna-rpkg-validation"
    )),
    seed = 123L) {
  list(
    data_root = normalizePath(data_root, winslash = "/", mustWork = FALSE),
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    seed = as.integer(seed),
    subject_col = "subject_id",
    sample_col = "sample_id",
    group_col = "group",
    celltype_col = "cell_type",
    min_subjects_per_group = 3L,
    infer_groups = FALSE,
    run_pseudobulk = TRUE,
    run_cellchat = FALSE,
    run_pseudotime = FALSE,
    run_hdwgcna = FALSE,
    force_recompute = FALSE
  )
}

#' Resolve project paths without writing to the repository
#' @param config Output from [pah_default_config()].
#' @param create Create writable external output folders.
#' @return Named character vector.
#' @export
pah_project_paths <- function(config = pah_default_config(), create = FALSE) {
  paths <- c(
    raw = file.path(config$data_root, "raw"),
    interim = file.path(config$data_root, "interim"),
    processed = file.path(config$data_root, "processed"),
    output = config$output_root,
    checkpoints = file.path(config$output_root, "checkpoints"),
    logs = file.path(config$output_root, "logs")
  )
  paths <- vapply(paths, normalizePath, character(1), winslash = "/", mustWork = FALSE)
  if (isTRUE(create)) {
    writable <- paths[names(paths) != "raw"]
    invisible(vapply(writable, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  }
  paths
}

