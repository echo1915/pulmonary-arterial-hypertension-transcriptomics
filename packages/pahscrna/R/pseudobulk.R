#' Aggregate cell-level raw counts to subject-by-cell-type pseudobulk counts
#' @param counts Gene-by-cell numeric matrix, preferably sparse.
#' @param metadata Cell metadata in the same order as columns of `counts`.
#' @param subject_col Subject identifier.
#' @param group_col Group identifier.
#' @param celltype_col Cell-type identifier.
#' @param min_cells Minimum cells in a subject/cell-type pseudobulk library.
#' @param min_subjects_per_group Minimum subjects per group after filtering.
#' @return List containing counts, library metadata, and excluded libraries.
#' @export
pah_pseudobulk_counts <- function(
    counts, metadata,
    subject_col = "subject_id", group_col = "group", celltype_col = "cell_type",
    min_cells = 20L, min_subjects_per_group = 3L) {
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) stop("counts must be matrix-like", call. = FALSE)
  if (ncol(counts) != nrow(metadata)) stop("ncol(counts) must equal nrow(metadata)", call. = FALSE)
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) stop("Gene row names must be unique", call. = FALSE)
  if (any(counts < 0) || any(abs(counts - round(counts)) > .Machine$double.eps^0.5)) {
    stop("pseudobulk requires non-negative raw integer counts", call. = FALSE)
  }
  pah_validate_metadata(metadata, subject_col, subject_col, group_col, celltype_col,
                        min_subjects_per_group)

  key <- interaction(metadata[[subject_col]], metadata[[celltype_col]], drop = TRUE, sep = "||")
  members <- split(seq_len(ncol(counts)), key)
  cell_n <- lengths(members)
  keep <- cell_n >= as.integer(min_cells)
  if (!any(keep)) stop("No pseudobulk libraries pass min_cells", call. = FALSE)
  aggregate_one <- function(i) {
    block <- counts[, i, drop = FALSE]
    if (inherits(block, "Matrix")) Matrix::rowSums(block) else base::rowSums(block)
  }
  bulk <- do.call(cbind, lapply(members[keep], aggregate_one))
  rownames(bulk) <- rownames(counts)
  colnames(bulk) <- names(members)[keep]

  first <- vapply(members[keep], `[`, integer(1), 1L)
  lib_meta <- metadata[first, c(subject_col, group_col, celltype_col), drop = FALSE]
  lib_meta$n_cells <- unname(cell_n[keep])
  rownames(lib_meta) <- colnames(bulk)

  by_ct <- split(lib_meta, lib_meta[[celltype_col]])
  valid_ct <- vapply(by_ct, function(x) {
    all(table(unique(x[c(subject_col, group_col)])[[group_col]]) >= min_subjects_per_group)
  }, logical(1))
  keep_lib <- lib_meta[[celltype_col]] %in% names(valid_ct)[valid_ct]
  if (!any(keep_lib)) stop("No cell type retains sufficient subjects per group", call. = FALSE)
  excluded <- lib_meta[!keep_lib, , drop = FALSE]
  list(counts = bulk[, keep_lib, drop = FALSE], metadata = lib_meta[keep_lib, , drop = FALSE], excluded = excluded)
}
