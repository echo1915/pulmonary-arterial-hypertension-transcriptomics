#' Validate subject-level single-cell metadata
#' @param metadata Data frame with one row per cell.
#' @param subject_col Subject identifier column.
#' @param sample_col Sample identifier column.
#' @param group_col Biological comparison column.
#' @param celltype_col Cell-type annotation column.
#' @param min_subjects_per_group Minimum independent subjects in every group.
#' @return Invisibly returns a compact audit table.
#' @export
pah_validate_metadata <- function(
    metadata,
    subject_col = "subject_id",
    sample_col = "sample_id",
    group_col = "group",
    celltype_col = "cell_type",
    min_subjects_per_group = 3L) {
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame", call. = FALSE)
  required <- c(subject_col, sample_col, group_col, celltype_col)
  missing_cols <- setdiff(required, names(metadata))
  if (length(missing_cols)) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!nrow(metadata)) stop("metadata has no cells", call. = FALSE)
  bad <- vapply(metadata[required], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(bad)) stop("Missing/blank values in: ", paste(names(bad)[bad], collapse = ", "), call. = FALSE)

  subject_group <- unique(metadata[c(subject_col, group_col)])
  split_groups <- split(subject_group[[group_col]], subject_group[[subject_col]])
  crossed <- names(Filter(function(x) length(unique(x)) != 1L, split_groups))
  if (length(crossed)) {
    stop("Subjects assigned to multiple groups: ", paste(crossed, collapse = ", "), call. = FALSE)
  }

  sample_subject <- unique(metadata[c(sample_col, subject_col)])
  split_subjects <- split(sample_subject[[subject_col]], sample_subject[[sample_col]])
  ambiguous <- names(Filter(function(x) length(unique(x)) != 1L, split_subjects))
  if (length(ambiguous)) {
    stop("Samples assigned to multiple subjects: ", paste(ambiguous, collapse = ", "), call. = FALSE)
  }

  n_subjects <- table(subject_group[[group_col]])
  if (any(n_subjects < min_subjects_per_group)) {
    detail <- paste(names(n_subjects), as.integer(n_subjects), sep = "=", collapse = ", ")
    stop("Insufficient subject-level replication (need >= ", min_subjects_per_group,
         " per group): ", detail, call. = FALSE)
  }
  audit <- data.frame(group = names(n_subjects), subjects = as.integer(n_subjects), row.names = NULL)
  invisible(audit)
}

