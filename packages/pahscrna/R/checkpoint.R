#' Run an expression with a versioned RDS checkpoint
#' @param path RDS path under an external output directory.
#' @param code Expression to evaluate when no valid checkpoint exists.
#' @param fingerprint String identifying inputs and parameters.
#' @param force Ignore an existing checkpoint.
#' @return The computed or cached value.
#' @export
pah_checkpoint <- function(path, code, fingerprint, force = FALSE) {
  if (!is.character(fingerprint) || length(fingerprint) != 1L || !nzchar(fingerprint)) {
    stop("fingerprint must be one non-empty string", call. = FALSE)
  }
  if (!force && file.exists(path)) {
    saved <- readRDS(path)
    if (is.list(saved) && identical(saved$fingerprint, fingerprint) && "value" %in% names(saved)) {
      return(saved$value)
    }
  }
  value <- eval.parent(substitute(code))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "checkpoint-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(list(fingerprint = fingerprint, created_at = Sys.time(), value = value), tmp)
  if (!file.rename(tmp, path)) stop("Could not atomically install checkpoint: ", path, call. = FALSE)
  value
}

#' Report availability of optional analysis modules
#' @return Data frame; this function never installs packages.
#' @export
pah_optional_modules <- function() {
  packages <- c(
    cellchat = "CellChat", pseudotime_monocle3 = "monocle3",
    pseudotime_slingshot = "slingshot", hdwgcna = "hdWGCNA",
    differential_expression = "edgeR"
  )
  data.frame(
    module = names(packages), package = unname(packages),
    available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    row.names = NULL
  )
}

