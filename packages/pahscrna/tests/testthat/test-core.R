test_that("metadata validation enforces subject replication", {
  md <- data.frame(
    subject_id = rep(paste0("S", 1:6), each = 2),
    sample_id = rep(paste0("X", 1:6), each = 2),
    group = rep(c("Control", "PAH"), each = 6),
    cell_type = "SMC"
  )
  expect_silent(pah_validate_metadata(md))
  expect_error(pah_validate_metadata(md[md$subject_id != "S6", ]), "Insufficient")
})

test_that("pseudobulk sums raw counts by subject and cell type", {
  md <- data.frame(
    subject_id = rep(paste0("S", 1:6), each = 2),
    group = rep(c("Control", "PAH"), each = 6),
    cell_type = "SMC"
  )
  x <- matrix(1L, nrow = 2, ncol = nrow(md), dimnames = list(c("PDE8B", "ACTA2"), NULL))
  out <- pah_pseudobulk_counts(x, md, min_cells = 2)
  expect_equal(dim(out$counts), c(2L, 6L))
  expect_true(all(out$counts == 2L))
  expect_equal(rownames(out$counts), c("PDE8B", "ACTA2"))

  sparse <- Matrix::Matrix(x, sparse = TRUE)
  sparse_out <- pah_pseudobulk_counts(sparse, md, min_cells = 2)
  expect_equal(sparse_out$counts, out$counts)
})
