args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript script.R <input_dir> <output_dir> <library_dir>")
input_dir <- args[[1]]
output_dir <- args[[2]]
library_dir <- args[[3]]
.libPaths(c(library_dir, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldKnk)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260813)
counts <- readMM(file.path(input_dir, "counts.mtx"))
genes <- readLines(file.path(input_dir, "genes.tsv"))
cells <- read.csv(file.path(input_dir, "cells.csv"), check.names = FALSE)
colnames(counts) <- cells[[1]]

# Collapse duplicate symbols, retain genes detected in >= 1% of cells, then select
# high-dispersion genes. PDE8B is explicitly retained for a feasibility test.
groups <- split(seq_along(genes), genes)
collapsed <- do.call(rbind, lapply(groups, function(ii) Matrix::colSums(counts[ii, , drop = FALSE])))
rownames(collapsed) <- names(groups)
detected <- Matrix::rowSums(collapsed > 0)
keep <- detected >= max(5, ceiling(ncol(collapsed) * 0.01))
keep["PDE8B"] <- TRUE
x <- collapsed[keep, , drop = FALSE]
mu <- Matrix::rowMeans(x)
v <- apply(as.matrix(x), 1, var)
disp <- (v - mu) / pmax(mu^2, 1e-8)
selected <- names(sort(disp, decreasing = TRUE))[seq_len(min(119, length(disp)))]
selected <- unique(c("PDE8B", selected))
x <- x[selected, , drop = FALSE]

audit <- data.frame(
  cells = ncol(x), genes = nrow(x),
  pde8b_nonzero_cells = sum(x["PDE8B", ] > 0),
  pde8b_total_counts = sum(x["PDE8B", ]),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(output_dir, "smoke_input_audit.csv"), row.names = FALSE)
writeLines(rownames(x), file.path(output_dir, "smoke_selected_genes.txt"))

result <- scTenifoldKnk(
  countMatrix = as.matrix(x), gKO = "PDE8B",
  qc_mtThreshold = 1, qc_minLSize = 0,
  nc_nNet = 1, nc_nCells = min(100, ncol(x)), nc_nComp = 3,
  nc_q = 0.9, td_K = 1, td_maxIter = 50, ma_nDim = 2
)
saveRDS(result, file.path(output_dir, "smoke_result.rds"))
write.csv(result$diffRegulation, file.path(output_dir, "smoke_diff_regulation.csv"), row.names = TRUE)
wt <- result$tensorNetworks$WT
ko <- result$tensorNetworks$KO
pde_out <- sum(abs(wt["PDE8B", ]))
network_audit <- data.frame(
  wt_genes = nrow(wt), pde8b_wt_outdegree_weight = pde_out,
  pde8b_ko_outdegree_weight = sum(abs(ko["PDE8B", ])),
  significant_dr_genes = sum(result$diffRegulation$FDR < 0.05, na.rm = TRUE)
)
write.csv(network_audit, file.path(output_dir, "smoke_network_audit.csv"), row.names = FALSE)
print(audit)
print(network_audit)
