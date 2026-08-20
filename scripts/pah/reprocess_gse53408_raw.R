options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(AnnotationDbi)
  library(hugene10sttranscriptcluster.db)
})

project <- "project-001-pulmonary-arterial-hypertension-transcriptomics"
data_root <- Sys.getenv("PAH_DATA_ROOT", file.path("G:/workdata/projects", project))
cel_dir <- file.path(data_root, "raw", "gse53408", "cel_gz")
out_dir <- file.path(data_root, "interim", "gse53408_rma")
audit_dir <- file.path(data_root, "outputs", "current-results", "pah_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(cel_dir, pattern = "[.]CEL[.]gz$", full.names = TRUE, ignore.case = TRUE))
stopifnot(length(files) == 23L)
gsm <- sub("_.*$", "", basename(files))
gsm_num <- as.integer(sub("GSM", "", gsm))
group <- ifelse(gsm_num <= 1290999L, "PAH", "Control")
stopifnot(sum(group == "PAH") == 12L, sum(group == "Control") == 11L)

message("Reading 23 CEL files")
raw <- read.celfiles(files)
sampleNames(raw) <- gsm
pData(raw)$gsm <- gsm
pData(raw)$group <- factor(group, levels = c("Control", "PAH"))

message("Running RMA at transcript-cluster core level")
eset <- rma(raw, target = "core")
expr_probe <- exprs(eset)
probe_id <- rownames(expr_probe)
symbol <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = probe_id,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)
keep <- !is.na(symbol) & nzchar(symbol)
expr_gene <- avereps(expr_probe[keep, , drop = FALSE], ID = unname(symbol[keep]))

design <- model.matrix(~ 0 + pData(eset)$group)
colnames(design) <- c("Control", "PAH")
fit <- lmFit(expr_gene, design)
fit <- contrasts.fit(fit, makeContrasts(PAH_vs_Control = PAH - Control, levels = design))
fit <- eBayes(fit, robust = TRUE)
de <- topTable(fit, coef = "PAH_vs_Control", number = Inf, sort.by = "P")
de$symbol <- rownames(de)

hedges_g <- function(x_case, x_control) {
  n1 <- length(x_case); n0 <- length(x_control)
  sp <- sqrt(((n1 - 1) * var(x_case) + (n0 - 1) * var(x_control)) / (n1 + n0 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  d <- (mean(x_case) - mean(x_control)) / sp
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  correction * d
}
g <- apply(expr_gene, 1, function(z) hedges_g(z[group == "PAH"], z[group == "Control"]))
de$hedges_g <- g[de$symbol]
de <- de[, c("symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "hedges_g")]

candidate <- c("PDE8B", "PIEZO2", "SLC16A12")
candidate_de <- de[match(candidate, de$symbol), , drop = FALSE]
candidate_de$present <- !is.na(candidate_de$symbol)

sample_meta <- data.frame(gsm = gsm, group = group, file = basename(files))
sample_meta$median_raw_log2 <- apply(log2(exprs(raw) + 1), 2, median, na.rm = TRUE)
sample_meta$median_rma <- apply(expr_probe, 2, median, na.rm = TRUE)
probe_medians <- apply(expr_probe, 1, median, na.rm = TRUE)
rle <- sweep(expr_probe, 1, probe_medians, '-')
sample_meta$rle_median <- apply(rle, 2, median, na.rm = TRUE)
sample_meta$rle_iqr <- apply(rle, 2, IQR, na.rm = TRUE)

pca <- prcomp(t(expr_gene), scale. = FALSE)
pca_out <- data.frame(
  gsm = gsm,
  group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC1_variance = summary(pca)$importance[2, 1],
  PC2_variance = summary(pca)$importance[2, 2]
)

write.csv(sample_meta, file.path(out_dir, "gse53408_sample_qc.csv"), row.names = FALSE)
write.csv(pca_out, file.path(out_dir, "gse53408_pca_coordinates.csv"), row.names = FALSE)
write.csv(de, gzfile(file.path(out_dir, "gse53408_full_differential_expression.csv.gz")), row.names = FALSE)
write.csv(candidate_de, file.path(audit_dir, "gse53408_raw_cel_candidate_validation.csv"), row.names = FALSE)
write.csv(data.frame(symbol = rownames(expr_gene), expr_gene), gzfile(file.path(out_dir, "gse53408_rma_gene_expression.csv.gz")), row.names = FALSE)

png(file.path(out_dir, "gse53408_rma_qc_draft.png"), width = 1600, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
boxplot(expr_probe, names = gsm, las = 2, cex.axis = 0.6, ylab = "RMA log2 expression", main = "RMA distributions")
cols <- ifelse(group == "PAH", "#C94C4C", "#2F6B9A")
plot(pca$x[, 1], pca$x[, 2], pch = 19, col = cols,
     xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]),
     main = "PCA of gene-level RMA expression")
text(pca$x[, 1], pca$x[, 2], labels = gsm, pos = 3, cex = 0.55)
legend("topright", legend = c("Control", "PAH"), col = c("#2F6B9A", "#C94C4C"), pch = 19, bty = "n")
dev.off()

summary_out <- list(
  accession = "GSE53408",
  platform = "GPL6244 HuGene-1_0-st",
  samples = ncol(expr_gene),
  PAH = sum(group == "PAH"),
  controls = sum(group == "Control"),
  transcript_clusters = nrow(expr_probe),
  mapped_gene_symbols = nrow(expr_gene),
  normalization = "oligo::rma(target='core')",
  differential_model = "limma PAH vs Control; robust empirical Bayes",
  candidate_results = candidate_de
)
saveRDS(summary_out, file.path(out_dir, "gse53408_raw_reprocessing_summary.rds"))
writeLines(capture.output(str(summary_out)), file.path(out_dir, "gse53408_raw_reprocessing_summary.txt"))
print(candidate_de)


