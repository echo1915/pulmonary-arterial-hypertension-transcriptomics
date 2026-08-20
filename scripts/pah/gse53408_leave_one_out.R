options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(limma))

project <- "project-001-pulmonary-arterial-hypertension-transcriptomics"
data_root <- Sys.getenv("PAH_DATA_ROOT", file.path("G:/workdata/projects", project))
input_dir <- file.path(data_root, "interim", "gse53408_rma")
out_dir <- file.path(data_root, "outputs", "gse53408-sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_df <- read.csv(gzfile(file.path(input_dir, "gse53408_rma_gene_expression.csv.gz")), check.names = FALSE)
symbols <- expr_df[[1]]
expr <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr) <- symbols
gsm <- colnames(expr)
gsm_num <- as.integer(sub("GSM", "", gsm))
group <- ifelse(gsm_num <= 1290999L, "PAH", "Control")
candidate <- c("PDE8B", "PIEZO2", "SLC16A12")

hedges_g <- function(x_case, x_control) {
  n1 <- length(x_case); n0 <- length(x_control)
  sp <- sqrt(((n1 - 1) * var(x_case) + (n0 - 1) * var(x_control)) / (n1 + n0 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  d <- (mean(x_case) - mean(x_control)) / sp
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  correction * d
}

fit_subset <- function(keep, omitted_gsm, omitted_group) {
  sub_group <- factor(group[keep], levels = c("Control", "PAH"))
  design <- model.matrix(~ 0 + sub_group)
  colnames(design) <- c("Control", "PAH")
  fit <- lmFit(expr[, keep, drop = FALSE], design)
  fit <- contrasts.fit(fit, makeContrasts(PAH_vs_Control = PAH - Control, levels = design))
  fit <- eBayes(fit, robust = TRUE)
  de <- topTable(fit, coef = "PAH_vs_Control", number = Inf, sort.by = "none")
  de$symbol <- rownames(de)
  rows <- de[match(candidate, de$symbol), , drop = FALSE]
  rows$hedges_g <- vapply(candidate, function(gene) {
    z <- expr[gene, keep]
    hedges_g(z[sub_group == "PAH"], z[sub_group == "Control"])
  }, numeric(1))
  data.frame(
    omitted_gsm = omitted_gsm,
    omitted_group = omitted_group,
    n_PAH = sum(sub_group == "PAH"),
    n_Control = sum(sub_group == "Control"),
    symbol = candidate,
    logFC = rows$logFC,
    t = rows$t,
    P.Value = rows$P.Value,
    adj.P.Val = rows$adj.P.Val,
    hedges_g = rows$hedges_g,
    stringsAsFactors = FALSE
  )
}

full <- fit_subset(rep(TRUE, length(gsm)), "None", "None")
loo <- do.call(rbind, lapply(seq_along(gsm), function(i) {
  keep <- seq_along(gsm) != i
  fit_subset(keep, gsm[i], group[i])
}))
full_g <- setNames(full$hedges_g, full$symbol)
loo$full_hedges_g <- full_g[loo$symbol]
loo$delta_g_from_full <- loo$hedges_g - loo$full_hedges_g
loo$direction_concordant <- sign(loo$hedges_g) == sign(loo$full_hedges_g)
loo$fdr_lt_0.05 <- loo$adj.P.Val < 0.05

summary_rows <- do.call(rbind, lapply(candidate, function(gene) {
  z <- loo[loo$symbol == gene, ]
  data.frame(
    symbol = gene,
    full_logFC = full$logFC[full$symbol == gene],
    full_hedges_g = full$hedges_g[full$symbol == gene],
    full_fdr = full$adj.P.Val[full$symbol == gene],
    loo_min_g = min(z$hedges_g),
    loo_max_g = max(z$hedges_g),
    loo_median_g = median(z$hedges_g),
    direction_concordance_n = sum(z$direction_concordant),
    direction_concordance_fraction = mean(z$direction_concordant),
    fdr_lt_0.05_n = sum(z$fdr_lt_0.05),
    fdr_lt_0.05_fraction = mean(z$fdr_lt_0.05),
    max_abs_delta_g = max(abs(z$delta_g_from_full)),
    most_influential_sample = z$omitted_gsm[which.max(abs(z$delta_g_from_full))],
    most_influential_group = z$omitted_group[which.max(abs(z$delta_g_from_full))],
    worst_fdr = max(z$adj.P.Val),
    stringsAsFactors = FALSE
  )
}))

write.csv(full, file.path(out_dir, "gse53408_full_model_candidates.csv"), row.names = FALSE)
write.csv(loo, file.path(out_dir, "gse53408_candidate_leave_one_sample_out.csv"), row.names = FALSE)
write.csv(summary_rows, file.path(out_dir, "gse53408_candidate_loo_summary.csv"), row.names = FALSE)
print(summary_rows)
