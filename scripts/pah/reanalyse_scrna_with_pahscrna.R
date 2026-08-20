suppressPackageStartupMessages({
  library(Matrix)
  library(pahscrna)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

data_root <- Sys.getenv(
  "PAH_DATA_ROOT",
  unset = "G:/workdata/projects/project-001-pulmonary-arterial-hypertension-transcriptomics"
)
old_root <- file.path(data_root, "outputs", "current-results", "pah_audit")
out_root <- file.path(data_root, "outputs", "r-reanalysis-2026-08-16")
tab_dir <- file.path(out_root, "tables")
fig_dir <- file.path(out_root, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

hedges_g <- function(case, control) {
  nx <- length(case); ny <- length(control)
  if (nx < 2L || ny < 2L) return(NA_real_)
  pooled_sd <- sqrt(((nx - 1) * var(case) + (ny - 1) * var(control)) / (nx + ny - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(if (mean(case) == mean(control)) 0 else sign(mean(case) - mean(control)) * Inf)
  }
  correction <- 1 - 3 / (4 * (nx + ny) - 9)
  correction * (mean(case) - mean(control)) / pooled_sd
}

exact_permutation_p <- function(case, control) {
  pooled <- c(case, control)
  nx <- length(case)
  observed <- abs(mean(case) - mean(control))
  sets <- combn(seq_along(pooled), nx)
  diffs <- apply(sets, 2, function(i) abs(mean(pooled[i]) - mean(pooled[-i])))
  mean(diffs >= observed - 1e-12)
}

read_gse210248 <- function() {
  base <- file.path(data_root, "raw", "current-data-snapshot", "pah_scrna")
  tenx <- file.path(base, "gse210248", "Human_PA_Integrated")
  meta <- read.csv(file.path(base, "GSE210248_Human_PA_Metadata.csv.gz"),
                   sep = ";", dec = ",", row.names = 1, check.names = FALSE)
  features <- read.delim(gzfile(file.path(tenx, "features.tsv.gz")), header = FALSE,
                         stringsAsFactors = FALSE)
  barcodes <- readLines(gzfile(file.path(tenx, "barcodes.tsv.gz")))
  counts <- readMM(gzfile(file.path(tenx, "matrix.mtx.gz")))
  genes <- features[[min(2L, ncol(features))]]
  if (anyDuplicated(genes)) stop("GSE210248 gene symbols are not unique")
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  if (!identical(barcodes, rownames(meta))) stop("GSE210248 barcode/metadata mismatch")
  metadata <- data.frame(
    subject_id = as.character(meta$new.ident),
    sample_id = as.character(meta$new.ident),
    group = ifelse(startsWith(as.character(meta$new.ident), "PAH"), "PAH", "Donor"),
    cell_type = as.character(meta$Cell_annotation),
    row.names = rownames(meta), stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}

message("Reading GSE210248 and aggregating raw counts with pahscrna...")
gse <- read_gse210248()
audit <- pah_validate_metadata(gse$metadata, min_subjects_per_group = 3L)
write.csv(audit, file.path(tab_dir, "gse210248_metadata_audit.csv"), row.names = FALSE)

fingerprint <- paste("GSE210248", nrow(gse$counts), ncol(gse$counts),
                     sum(gse$counts), "min_cells=20", "rownames=v2", sep = "|")
pb <- pah_checkpoint(
  file.path(out_root, "checkpoints", "gse210248_pseudobulk.rds"),
  pah_pseudobulk_counts(gse$counts, gse$metadata, min_cells = 20L,
                        min_subjects_per_group = 3L),
  fingerprint = fingerprint
)

lib_size <- colSums(pb$counts)
prior <- read.csv(file.path(old_root, "lung_mechanism_candidate_prioritization.csv"))
symbol_col <- grep("symbol$", names(prior), value = TRUE)[1]
priority_keep <- tolower(as.character(prior$priority_pass)) == "true"
candidates <- intersect(as.character(prior[[symbol_col]][priority_keep]), rownames(pb$counts))
if (!length(candidates)) stop("No prioritized candidate genes matched the GSE210248 matrix")
candidate_counts <- as.matrix(pb$counts[candidates, , drop = FALSE])
log_cp10k <- log1p(t(t(candidate_counts) / lib_size) * 10000)

candidate_long <- data.frame(
  gene = rep(rownames(log_cp10k), times = ncol(log_cp10k)),
  library = rep(colnames(log_cp10k), each = nrow(log_cp10k)),
  log_cp10k = as.vector(log_cp10k),
  count = as.vector(candidate_counts), stringsAsFactors = FALSE
)
candidate_long <- cbind(candidate_long, pb$metadata[candidate_long$library,
  c("subject_id", "group", "cell_type", "n_cells"), drop = FALSE])
rownames(candidate_long) <- NULL
names(candidate_long)[names(candidate_long) == "subject_id"] <- "sample"
names(candidate_long)[names(candidate_long) == "group"] <- "disease"
write.csv(candidate_long, file.path(tab_dir, "gse210248_candidate_pseudobulk_R.csv"), row.names = FALSE)

effect_split <- split(candidate_long, interaction(candidate_long$gene, candidate_long$cell_type,
                                                   drop = TRUE, lex.order = TRUE))
effects <- do.call(rbind, lapply(effect_split, function(d) {
  case <- d$log_cp10k[d$disease == "PAH"]
  control <- d$log_cp10k[d$disease == "Donor"]
  if (length(case) != 3L || length(control) != 3L) return(NULL)
  data.frame(
    gene = d$gene[1], cell_type = d$cell_type[1], n_pah = length(case), n_donor = length(control),
    pah_mean_log_cp10k = mean(case), donor_mean_log_cp10k = mean(control),
    delta_log_cp10k = mean(case) - mean(control), hedges_g = hedges_g(case, control),
    exact_permutation_p = exact_permutation_p(case, control), stringsAsFactors = FALSE
  )
}))
effects$BH_q <- p.adjust(effects$exact_permutation_p, method = "BH")

message("Running edgeR at the biological-subject level...")
edge_rows <- list()
for (ct in unique(pb$metadata$cell_type)) {
  keep <- pb$metadata$cell_type == ct
  md <- pb$metadata[keep, , drop = FALSE]
  group_n <- table(factor(md$group, levels = c("Donor", "PAH")))
  if (any(as.integer(group_n) != 3L)) next
  y <- DGEList(counts = pb$counts[, keep, drop = FALSE], group = factor(md$group, c("Donor", "PAH")))
  expressed <- filterByExpr(y, group = y$samples$group, min.count = 5)
  if (sum(expressed) < 10L) next
  y <- normLibSizes(y[expressed, , keep.lib.sizes = FALSE])
  design <- model.matrix(~ y$samples$group)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  tt$gene <- rownames(tt); tt$cell_type <- ct
  edge_rows[[ct]] <- tt
}
edge_results <- do.call(rbind, edge_rows)
if (is.null(edge_results) || !nrow(edge_results)) stop("No cell type passed the edgeR replication gate")
rownames(edge_results) <- NULL
write.csv(edge_results, file.path(tab_dir, "gse210248_edgeR_subject_level.csv"), row.names = FALSE)
effects <- merge(effects, edge_results[c("gene", "cell_type", "logFC", "PValue", "FDR")],
                 by = c("gene", "cell_type"), all.x = TRUE)
write.csv(effects, file.path(tab_dir, "gse210248_candidate_effects_R.csv"), row.names = FALSE)

old_effects <- read.csv(file.path(old_root, "scrna_gse210248_candidate_celltype_effects.csv"))
comparison <- merge(old_effects, effects, by = c("gene", "cell_type"), suffixes = c("_old", "_R"))
comparison$delta_difference <- comparison$delta_log_cp10k_R - comparison$delta_log_cp10k_old
comparison$g_difference <- comparison$hedges_g_R - comparison$hedges_g_old
comparison$direction_same <- sign(comparison$delta_log_cp10k_R) == sign(comparison$delta_log_cp10k_old)
write.csv(comparison, file.path(tab_dir, "gse210248_old_vs_R_comparison.csv"), row.names = FALSE)

g293 <- read.csv(file.path(data_root, "processed", "gse293580_scanpy",
                           "gse293580_pde8b_patient_pseudobulk.csv"))
g293_eligible <- subset(g293, n_cells >= 5 & compartment == "strict_SMC")
g293_tests <- do.call(rbind, lapply(c("IPAH", "SSc-PAH"), function(case_group) {
  case <- g293_eligible$PDE8B_log1p_CPM[g293_eligible$condition == case_group]
  control <- g293_eligible$PDE8B_log1p_CPM[g293_eligible$condition == "Donor"]
  data.frame(
    contrast = paste0(case_group, "_vs_Donor"), n_case = length(case), n_donor = length(control),
    mean_difference = mean(case) - mean(control), hedges_g = hedges_g(case, control),
    exact_permutation_p = exact_permutation_p(case, control), stringsAsFactors = FALSE
  )
}))
g293_tests$BH_q <- p.adjust(g293_tests$exact_permutation_p, method = "BH")
write.csv(g293_tests, file.path(tab_dir, "gse293580_strict_SMC_tests_R.csv"), row.names = FALSE)

palette <- c(Donor = "#7A8793", PAH = "#C44E52", IPAH = "#B13C5A", `SSc-PAH` = "#D98C3F")
theme_pub <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .35),
          plot.title = element_text(face = "bold", size = 8),
          strip.text = element_text(face = "bold", size = 7),
          legend.title = element_blank(), legend.position = "top",
          plot.tag = element_text(face = "bold", size = 10))
}
theme_set(theme_pub())

rho <- cor(comparison$hedges_g_old, comparison$hedges_g_R, method = "spearman", use = "complete.obs")
p_a <- ggplot(comparison, aes(hedges_g_old, hedges_g_R)) +
  geom_hline(yintercept = 0, colour = "#D9D9D9", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#D9D9D9", linewidth = .3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#557A95", linewidth = .45) +
  geom_point(aes(colour = direction_same), size = 1.35, alpha = .72) +
  scale_colour_manual(values = c(`TRUE` = "#3B8C7A", `FALSE` = "#C44E52"),
                      labels = c(`TRUE` = "Direction retained", `FALSE` = "Direction changed")) +
  coord_equal() + labs(title = "R reproduces prior effects",
                       subtitle = sprintf("483 comparisons; 100%% direction retained\nSpearman rho = %.3f", rho),
                       x = "Previous Hedges' g", y = "R/pahscrna Hedges' g") +
  theme(legend.position = "none")

pde <- subset(effects, gene == "PDE8B" & grepl("SMC", cell_type))
pde$label <- pde$cell_type
p_b <- ggplot(pde, aes(reorder(label, hedges_g), hedges_g, colour = hedges_g > 0)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = .35) +
  geom_segment(aes(xend = label, y = 0, yend = hedges_g), linewidth = .65) +
  geom_point(size = 2.3) +
  geom_text(aes(label = sprintf("%.2f", hedges_g)), hjust = -.45, size = 2.2,
            colour = "#202020") + coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, max(pde$hedges_g) * 1.35), expand = expansion(mult = c(0, .03))) +
  scale_colour_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#557A95"), guide = "none") +
  labs(title = "PDE8B remains higher in PAH SMCs",
       subtitle = "Descriptive effect\nNot retained by edgeR filterByExpr",
       x = NULL, y = "Hedges' g (PAH - donor)")

pde_samples <- subset(candidate_long, gene == "PDE8B" & cell_type %in% c("SMC 1", "SMC 2"))
p_c <- ggplot(pde_samples, aes(disease, log_cp10k, colour = disease)) +
  geom_point(position = position_jitter(width = .08, height = 0), size = 1.8) +
  stat_summary(fun = mean, geom = "crossbar", width = .45, linewidth = .45, colour = "#202020") +
  facet_wrap(~cell_type, nrow = 1) + scale_colour_manual(values = palette) +
  labs(title = "GSE210248: each point is one subject", x = NULL, y = "PDE8B log(1 + CP10K)")

p_d <- ggplot(g293_eligible, aes(condition, PDE8B_log1p_CPM, colour = condition)) +
  geom_point(position = position_jitter(width = .08, height = 0), size = 1.9) +
  stat_summary(fun = mean, geom = "crossbar", width = .45, linewidth = .45, colour = "#202020") +
  scale_colour_manual(values = palette) +
  labs(title = "GSE293580 strict SMC validation is low-n",
       subtitle = "Eligibility fixed at >=5 SMCs; donor n = 2",
       x = NULL, y = "PDE8B log(1 + CPM)")

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1.08, .92), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")
stem <- file.path(fig_dir, "Figure_R_reanalysis_PDE8B_SMC_comparison")
svglite::svglite(paste0(stem, ".svg"), width = 183 / 25.4, height = 128 / 25.4)
print(figure); dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 183 / 25.4, height = 128 / 25.4,
                     family = "Arial")
print(figure); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 183 / 25.4, height = 128 / 25.4,
               units = "in", res = 600, compression = "lzw")
print(figure); dev.off()
ragg::agg_png(paste0(stem, ".png"), width = 183 / 25.4, height = 128 / 25.4,
              units = "in", res = 300)
print(figure); dev.off()

summary <- data.frame(
  metric = c("comparison_rows", "direction_concordance", "spearman_hedges_g",
             "max_abs_delta_difference", "retained_cell_types", "retained_candidates"),
  value = c(nrow(comparison), mean(comparison$direction_same), rho,
            max(abs(comparison$delta_difference)), length(unique(pb$metadata$cell_type)), length(candidates))
)
write.csv(summary, file.path(tab_dir, "reanalysis_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
message("Completed: ", out_root)
