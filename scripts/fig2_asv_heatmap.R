#!/usr/bin/env Rscript
# fig2_asv_heatmap.R
#
# ASV-level family heatmap, alpha/beta diversity, and DESeq2 differential
# abundance analysis. Paths come from a config .ini file; sample order,
# treatment groups, comparisons, and colors come from a per-study R file
# referenced inside that config (see study_params_template.R).
#
# Usage: Rscript fig2_asv_heatmap.R <path_to_config_file>

library(dada2)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(phyloseq)
library(vegan)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(scales)
library(forcats)
library(rstatix)
library(DESeq2)
library(apeglm)

source("utils.R")  # also provides read_config()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript fig2_asv_heatmap.R <path_to_config_file>")
}

config <- read_config(args[1], required_fields = c(
  "data_output_dir", "figure_output_dir", "metadata_path", "study_params"
))

dada2_data_out <- config$data_output_dir
dada2_figure_out <- config$figure_output_dir
dir.create(dada2_figure_out, recursive = TRUE, showWarnings = FALSE)

# Per-study parameters: sample_order, treatment_levels, comparisons,
# samples_to_remove, family_colors
source(config$study_params)

meta_data <- read.table(config$metadata_path, header = TRUE,
                        sep = ",", stringsAsFactors = FALSE)
rownames(meta_data) <- meta_data$sample_id

seqtab.nochim <- readRDS(file.path(dada2_data_out, "seqtab_final.rds"))
seqtab.filt <- seqtab.nochim

# Confirm seqtab and metadata have the same samples
setequal(rownames(seqtab.filt), rownames(meta_data))
meta_data <- meta_data[rownames(seqtab.filt), ]

# assign ASVs to taxa (family or genus)
taxa <- read.table(
  file.path(dada2_data_out, "taxa.csv"),
  sep = " ",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE
)

#
# Filtering ASVs
#

taxa_filtered <- taxa[!taxa$label %in% c("NA_NA_NA", NA), ]
seqtab.filt2 <- seqtab.filt[, colnames(seqtab.filt) %in% rownames(taxa_filtered)]

# get rid of all-zero ASVs
seqtab.filt2 <- seqtab.filt2[, colSums(seqtab.filt2) > 0]
taxa_filtered <- taxa_filtered[colnames(seqtab.filt2), , drop = FALSE]

p <- plot_qc_filtering_subjects(seqtab.filt2, meta_data, "sample_id",
                                relabund_thresh_range = seq(0, 0.01, length.out = 50),
                                min_subjects = c(1, 2, 3))
ggsave(file.path(dada2_figure_out, "dada2_filtering_elbow.pdf"), plot = p, width = 6, height = 4)

seqtab.filt3 <- rel_abund_filter_subject(seqtab.filt2, meta_data, "sample_id",
                                        relabund_thresh = 0.001,
                                        n_subjects_thresh = 2)
print(paste("num samples:", nrow(seqtab.filt3),
            "num ASVs:", ncol(seqtab.filt3)))

taxa_combined <- taxa_filtered[rownames(taxa_filtered) %in% colnames(seqtab.filt3), , drop = FALSE]

asv_mapping_res <- assign_asv_ids(seqtab.filt3, taxa_combined)
seqtab_combined <- asv_mapping_res$seqtab
asv_mapping <- asv_mapping_res$mapping
rownames(asv_mapping) <- asv_mapping$ASV
taxa_combined$ASV_ID <- asv_mapping$ASV_ID[match(rownames(taxa_combined), asv_mapping$Sequence)]
rownames(taxa_combined) <- taxa_combined$ASV_ID
write.csv(asv_mapping,
          file.path(dada2_data_out, "seqtab_combined_asv_mapping.csv"),
          row.names = FALSE)

#
# Family-level heatmap
# (sample_order and family_colors come from the study params file)
#

setequal(rownames(meta_data), sample_order)

result <- plot_family_heatmap_v2(
  seqtab = seqtab_combined,
  taxa = taxa_combined,
  top_n = 15,
  sample_order = sample_order,
  family_colors = family_colors,
  min_rel_abund = 0.01,
  out_file = file.path(dada2_figure_out, "fig2_family_top15_stacked_bars_v2.pdf")
)

### use seqtab.nochim and taxa (no filtering, include singletons)
colnames(seqtab_combined) <- names(colnames(seqtab_combined))
taxa_filtered <- taxa[colnames(seqtab_combined), , drop = FALSE]
all(rownames(taxa_filtered) == colnames(seqtab_combined))  # sanity check

meta_data <- meta_data[match(rownames(seqtab_combined), meta_data$sample_id), ]
taxa_mat <- as.matrix(taxa_filtered)
ps <- phyloseq(
  otu_table(seqtab_combined, taxa_are_rows = FALSE),
  sample_data(meta_data),
  tax_table(taxa_mat)
)

# remove technical replicates / excluded samples (from study params)
ps <- subset_samples(ps, !(sample_id %in% samples_to_remove))

# ---------- Alpha Diversity Analysis ----------

run_alpha_diversity <- function(ps, groups, out_name) {
  num_species <- ntaxa(ps) # for normalized Shannon entropy
  alpha_div <- estimate_richness(ps, measures = c("Shannon")) / log(num_species)
  alpha_div$sample_id <- rownames(alpha_div)
  alpha_div_merged <- merge(alpha_div, meta_data, by.x = "sample_id", by.y = "sample_id", all.x = TRUE)
  alpha_div_merged$treatment <- factor(alpha_div_merged$treatment, levels = groups)

  p <- ggplot(alpha_div_merged, aes(x = treatment, y = Shannon)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_continuous(limits = range(alpha_div_merged$Shannon) + c(-0.1, 0.1)) +
    labs(x = "Sample Type", y = "Normalized Shannon Diversity Index") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(dada2_figure_out, sprintf("alpha_diversity_Shannon_boxplot_vs_treatment%s.pdf", out_name)),
          plot = p, width = 6, height = 6, dpi = 300)

  alpha_div_wilcox <- alpha_div_merged %>%
    pairwise_wilcox_test(Shannon ~ treatment, p.adjust.method = "BH") %>%
    add_significance()
  write.csv(alpha_div_wilcox, file.path(dada2_data_out, sprintf("alpha_diversity_Shannon_wilcox_test_results%s.csv", out_name)))

  alpha_div_IQR <- alpha_div_merged %>% group_by(treatment) %>%
    summarise(median = median(Shannon), IQR = IQR(Shannon))
  write.csv(alpha_div_IQR, file.path(dada2_data_out, sprintf("alpha_diversity_Shannon_IQR_results%s.csv", out_name)))

  write.csv(alpha_div_merged, file.path(dada2_data_out, sprintf("alpha_diversity_Shannon_per_sample%s.csv", out_name)), row.names = FALSE)
}

# treatment_levels comes from the study params file
run_alpha_diversity(ps, treatment_levels, "")

# ---------- Beta Diversity Analysis ----------

bc_dist <- phyloseq::distance(ps, method = "bray")
meta <- data.frame(sample_data(ps))

set.seed(101)
adonis_res <- adonis2(bc_dist ~ treatment, permutations = 9999, data = meta, by = "terms")
print(adonis_res)
write.csv(adonis_res, file.path(dada2_data_out, "beta_diversity_bray_curtis_adonis_results.csv"))

disp <- betadisper(bc_dist, meta$treatment)
perm_res <- permutest(disp)
write.csv(as.data.frame(perm_res$tab), file.path(dada2_data_out, "beta_diversity_bray_curtis_betadisper_results.csv"))

p1 <- plot_beta_diversity(ps,
                         color_var = "treatment",
                         label_var = "sample_id",
                         title = "PCoA of Bray-Curtis Distances by Treatment")
ggsave(file.path(dada2_figure_out, "beta_diversity_pcoa_treatment.pdf"), plot = p1)

#
# ---------- Differential abundance analysis ----------
#

seqtab_combined <- asv_mapping_res$seqtab
names(colnames(seqtab_combined)) <- NULL
taxa_combined <- taxa_combined[colnames(seqtab_combined), , drop = FALSE]
all(rownames(taxa_combined) == colnames(seqtab_combined))  # sanity check
taxa_mat <- as.matrix(taxa_combined)
ps <- phyloseq(
  otu_table(seqtab_combined, taxa_are_rows = FALSE),
  sample_data(meta_data),
  tax_table(taxa_mat)
)
ps <- subset_samples(ps, !(sample_id %in% samples_to_remove))

run_deseq2_phyloseq <- function(ps, group_var, group_A, group_B,
                               min_count = 10,
                               shrink = TRUE,
                               plot_disp = FALSE,
                               fit_type = "parametric") {
  require(DESeq2)
  require(apeglm)
  require(phyloseq)

  sample_df <- data.frame(sample_data(ps))

  keep_samples <- rownames(sample_df)[sample_df[[group_var]] %in% c(group_A, group_B)]
  ps_sub <- prune_samples(keep_samples, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  if (!is.null(min_count)) {
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > min_count, ps_sub)
  }

  design_formula <- as.formula(paste("~", group_var))
  dds <- phyloseq_to_deseq2(ps_sub, design_formula)

  colData(dds)[[group_var]] <- relevel(
    factor(colData(dds)[[group_var]]),
    ref = group_A
  )

  if (fit_type %in% c("parametric", "local")) {
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds, fitType = fit_type)
    dds <- nbinomWaldTest(dds)
  } else {
    stop("Invalid fit_type. Use 'parametric' or 'local'.")
  }

  if (plot_disp) {
    pdf(file.path(dada2_figure_out, paste0("deseq2_dispersion_plot_", group_B, "_vs_", group_A, ".pdf")))
    plotDispEsts(dds)
    dev.off()
  }

  res <- results(dds, contrast = c(group_var, group_B, group_A))

  if (shrink) {
    coef_name <- paste0(group_var, "_", group_B, "_vs_", group_A)
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    res_shrunk <- NULL
  }

  return(list(dds = dds, results = res, results_shrunk = res_shrunk))
}

make_df <- function(res, label) {
  df <- as.data.frame(res)
  df$ASV <- rownames(df)
  df$comparison <- label
  return(df)
}

# Run every comparison defined in the study params file
df_list <- list()
for (cmp in comparisons) {
  res_list <- run_deseq2_phyloseq(
    ps,
    group_var = "treatment",
    group_A = cmp$group_A,
    group_B = cmp$group_B,
    min_count = 10,
    shrink = TRUE,
    plot_disp = TRUE,
    fit_type = "local"
  )
  df_list[[cmp$label]] <- make_df(res_list$results_shrunk, cmp$label)
}

df_all <- bind_rows(df_list)

df_all$padj_global <- p.adjust(df_all$pvalue, method = "BH")
df_all$padj_global[is.na(df_all$padj_global)] <- 1

df_all <- df_all %>%
  mutate(
    sig_global = padj_global < 0.05,
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )

# order ASVs by effect size
df_all$ASV <- factor(df_all$ASV,
                     levels = unique(df_all$ASV[order(df_all$log2FoldChange)]))

# join with taxonomic information
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

df_all <- left_join(df_all, tax[, c("ASV", "label")], by = "ASV")
df_all$label <- gsub("_NA", "", df_all$label)
df_all$label <- paste0(df_all$label, "|", df_all$ASV)

write.csv(df_all,
          file.path(dada2_data_out, "deseq2_results_all_comparisons.csv"),
          row.names = FALSE)

# only keep ASVs which are sig in at least one comparison
all_sig_asvs <- unique(df_all$ASV[df_all$sig_global == TRUE])
df_all <- df_all[df_all$ASV %in% all_sig_asvs, ]

write.csv(df_all,
          file.path(dada2_data_out, "deseq2_results_all_comparisons_sig.csv"),
          row.names = FALSE)

# Facet order follows the order comparisons are listed in study_params.R
comparison_levels <- sapply(comparisons, function(x) x$label)
df_all$comparison <- factor(df_all$comparison, levels = comparison_levels)

# +/- direction pattern across all comparisons. Generalized to work with
# any number of comparisons (the original version hardcoded exactly 3).
df_sign <- df_all %>%
  mutate(sign = ifelse(log2FoldChange > 0, "+", "-")) %>%
  select(ASV, comparison, sign) %>%
  distinct()

df_pattern <- df_sign %>%
  pivot_wider(names_from = comparison, values_from = sign)

pattern_cols <- comparison_levels[comparison_levels %in% colnames(df_pattern)]
df_pattern$pattern <- apply(df_pattern[, pattern_cols, drop = FALSE], 1,
                            function(row) {
                              if (any(is.na(row))) return(NA)
                              paste(row, collapse = "")
                            })

df_all2 <- df_all %>%
  left_join(df_pattern[, c("ASV", "pattern")], by = "ASV")

p <- ggplot(df_all2, aes(x = log2FoldChange, y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "black") +
  geom_point(aes(size = -log10(padj_global),
                 color = sig_global)) +
  facet_grid(pattern ~ comparison,
             scales = "free_y",
             space = "free_y") +
  scale_color_manual(
    values = c("FALSE" = "grey70",
               "TRUE"  = "red")
  ) +
  theme_bw() +
  labs(
    x = "log2 fold change",
    y = "ASV",
    size = "-log10(global FDR)",
    color = "Significant"
  )

ggsave(file.path(dada2_figure_out, "deseq2_dotplot_all_comparisons.pdf"),
       plot = p, device = "pdf",
       width = 12, height = 15)