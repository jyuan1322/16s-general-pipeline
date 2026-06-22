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

source("utils.R")

working_dir <- "/data/local/jy1008/project_folder"
meta_data <- read.table(file.path(working_dir, "metadata",
                "pasricha_metadata_from_filelabels.csv"), header = TRUE,
                sep = ",", stringsAsFactors = FALSE)
rownames(meta_data) <- meta_data$sample_id

# data and figure output dirs from dada2
dada2_data_out <- "test_data_output"
dada2_figure_out <- "test_figure_output"

seqtab.nochim <- readRDS(file.path(working_dir,
                    dada2_data_out, "seqtab_final.rds"))
seqtab.filt <- seqtab.nochim

# Confirm seqtab and metadata have the same samples
setequal(rownames(seqtab.filt), rownames(meta_data))
meta_data <- meta_data[rownames(seqtab.filt), ]


# assign ASVs to taxa (family or genus)
taxa <- read.table(
  file.path(working_dir, dada2_data_out, "taxa.csv"),
  sep = " ",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE
)
# colnames(taxa)[1] <- "ASV"

#
# Filtering ASVs
#

# filter out ASVs which are not assigned to any taxa
taxa_filtered <- taxa[!taxa$label %in% c("NA_NA_NA", NA), ]
seqtab.filt2 <- seqtab.filt[, colnames(seqtab.filt) %in%
                    rownames(taxa_filtered)]


# get rid of all-zero ASVs
seqtab.filt2 <- seqtab.filt2[, colSums(seqtab.filt2) > 0]
taxa_filtered <- taxa_filtered[colnames(seqtab.filt2), , drop = FALSE]

# plot elbow plots (from utils.R)
# plot_qc_elbow(seqtab.filt2)
p <- plot_qc_filtering_subjects(seqtab.filt2, meta_data, "sample_id",
                                relabund_thresh_range = seq(0, 0.01, length.out = 50),
                                min_subjects = c(1, 2, 3))
ggsave("dada2_filtering_elbow.pdf", plot = p, width = 6, height = 4)

# TODO: try a filter of 0.1% relative abundance in at least 2 subjects
seqtab.filt3 <- rel_abund_filter_subject(seqtab.filt2, meta_data, "sample_id",
                                        relabund_thresh = 0.001,
                                        n_subjects_thresh = 2)
print(paste("num samples:", nrow(seqtab.filt3),
            "num ASVs:", ncol(seqtab.filt3)))

# combined <- unify_data(seqtab_list, meta_list)
# seqtab_combined <- combined$seqtab
# meta_combined <- combined$meta
taxa_combined <- taxa_filtered[rownames(taxa_filtered) %in% colnames(seqtab.filt3), , drop = FALSE]

asv_mapping_res <- assign_asv_ids(seqtab.filt3, taxa_combined)
seqtab_combined <- asv_mapping_res$seqtab
asv_mapping <- asv_mapping_res$mapping
rownames(asv_mapping) <- asv_mapping$ASV
taxa_combined$ASV_ID <- asv_mapping$ASV_ID[match(rownames(taxa_combined), asv_mapping$Sequence)]
rownames(taxa_combined) <- taxa_combined$ASV_ID
write.csv(asv_mapping,
          file.path(working_dir, dada2_data_out, "seqtab_combined_asv_mapping.csv"),
          row.names = FALSE)

# Get label for a sequence
# asv_mapping_res$get_asv_label("ACGTGTA...")  # replace with ASV sequence

#
# Family-level heatmap
#

family_colors <- c(
  "Lachnospiraceae" = "#0074D9",
  "Lactobacillaceae" = "#FF851B",
  "Muribaculaceae" = "#2ECC40",
  "Prevotellaceae" = "#FF4136",
  "Oscillospiraceae" = "#B10DC9",
  "Ruminococcaceae" = "#FFDC00",
  "Bacteroidaceae" = "#39CCCC",
  "Erysipelotrichaceae" = "#85144B",
  "Peptostreptococcaceae" = "#001F3F",
  "Akkermansiaceae" = "#01FF70",
  "[Eubacterium] coprostanoligenes group" = "#F012BE",
  "Peptococcaceae" = "#3D9970",
  "Eggerthellaceae" = "#7FDBFF",
  "Bifidobacteriaceae" = "#A52A2A",
  "Butyricicoccaceae" = "#FF7F50",
  "Other"   = "#808080"
)

# sample_order <- meta_combined$MHMC.sampleID
sample_order <- c(
  "Control1", "Control2", "Control3", "Control4", "Control5",
  "ASS1", "ASS3", "ASS4",
  "DOMP4", "DOMP5",
  "GLP11", "GLP12", "GLP12filtered", "GLP13", "GLP14",
  "ASSDOMP1", "ASSDOMP2", "ASSDOMP3", "ASSDOMP4", "ASSDOMP5",
  "GLP1ASS1", "GLP1ASS2", "GLP1ASS2filtered", "GLP1ASS3", "GLP1ASS4", "GLP1ASS5",
  "GLPDOMP1", "GLPDOMP2", "GLP1DOMP3", "GLP1DOMP4", "GLP1DOMP5",
  "GLPASSDOMP2", "GLPASSDOMP3"
)
setequal(rownames(meta_data), sample_order)
# result <- plot_family_heatmap(
#   seqtab = seqtab_combined,
#   taxa = taxa_combined,
#   top_n = 15,
#   sample_order = sample_order,
#   family_colors = family_colors,
#   out_file = file.path(working_dir, dada2_figure_out, "fig2_family_top15_stacked_bars.pdf")
# )
result <- plot_family_heatmap_v2(
  seqtab = seqtab_combined,
  taxa = taxa_combined,
  top_n = 15,
  sample_order = sample_order,
  family_colors = family_colors,
  min_rel_abund = 0.01,
  out_file = file.path(working_dir, dada2_figure_out, "fig2_family_top15_stacked_bars_v2.pdf")
)


### use seqtab.nochim and taxa (no filtering, include singletons)
# 5/2/2026: use filtered taxa for alpha/beta diversity instead
# all(rownames(taxa) == colnames(seqtab_combined))

# Keep only taxa present in seqtab_combined
colnames(seqtab_combined) <- names(colnames(seqtab_combined))
taxa_filtered <- taxa[colnames(seqtab_combined), , drop = FALSE]
# sanity check
all(rownames(taxa_filtered) == colnames(seqtab_combined))

meta_data <- meta_data[match(rownames(seqtab_combined), meta_data$sample_id), ]
# taxa_mat <- as.matrix(taxa)
taxa_mat <- as.matrix(taxa_filtered)
ps <- phyloseq(
  otu_table(seqtab_combined, taxa_are_rows = FALSE),
  sample_data(meta_data),
  tax_table(taxa_mat)
)

# remove technical replicates so one subject is represented by one sample
samples_to_remove <- c("GLP12", "GLP1ASS2")
ps <- subset_samples(ps, !(sample_id %in% samples_to_remove))

# ---------- Alpha Diversity Analysis ----------

run_alpha_diversity <- function(ps, groups, out_name) {
  num_species <- ntaxa(ps) # for normalized Shannon entropy
  alpha_div <- estimate_richness(ps, measures = c("Shannon")) / log(num_species)
  alpha_div$sample_id <- rownames(alpha_div)
  alpha_div_merged <- merge(alpha_div, meta_data, by.x = "sample_id", by.y = "sample_id", all.x = TRUE)
  alpha_div_merged$treatment <- factor(alpha_div_merged$treatment,
                                      levels = groups)

  p <- ggplot(alpha_div_merged, aes(x = treatment, y = Shannon)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_continuous(limits = range(alpha_div_merged$Shannon) + c(-0.1, 0.1)) +
    labs(
      x = "Sample Type",
      y = "Normalized Shannon Diversity Index"
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(file.path(working_dir, dada2_figure_out, sprintf("alpha_diversity_Shannon_boxplot_vs_treatment%s.pdf", out_name)),
          plot = p, width = 6, height = 6, dpi = 300)

  alpha_div_wilcox <- alpha_div_merged %>%
    pairwise_wilcox_test(Shannon ~ treatment, p.adjust.method = "BH") %>%
    add_significance()
  write.csv(alpha_div_wilcox, file.path(working_dir, dada2_data_out, sprintf("alpha_diversity_Shannon_wilcox_test_results%s.csv", out_name)))

  alpha_div_IQR <- alpha_div_merged %>% group_by(treatment) %>%
    summarise(median = median(Shannon), IQR = IQR(Shannon))
  write.csv(alpha_div_IQR, file.path(working_dir, dada2_data_out, sprintf("alpha_diversity_Shannon_IQR_results%s.csv", out_name)))

  write.csv(alpha_div_merged, file.path(working_dir, dada2_data_out, sprintf("alpha_diversity_Shannon_per_sample%s.csv", out_name)), row.names = FALSE)
}
# all groups
run_alpha_diversity(ps, c("Control", "ASS", "GLP1_ASS", "ASS_DOMP", "DOMP", "GLP1", "GLP1_DOMP", "GLP1_ASS_DOMP"), "")

# ---------- Beta Diversity Analysis ----------

# Compute Bray-Curtis distance
bc_dist <- phyloseq::distance(ps, method = "bray")

meta <- data.frame(sample_data(ps))

set.seed(101)
adonis_res <- adonis2(bc_dist ~ treatment, permutations = 9999, data = meta, by = "terms")
print(adonis_res)
write.csv(adonis_res, file.path(working_dir, dada2_data_out, "beta_diversity_bray_curtis_adonis_results.csv"))

disp <- betadisper(bc_dist, meta$treatment)
perm_res <- permutest(disp)
perm_df <- as.data.frame(perm_res$tab)
write.csv(perm_df, file.path(working_dir, dada2_data_out, "beta_diversity_bray_curtis_betadisper_results.csv"))

p1 <- plot_beta_diversity(ps,
                         color_var = "treatment",
                         label_var = "sample_id",
                         title = "PCoA of Bray-Curtis Distances by Treatment")
ggsave(file.path(working_dir, dada2_figure_out, "beta_diversity_pcoa_treatment.pdf"), plot=p1)

# Compute Bray-Curtis distance, but only for selected groups
ps_sub <- subset_samples(ps, treatment %in% c("Control", "ASS", "GLP1_ASS", "ASS_DOMP"))
bc_dist <- phyloseq::distance(ps_sub, method = "bray")

meta <- data.frame(sample_data(ps_sub))

set.seed(101)
adonis_res <- adonis2(bc_dist ~ treatment, permutations = 9999, data = meta, by = "terms")
print(adonis_res)
write.csv(adonis_res, file.path(working_dir, dada2_data_out, "beta_diversity_bray_curtis_adonis_results_subset.csv"))

disp <- betadisper(bc_dist, meta$treatment)
perm_res <- permutest(disp)
perm_df <- as.data.frame(perm_res$tab)
write.csv(perm_df, file.path(working_dir, dada2_data_out, "beta_diversity_bray_curtis_betadisper_results_subset.csv"))

p1 <- plot_beta_diversity(ps_sub,
                         color_var = "treatment",
                         label_var = "sample_id",
                         title = "PCoA of Bray-Curtis Distances by Treatment")
ggsave(file.path(working_dir, dada2_figure_out, "beta_diversity_pcoa_treatment_subset.pdf"), plot=p1)


# CLR + Aitchison distance
# library(microbiome)
# 
# # CLR transform
# ps_clr <- microbiome::transform(ps, "clr")  # handles zeros
# 
# # Plot using Euclidean distance = Aitchison
# p <- plot_beta_diversity(
#   ps_clr,
#   color_var = "treatment",
#   dist_method = "euclidean",  # Euclidean on CLR = Aitchison
#   ordinate_method = "PCoA",
#   label_var = "sample_id",
#   output_path = "beta_diversity_clr_aitchison.pdf",
#   title = "Beta Diversity - CLR/Aitchison"
# )



#
# ---------- Differential abundance analysis ----------
#


library(DESeq2)
library(apeglm)

# filtered dataset to phyloseq
seqtab_combined <- asv_mapping_res$seqtab
names(colnames(seqtab_combined)) <- NULL
taxa_combined <- taxa_combined[colnames(seqtab_combined), , drop = FALSE]
# Optional sanity check
all(rownames(taxa_combined) == colnames(seqtab_combined))
taxa_mat <- as.matrix(taxa_combined) # convert to matrix for phyloseq
ps <- phyloseq(
  otu_table(seqtab_combined, taxa_are_rows = FALSE),
  sample_data(meta_data),
  tax_table(taxa_mat)
)
# remove technical replicates so one subject is represented by one sample
samples_to_remove <- c("GLP12", "GLP1ASS2")
ps <- subset_samples(ps, !(sample_id %in% samples_to_remove))

# NOTE: for ASV mapping to work, make sure phyloseq object
# is generated from the same as asv_mapping

# Treatments
# Control, ASS, DOMP, GLP-1
# ASS_DOMP, GLP1_ASS, GLP1_DOMP, GLP1_ASS_DOMP

# Comparisons
# ASS vs Control
# GLP1_ASS vs ASS
# DOMP vs ASS_DOMP

run_deseq2_phyloseq <- function(ps, group_var, group_A, group_B,
                               min_count = 10,
                               shrink = TRUE,
                               plot_disp = FALSE,
                               fit_type = "parametric") {
  require(DESeq2)
  require(apeglm)
  require(phyloseq)
  
  # Extract sample metadata safely
  sample_df <- data.frame(sample_data(ps))
  
  # 1. Subset samples
  keep_samples <- rownames(sample_df)[sample_df[[group_var]] %in% c(group_A, group_B)]
  ps_sub <- prune_samples(keep_samples, ps)
  
  # Drop taxa with zero counts
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  # Optional filtering
  if (!is.null(min_count)) {
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > min_count, ps_sub)
  }
  
  # 2. Convert to DESeq2
  design_formula <- as.formula(paste("~", group_var))
  dds <- phyloseq_to_deseq2(ps_sub, design_formula)

  # 3. Set reference
  colData(dds)[[group_var]] <- relevel(
    factor(colData(dds)[[group_var]]),
    ref = group_A
  )
  
  # 4. Run DESeq2
  if(fit_type %in% c("parametric", "local")) {
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds, fitType = fit_type)
    dds <- nbinomWaldTest(dds)
  } else {
    stop("Invalid fit_type. Use 'parametric' or 'local'.")
  }

  if(plot_disp) {
    pdf(file.path(working_dir, dada2_figure_out, paste0("deseq2_dispersion_plot_", group_B, "_vs_", group_A, ".pdf")))
    plotDispEsts(dds)
    dev.off()
  }

  # 5. Results
  res <- results(dds, contrast = c(group_var, group_B, group_A))
  
  # 6. Shrinkage
  if (shrink) {
    coef_name <- paste0(group_var, "_", group_B, "_vs_", group_A)
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    res_shrunk <- NULL
  }
  
  return(list(dds = dds, results = res, results_shrunk = res_shrunk))
}

res_list <- run_deseq2_phyloseq(
  ps,
  group_var = "treatment",
  group_A = "Control",
  group_B = "ASS",
  min_count = 10,
  shrink = TRUE,
  plot_disp = TRUE,
  fit_type = "local"
)
res <- res_list$results
res_shrunk_ASS_CONT <- res_list$results_shrunk

res_list <- run_deseq2_phyloseq(
  ps,
  group_var = "treatment",
  group_A = "ASS",
  group_B = "GLP1_ASS",
  min_count = 10,
  shrink = TRUE,
  plot_disp = TRUE,
  fit_type = "local"
)
res <- res_list$results
res_shrunk_GLP1ASS_ASS <- res_list$results_shrunk

res_list <- run_deseq2_phyloseq(
  ps,
  group_var = "treatment",
  group_A = "ASS",
  group_B = "ASS_DOMP",
  min_count = 10,
  shrink = TRUE,
  plot_disp = TRUE,
  fit_type = "local"
)
res <- res_list$results
res_shrunk_ASSDOMP_ASS <- res_list$results_shrunk

res1 <- res_shrunk_ASS_CONT
res2 <- res_shrunk_GLP1ASS_ASS
res3 <- res_shrunk_ASSDOMP_ASS

# sig1 <- rownames(res1)[which(res1$padj < 0.05)]
# sig2 <- rownames(res2)[which(res2$padj < 0.05)]
# sig3 <- rownames(res3)[which(res3$padj < 0.05)]
# all_sig_asvs <- unique(c(sig1, sig2, sig3))

library(dplyr)

make_df <- function(res, label) {
  df <- as.data.frame(res)
  df$ASV <- rownames(df)
  df$comparison <- label
  return(df)
}

df_all <- bind_rows(
  make_df(res1, "ASS_vs_Control"),
  make_df(res2, "GLP1-ASS_vs_ASS"),
  make_df(res3, "ASS-DOMP_vs_ASS")
)

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
# only show significant points
# df_all <- df_all %>% filter(sig_global)
# join with taxonomic information
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

df_all <- left_join(df_all, tax[, c("ASV", "label")], by = "ASV")
df_all$label <- gsub("_NA", "", df_all$label)
df_all$label <- paste0(df_all$label, "|", df_all$ASV)

write.csv(df_all,
          file.path(working_dir, dada2_data_out, "deseq2_results_all_comparisons.csv"),
          row.names = FALSE)

# only keep ASVs which are sig in at least one comparison
all_sig_asvs <- unique(df_all$ASV[df_all$sig_global == TRUE])
df_all <- df_all[df_all$ASV %in% all_sig_asvs, ]

write.csv(df_all,
          file.path(working_dir, dada2_data_out, "deseq2_results_all_comparisons_sig.csv"),
          row.names = FALSE)


p <- ggplot(df_all, aes(x = comparison, y = label)) +
  geom_point(aes(size = -log10(padj),
                 color = log2FoldChange)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() +
  labs(size = "-log10(padj)", color = "log2FC")


p <- ggplot(df_all, aes(x = log2FoldChange, y = label)) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.6, color = "black") +
  geom_point(aes(size = -log10(padj_global), color = sig_global)) +
  facet_grid(. ~ comparison) +
  theme_bw() +
  labs(
    x = "log2 fold change",
    y = "ASV",
    size = "-log10(global FDR)",
    color = "log2FC"
  )


# Group by direction
# Make sure the order of facet matches the order of the +/- pattern
df_all$comparison <- factor(df_all$comparison,
                             levels = c("ASS_vs_Control", "GLP1-ASS_vs_ASS", "ASS-DOMP_vs_ASS"))

df_sign <- df_all %>%
  mutate(sign = ifelse(log2FoldChange > 0, "+", "-")) %>%
  select(ASV, comparison, sign) %>%
  distinct()

df_pattern <- df_sign %>%
  pivot_wider(names_from = comparison,
              values_from = sign)

df_pattern$pattern <- paste(df_pattern$ASS_vs_Control,
                            df_pattern$`GLP1-ASS_vs_ASS`,
                            df_pattern$`ASS-DOMP_vs_ASS`,
                            sep = "")
# convert all patterns like -NA+ to NA
df_pattern$pattern <- ifelse(
  is.na(df_pattern$ASS_vs_Control) |
  is.na(df_pattern$`GLP1-ASS_vs_ASS`) |
  is.na(df_pattern$`ASS-DOMP_vs_ASS`),
  NA,
  df_pattern$pattern
)

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

ggsave(file.path(working_dir, dada2_figure_out, "deseq2_dotplot_all_comparisons.pdf"),
       plot = p, device = "pdf",
       width = 12, height = 15)