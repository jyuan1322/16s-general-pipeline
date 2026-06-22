#!/usr/bin/env Rscript
# dada2_additional_figures.R
#
# Supplementary QC plots: read counts through the pipeline, and quality
# scores across read cycles. Paths come from a config .ini file; sample
# order and QC thresholds come from a per-study R file referenced inside
# that config (see study_params_template.R).
#
# Usage: Rscript dada2_additional_figures.R <path_to_config_file>

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)

source("utils.R")  # provides read_config()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript dada2_additional_figures.R <path_to_config_file>")
}

config <- read_config(args[1], required_fields = c(
  "data_output_dir", "figure_output_dir", "study_params"
))

dada2_data_out <- config$data_output_dir
dada2_figure_out <- config$figure_output_dir
dir.create(dada2_figure_out, recursive = TRUE, showWarnings = FALSE)

# Per-study parameters: sample_order, read_count_threshold, quality_score_threshold
source(config$study_params)

df <- read.table(file.path(dada2_data_out, "dada2_stats.tsv"), header = TRUE, sep = "\t")
df$sample_name <- rownames(df)

# bar plot with more metrics
df_long <- df %>%
  pivot_longer(
    cols = c(input, filtered, merged, nonchim),
    names_to = "metric",
    values_to = "count"
  )

df_long$metric <- factor(
  df_long$metric,
  levels = c("input", "filtered", "merged", "nonchim")
)

# Order samples by the order given in study_params.R
df_long <- df_long %>%
  arrange(factor(sample_name, levels = sample_order)) %>%
  mutate(sample_name = factor(sample_name, levels = unique(sample_name)))

# Stacked bar plot
p <- ggplot(df_long, aes(x = sample_name, y = count, fill = metric)) +
  geom_col(position = "identity", alpha = 1.0) +
  geom_hline(yintercept = read_count_threshold, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  labs(x = "Sample", y = "Read Count", fill = "Processing Step") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dada2_figure_out, "dada2_filtering_read_count_bars.pdf"), plot = p, width = 12, height = 4)

# only show input and nonchim
df_long_filtered <- df_long %>%
  filter(metric %in% c("input", "nonchim"))
df_long_filtered$metric_label <- ifelse(df_long_filtered$metric == "input", "Input Reads", "post-QC Reads")
p2 <- ggplot(df_long_filtered, aes(x = sample_name, y = count, fill = metric_label)) +
  geom_col(position = "identity", alpha = 1.0) +
  geom_hline(yintercept = read_count_threshold, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  labs(x = "Sample", y = "Read Count", fill = "Processing Step") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dada2_figure_out, "dada2_filtering_read_count_bars_v2.pdf"), plot = p2, width = 12, height = 4)

qc_df <- read.csv(file.path(dada2_data_out, "read_quality_summary.csv"),
                  header = TRUE, stringsAsFactors = FALSE)
p <- ggplot(qc_df, aes(x = Cycle)) +
  geom_line(aes(y = Mean_All), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Mean_All - SD_All, ymax = Mean_All + SD_All),
              fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = quality_score_threshold, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  expand_limits(y = 20) +
  labs(x = "Cycle", y = "Quality Score", title = "Mean Quality \u00b1 SD Across Cycles")
ggsave(file.path(dada2_figure_out, "dada2_read_quality_across_cycles.pdf"), plot = p, width = 6, height = 4)