library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)

df <- read.table("/data/local/jy1008/project_folder/test_data_output/dada2_stats.tsv", header = TRUE, sep = "\t")
df$sample_name = rownames(df)

# bar plot with more metrics
# Pivot longer: convert counts into "metric"/"value"
df_long <- df %>%
  pivot_longer(
    cols = c(input, filtered, merged, nonchim), 
    names_to = "metric", 
    values_to = "count"
  )

# Make sure df_long$metric is ordered
df_long$metric <- factor(
  df_long$metric,
  levels = c("input", "filtered", "merged", "nonchim")
)

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

# Order samples by subject, then time
# df_long <- df_long %>%
  # Extract leading number as numeric
  # mutate(sample_num = as.numeric(str_extract(sample_name, "^\\d+"))) %>%
  # Arrange by that number
#   arrange(sample_order)
df_long <- df_long %>%
  arrange(factor(sample_name, levels = sample_order))
df_long <- df_long %>%
  mutate(sample_name = factor(sample_name, levels = unique(sample_name)))

# Stacked bar plot
p <- ggplot(df_long, aes(x = sample_name, y = count, fill = metric)) +
  geom_col(position = "identity", alpha = 1.0) +
  geom_hline(yintercept = 30000, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  labs(x = "Sample", y = "Read Count", fill = "Processing Step") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("dada2_filtering_read_count_bars.pdf", plot = p, width = 12, height = 4)

# only show input and nonchim
df_long_filtered <- df_long %>%
  filter(metric %in% c("input", "nonchim"))
df_long_filtered$metric_label <- ifelse(df_long_filtered$metric == "input", "Input Reads", "post-QC Reads")
p2 <- ggplot(df_long_filtered, aes(x = sample_name, y = count, fill = metric_label)) +
  geom_col(position = "identity", alpha = 1.0) +
  geom_hline(yintercept = 30000, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  labs(x = "Sample", y = "Read Count", fill = "Processing Step") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("dada2_filtering_read_count_bars_v2.pdf", plot = p2, width = 12, height = 4)



qc_df <- read.csv("/data/local/jy1008/project_folder/test_data_output/read_quality_summary.csv",
                  header = TRUE, stringsAsFactors = FALSE)
p <- ggplot(qc_df, aes(x = Cycle)) +
  geom_line(aes(y = Mean_All), color = "blue", size = 1) +                # mean line
  geom_ribbon(aes(ymin = Mean_All - SD_All, ymax = Mean_All + SD_All),   # shaded area
              fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  expand_limits(y = 20) +  # force y-axis to start at 0
  labs(x = "Cycle", y = "Quality Score", title = "Mean Quality ± SD Across Cycles")
ggsave("dada2_read_quality_across_cycles.pdf", plot = p, width = 6, height = 4)