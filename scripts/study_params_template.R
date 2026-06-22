# study_params.R
#
# Per-study parameters for the downstream analysis/plotting scripts
# (fig2_asv_heatmap.R, dada2_additional_figures.R).
#
# Copy this file, rename it (e.g. study_params_<studyname>.R), fill in
# your study's actual values, and point to it via the `study_params`
# field in your config .ini file.

# Order in which samples should appear in plots (x-axis / heatmap columns).
# Must match the sample_id values in your metadata file.
sample_order <- c(
  "Sample1", "Sample2", "Sample3"
  # ...
)

# Factor levels for the treatment/condition column in your metadata,
# in the order you want them displayed (e.g. on alpha diversity boxplots).
treatment_levels <- c(
  "Control", "Treatment1", "Treatment2"
  # ...
)

# Pairwise DESeq2 comparisons to run, in the order you want them
# displayed across plot facets.
#   group_A = reference/baseline level
#   group_B = comparison level
#   label   = name used in output filenames, CSV "comparison" column,
#             and plot facet labels
comparisons <- list(
  list(group_A = "Control", group_B = "Treatment1", label = "Treatment1_vs_Control")
  # list(group_A = "Treatment1", group_B = "Treatment2", label = "Treatment2_vs_Treatment1"),
  # ...
)

# Sample IDs to exclude from diversity/differential-abundance analyses
# (e.g. technical replicates or failed samples). Leave as c() if none.
samples_to_remove <- c()

# Colors for taxa (e.g. bacterial families) in stacked bar / heatmap
# plots. Any taxon not listed here falls back to "Other".
family_colors <- c(
  "Other" = "#808080"
)

# QC thresholds plotted as dashed reference lines in dada2_additional_figures.R
read_count_threshold <- 30000   # read-count bar plot
quality_score_threshold <- 30   # quality-by-cycle plot