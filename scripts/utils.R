# utils.R

library(ShortRead)
library(scales)
library(dplyr)
library(ComplexHeatmap)
library(circlize)  # for colorRamp2

summarize_quality_by_position <- function(fastq_files, output_csv = "combined_quality_summary.csv") {
  quality_list <- list()
  
  for (fq in fastq_files) {
    cat("Processing:", fq, "\n")
    
    fq_reads <- readFastq(fq)

    qs_matrix <- as(quality(fq_reads), "matrix")
    mean_q <- colMeans(qs_matrix)

    # Clean base name
    fq_base <- tools::file_path_sans_ext(basename(fq))
    fq_base <- sub("\\.fastq$", "", fq_base)
    fq_base <- sub("\\.fq$", "", fq_base)
    fq_base <- sub("\\.gz$", "", fq_base)

    quality_list[[fq_base]] <- round(mean_q, 2)
  }

  # Combine into data frame
  quality_df <- as.data.frame(quality_list)
  quality_df$Cycle <- seq_len(nrow(quality_df))
  quality_df <- quality_df[, c("Cycle", setdiff(names(quality_df), "Cycle"))]

  # Mean and standard deviation across samples
  quality_matrix <- as.matrix(quality_df[, -1])
  quality_df$Mean_All <- round(rowMeans(quality_matrix), 2)
  quality_df$SD_All <- round(apply(quality_matrix, 1, sd), 2)

  # Write to CSV
  write.csv(quality_df, output_csv, row.names = FALSE)
  cat("Saved combined quality summary to:", output_csv, "\n")
  
  return(quality_df)
}



# print percent of reads assigned to each taxa level in each file
classification_rate_table <- function(seqtab, taxa) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Ensure taxa and seqtab match by ASV (column names)
  common_seqs <- intersect(colnames(seqtab), rownames(taxa))
  seqtab <- seqtab[, common_seqs]
  taxa <- taxa[common_seqs, ranks, drop = FALSE]
  
  n_samples <- nrow(seqtab)
  result <- matrix(NA, nrow = n_samples, ncol = length(ranks) + 1)
  rownames(result) <- rownames(seqtab)
  colnames(result) <- c("Total_Reads", ranks)

  for (i in seq_len(n_samples)) {
    sample_counts <- seqtab[i, ]
    total_reads <- sum(sample_counts)
    result[i, "Total_Reads"] <- total_reads
    
    for (rank in ranks) {
      assigned <- !is.na(taxa[, rank])
      assigned_reads <- sum(sample_counts[assigned])
      result[i, rank] <- round(100 * assigned_reads / total_reads, 2)
    }
  }

  return(as.data.frame(result))
}

# Function to compute percent composition at a taxa level
taxa_percent_table <- function(seqtab, taxa, rank = "Genus") {
  # Check the rank exists
  if (!rank %in% colnames(taxa)) {
    stop(paste("Rank", rank, "not found in taxonomy table."))
  }
  
  # Align ASVs between seqtab and taxa
  common_seqs <- intersect(colnames(seqtab), rownames(taxa))
  seqtab <- seqtab[, common_seqs]
  taxa <- taxa[common_seqs, , drop = FALSE]
  
  # Replace NA with "Unassigned" in the selected rank
  tax_level <- taxa[, rank]
  tax_level[is.na(tax_level)] <- "Unassigned"
  
  # Sum reads by taxon at the chosen rank
  tax_counts <- t(rowsum(t(seqtab), group = tax_level))
  
  # Convert to percent per sample
  tax_pct <- tax_counts / rowSums(tax_counts) * 100
  tax_pct <- round(tax_pct, 2)
  
  return(as.data.frame(tax_pct))
}

filter_taxa_by_prevalence <- function(pct_table, min_abundance = 1, min_prevalence = 0.2) {
  # Remove Sample column if present
  sample_col <- "Sample" %in% colnames(pct_table)
  if (sample_col) {
    samples <- pct_table$Sample
    pct_table <- pct_table[, !colnames(pct_table) %in% "Sample"]
  }
  
  # Compute in how many samples each taxon exceeds min_abundance
  prevalence <- colMeans(pct_table > min_abundance)

  # Keep taxa that are above threshold in at least min_prevalence fraction of samples
  keep_taxa <- names(prevalence[prevalence >= min_prevalence])

  # Reconstruct filtered table
  filtered <- pct_table[, keep_taxa, drop = FALSE]

  # Restore Sample column if it was present
  if (sample_col) {
    filtered$Sample <- samples
  }

  return(filtered)
}

plot_taxa_heatmap <- function(pct_table,
                              min_abundance = 1,
                              min_prevalence = 0.2,
                              scale_rows = FALSE,
                              tax_rank = "Genus",
                              cluster_rows = FALSE,
                              cluster_cols = TRUE) {

  # Extract and drop sample column if present
  has_sample <- "Sample" %in% colnames(pct_table)
  if (has_sample) {
    samples <- pct_table$Sample
    pct_table <- pct_table[, !colnames(pct_table) %in% "Sample"]
    rownames(pct_table) <- samples
  }

  # Filter taxa based on abundance and prevalence
  prevalence <- colMeans(pct_table > min_abundance)
  keep_taxa <- names(prevalence[prevalence >= min_prevalence])
  # Always include "Unassigned" if present
  if ("Unassigned" %in% colnames(pct_table)) {
    keep_taxa <- union(keep_taxa, "Unassigned")
  }
  filtered <- pct_table[, keep_taxa, drop = FALSE]

  # Optionally scale rows (z-score normalization)
  if (scale_rows) {
    filtered <- t(scale(t(filtered)))
  }

  # Melt for ggplot
  # df_long <- melt(filtered, varnames = c("Sample", "Taxon"), value.name = "Abundance")
  filtered$Sample <- rownames(filtered)
  filtered$Sample <- sub("^CP101", "", filtered$Sample)
  filtered$Sample <- sub("^(.{2})(.)", "\\1 \\2", filtered$Sample)

  df_long <- melt(filtered, id.vars = "Sample", variable.name = "Taxon", value.name = "Abundance")
  colnames(df_long)[2] <- tax_rank  # Rename for prettier axis label

  # Plot heatmap
  gg <- ggplot(df_long, aes(x = Sample, y = .data[[tax_rank]], fill = Abundance)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = if (scale_rows) 0 else median(df_long$Abundance, na.rm = TRUE),
                         name = if (scale_rows) "Z-score" else "Percent") +
    scale_x_discrete(position = "top") +
    labs(x = "Sample", y = tax_rank,
         title = paste("Heatmap of", tax_rank, "Abundance")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0),
          axis.text.y = element_text(size = 8))

  # Optional clustering (sort rows/columns by hierarchical clustering)
  if (cluster_rows) {
    filtered_clustering <- filtered %>% select(-Sample, Unassigned)
    taxa_order <- hclust(dist(t(filtered_clustering)))$order
    taxa_names_ordered <- colnames(filtered_clustering)[taxa_order]
    taxa_names_ordered <- c(taxa_names_ordered, "Unassigned")
    gg <- gg + scale_y_discrete(limits = taxa_names_ordered)
  }
  if (cluster_cols) {
    sample_order <- hclust(dist(filtered))$order
    gg <- gg + scale_x_discrete(limits = rownames(filtered)[sample_order])
  }

  return(gg)
}


# Function to plot ASV heatmap with optional filtering
# seqtab_filtered: ASV table with samples as rows and ASVs as columns
# meta_data: metadata for samples (partial mayo, treatment)
# min_abundance: minimum read count for an ASV to be included
# min_prevalence: minimum fraction of samples an ASV must be present in
# scale_rows: whether to z-score normalize rows (ASVs)
# cluster_rows: whether to cluster ASVs by hierarchical clustering
# cluster_cols: whether to cluster samples by hierarchical clustering
# Returns a ggplot object
# Note: Requires ggplot2, reshape2, ComplexHeatmap
# If ComplexHeatmap is not available, it falls back to ggplot2 heatmap
plot_asv_heatmap <- function(seqtab_filtered,
                             meta_data,
                             min_abundance = 1,
                             min_prevalence = 0.2,
                             rel_abund = TRUE,
                             cluster_rows = TRUE,
                             cluster_cols = FALSE) {
    library(reshape2)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)

    # Ensure ASV IDs are columns and sample names are rownames
    asv_table <- as.data.frame(seqtab_filtered)
    asv_table$Sample <- rownames(asv_table)

    # Filter ASVs based on prevalence
    prevalence <- colMeans(asv_table[, !colnames(asv_table) %in% "Sample"] > min_abundance)
    keep_asvs <- names(prevalence[prevalence >= min_prevalence])
    filtered <- asv_table[, c(keep_asvs, "Sample"), drop = FALSE]

    # Join with metadata
    meta <- meta_data[match(filtered$Sample, meta_data$sample_id), ]
    # Reorder by treatment group if clustering is off
    if (!cluster_cols && "treatment" %in% colnames(meta)) {
        treatment_order <- order(meta$treatment)
        filtered <- filtered[treatment_order, ]
        meta <- meta[treatment_order, ]
    }

    # Optionally scale rows (Z-score by ASV)
    mat <- as.matrix(filtered[, keep_asvs])
    # if (scale_rows) {
        # mat <- t(scale(t(mat)))
    if (rel_abund) {
        # mat <- sweep(mat, 1, rowSums(mat), FUN = "/")  # convert to relative abundance
        # Avoid division by zero
        row_totals <- rowSums(mat)
        safe_row_totals <- ifelse(row_totals == 0, 1, row_totals)
        mat <- sweep(mat, 1, safe_row_totals, FUN = "/")
        mat[row_totals == 0, ] <- 0
    }
    rownames(mat) <- filtered$Sample

    # If clustering is enabled, use ComplexHeatmap
    if (cluster_rows || cluster_cols) {
        col_fun <- circlize::colorRamp2(
            breaks = c(min(mat), max(mat)),
            colors = c("white", "red")
        )

        # Define treatment color
        treatment_colors <- c("placebo" = "#1f77b4", "active" = "#ff7f0e")

        # Define annotation
        bottom_annot <- HeatmapAnnotation(
            Treatment = meta$treatment,
            Mayo = meta$total_partial_mayo,
            col = list(
                Treatment = treatment_colors,
                Mayo = circlize::colorRamp2(
                    range(meta$total_partial_mayo, na.rm = TRUE),
                    c("white", "darkgreen")
                )
            ),
            which = "column",
            annotation_name_side = "left",
            annotation_legend_param = list(
                Treatment = list(title = "Treatment"),
                Mayo = list(title = "Mayo Score")
            )
        )

        ht <- Heatmap(
            t(mat),
            name = if (rel_abund) "Rel Abundance" else "Abundance",
            col = col_fun,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = cluster_rows,
            show_column_dend = cluster_cols,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 90,
            heatmap_legend_param = list(direction = "horizontal"),
            top_annotation = NULL,
            bottom_annotation = bottom_annot
        )
        return(ht)
    }

    # Otherwise, fall back to ggplot2 heatmap
    filtered$Sample <- sub("^CP101", "", filtered$Sample)
    filtered$Sample <- sub("^(.{2})(.)", "\\1 \\2", filtered$Sample)

    df_long <- melt(filtered, id.vars = "Sample",
                    variable.name = "ASV", value.name = "Abundance")

    gg <- ggplot(df_long, aes(x = Sample, y = ASV, fill = Abundance)) +
        geom_tile() +
        scale_fill_gradient(
            low = "white", high = "red",
            name = if (rel_abund) "Rel Abund" else "Abundance"
        ) +
        scale_x_discrete(position = "top") +
        labs(x = "Sample", y = "ASV", title = "Heatmap of ASV Abundance") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0),
            axis.text.y = element_text(size = 6)
        )

    return(gg)
}



plot_asv_heatmap2 <- function(seqtab_filtered,
                             meta_data,
                             min_abundance = 1,
                             min_prevalence = 0.2,
                             rel_abund = TRUE,
                             cluster_rows = TRUE,
                             cluster_cols = FALSE) {
    library(reshape2)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(stringr)

    # Ensure ASV IDs are columns and sample names are rownames
    asv_table <- as.data.frame(seqtab_filtered)
    asv_table$Sample <- rownames(asv_table)

    # Filter ASVs based on prevalence
    prevalence <- colMeans(asv_table[, !colnames(asv_table) %in% "Sample"] > min_abundance)
    keep_asvs <- names(prevalence[prevalence >= min_prevalence])
    filtered <- asv_table[, c(keep_asvs, "Sample"), drop = FALSE]
    filtered$SampleSort <- gsub("Day1", "Wk0", gsub("CP101", "", filtered$Sample))

    # Join with metadata
    meta <- meta_data[match(filtered$Sample, meta_data$sample_id), ]

    # Reorder samples: treatment first, then numeric order from sample name
    if (!cluster_cols && "treatment" %in% colnames(meta)) {
        # Extract numeric part from sample name and time point
        all_nums <- str_extract_all(filtered$SampleSort, "\\d+")
        sample_num <- sapply(all_nums, function(v) if(length(v) >= 1) as.numeric(v[1]) else NA)
        time_num <- sapply(all_nums, function(v) if(length(v) >= 2) as.numeric(v[2]) else NA)

        treatment_order <- order(meta$treatment, sample_num, time_num)
        filtered <- filtered[treatment_order, ]
        meta <- meta[treatment_order, ]
    }

    # Optionally scale to relative abundance
    mat <- as.matrix(filtered[, keep_asvs])
    if (rel_abund) {
        row_totals <- rowSums(mat)
        safe_row_totals <- ifelse(row_totals == 0, 1, row_totals)
        mat <- sweep(mat, 1, safe_row_totals, FUN = "/")
        mat[row_totals == 0, ] <- 0
    }
    rownames(mat) <- filtered$Sample

    # If clustering is enabled, use ComplexHeatmap
    if (cluster_rows || cluster_cols) {
        col_fun <- circlize::colorRamp2(
            breaks = c(min(mat), max(mat)),
            colors = c("white", "red")
        )

        treatment_colors <- c("placebo" = "#1f77b4", "active" = "#ff7f0e")

        bottom_annot <- HeatmapAnnotation(
            Treatment = meta$treatment,
            Mayo = meta$total_partial_mayo,
            col = list(
                Treatment = treatment_colors,
                Mayo = circlize::colorRamp2(
                    range(meta$total_partial_mayo, na.rm = TRUE),
                    c("white", "darkgreen")
                )
            ),
            which = "column",
            annotation_name_side = "left",
            annotation_legend_param = list(
                Treatment = list(title = "Treatment"),
                Mayo = list(title = "Mayo Score")
            )
        )

        ht <- Heatmap(
            t(mat),
            name = if (rel_abund) "Rel Abundance" else "Abundance",
            col = col_fun,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = cluster_rows,
            show_column_dend = cluster_cols,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 90,
            heatmap_legend_param = list(direction = "horizontal"),
            top_annotation = NULL,
            bottom_annotation = bottom_annot
        )
        return(ht)
    }

    # Fallback: ggplot2 heatmap
    df_long <- melt(filtered, id.vars = "Sample",
                    variable.name = "ASV", value.name = "Abundance")

    gg <- ggplot(df_long, aes(x = factor(Sample, levels = unique(filtered$Sample)),
                              y = ASV, fill = Abundance)) +
        geom_tile() +
        scale_fill_gradient(
            low = "white", high = "red",
            name = if (rel_abund) "Rel Abund" else "Abundance"
        ) +
        scale_x_discrete(position = "top") +
        labs(x = "Sample", y = "ASV", title = "Heatmap of ASV Abundance") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0),
            axis.text.y = element_text(size = 6)
        )

    return(gg)
}




plot_qc_elbow <- function(seqtab_filtered) {
  # Convert counts to relative abundance per sample
  rel_abund <- sweep(seqtab_filtered, 1, rowSums(seqtab_filtered), FUN = "/")

  # Get 2nd highest relative abundance per ASV
  second_highest <- apply(rel_abund, 2, function(x) {
    sort(x, decreasing = TRUE)[2]  # second largest value
  })

  # Build data frame for plotting
  df_elbow <- data.frame(
    ASV = names(second_highest),
    second_highest = second_highest
  ) %>%
    arrange(desc(second_highest)) %>%
    mutate(rank = row_number())

  # Elbow plot
  p1 <- ggplot(df_elbow, aes(x = rank, y = second_highest)) +
    geom_line() +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 0.001, linetype = "dashed", color = "red") +
    theme_bw() +
    labs(x = "ASV Rank", y = "Second-highest relative abundance",
        title = "Elbow plot for ASV filtering")
  ggsave("dada2_filtering_elbow_second_abund.pdf", plot = p1, width = 6, height = 4)
  p1 <- ggplot(df_elbow, aes(x = rank, y = second_highest)) +
    geom_line() +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 0.001, linetype = "dashed", color = "red") +
    scale_y_log10() +
    theme_bw() +
    labs(x = "ASV Rank", y = "Second-highest relative abundance (log scale)",
        title = "Elbow plot for ASV filtering")
  ggsave("dada2_filtering_elbow_second_abund_log.pdf", plot = p1, width = 6, height = 4)

  # Prevalence = number of samples with nonzero abundance
  prevalence <- colSums(seqtab_filtered > 0)

  # Total abundance = sum across all samples
  total_abundance <- colSums(seqtab_filtered)

  # ---- Prevalence elbow plot ----
  prev_df <- data.frame(
    ASV = colnames(seqtab_filtered),
    prevalence = prevalence
  ) %>%
    arrange(desc(prevalence)) %>%
    mutate(rank = row_number())

  p2 <- ggplot(prev_df, aes(x = rank, y = prevalence)) +
    geom_line() +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(x = "ASVs (ranked)", y = "Prevalence (# samples)")
  ggsave("dada2_filtering_elbow_prev.pdf", plot = p2, width = 6, height = 4)

  # ---- Total abundance elbow plot ----
  abund_df <- data.frame(
    ASV = colnames(seqtab_filtered),
    total_abundance = total_abundance
  ) %>%
    arrange(desc(total_abundance)) %>%
    mutate(rank = row_number())

  p3 <- ggplot(abund_df, aes(x = rank, y = total_abundance)) +
    geom_line() +
    geom_point(size = 0.5) +
    theme_bw() +
    labs(x = "ASVs (ranked)", y = "Total abundance (reads)")
  ggsave("dada2_filtering_elbow_total_count.pdf", plot = p3, width = 6, height = 4)
  p3 <- ggplot(abund_df, aes(x = rank, y = total_abundance)) +
    geom_line() +
    geom_point(size = 0.5) +
    scale_y_log10() +   # log-scale helps spread rare ASVs
    theme_bw() +
    labs(x = "ASVs (ranked)", y = "Total abundance (reads)")
  ggsave("dada2_filtering_elbow_total_count_log.pdf", plot = p3, width = 6, height = 4)

}


plot_qc_filtering_subjects <- function(seqtab, meta_data, subjects_col,
                                       relabund_thresh_range = seq(0, 0.05, length.out = 100),
                                       min_subjects = c(1, 2, 3)) {
  # seqtab: samples x ASVs
  # meta_data: dataframe with sample info, must have sample IDs
  # subjects_col: column in meta_data that identifies subject ID for each sample
  
  # Step 1: relative abundances per sample
  relabund <- sweep(seqtab, 1, rowSums(seqtab), FUN = "/")
  
  # Map sample -> subject
  subjects <- meta_data[[subjects_col]][match(rownames(relabund), meta_data[[subjects_col]])]
  
  # Prepare result data frame
  results <- expand.grid(
    rel_abund_thresh = relabund_thresh_range,
    min_subjects = min_subjects
  )
  results$n_ASVs <- 0
  
  # Step 2: loop over thresholds and min_subjects
  for (i in seq_len(nrow(results))) {
    thresh <- results$rel_abund_thresh[i]
    min_subj <- results$min_subjects[i]
    
    # For each ASV, check how many subjects have >= thresh in any sample
    asv_pass <- sapply(1:ncol(relabund), function(j) {
      abund_j <- relabund[, j]
      subj_pass <- tapply(abund_j, subjects, function(x) any(x >= thresh))
      sum(subj_pass) >= min_subj
    })
    
    results$n_ASVs[i] <- sum(asv_pass)
  }
  
  # Step 3: plot
  p <- ggplot(results, aes(x = rel_abund_thresh, y = n_ASVs, color = factor(min_subjects))) +
    geom_line(size = 1) +
    # scale_x_log10() +
    theme_bw() +
    labs(
      x = "Relative abundance threshold",
      y = "Number of ASVs passing filter",
      color = "Min subjects"
    ) # +
    # ggtitle("ASV filtering summary (subject-level)")
  
  # ggsave("dada2_filtering_subject_summary.pdf", plot = p, width = 6, height = 4)
  return(p)
}



rel_abund_filter_subject <- function(seqtab, meta_data, subjects_col,
                                     relabund_thresh = 0.001,
                                     n_subjects_thresh = 2) {
  # Check inputs
  if (is.null(rownames(seqtab))) stop("seqtab must have sample names as rownames")
  if (!(subjects_col %in% colnames(meta_data))) stop(paste("Column", subjects_col, "not found in meta_data"))
  
  # Step 1: compute relative abundance per sample
  relabund <- sweep(seqtab, 1, rowSums(seqtab), FUN = "/")
  
  # Step 2: map samples to subjects (same order as seqtab rows)
  subjects <- meta_data[match(rownames(relabund), rownames(meta_data)), subjects_col]
  
  # Handle missing or unmatched samples
  if (any(is.na(subjects))) {
    warning("Some samples in seqtab were not found in meta_data; those rows will be ignored")
    keep <- !is.na(subjects)
    relabund <- relabund[keep, , drop = FALSE]
    subjects <- subjects[keep]
  }

  # Step 3: split samples by subject
  samples_by_subject <- split(seq_len(nrow(relabund)), subjects)
  # This gives a list: names = subject IDs, values = row indices of that subject

  # Step 4: for each ASV, check which subjects have >= threshold in any sample
  asv_subject_pass <- matrix(FALSE, nrow = length(unique(subjects)), ncol = ncol(relabund),
                            dimnames = list(unique(subjects), colnames(relabund)))

  for (asv_idx in seq_len(ncol(relabund))) {
    for (subject_id in names(samples_by_subject)) {
      sample_indices <- samples_by_subject[[subject_id]]
      # any sample for this subject meets threshold?
      asv_subject_pass[subject_id, asv_idx] <- any(relabund[sample_indices, asv_idx] >= relabund_thresh)
    }
  }

  # Step 5: count subjects meeting threshold per ASV
  subjects_per_asv <- colSums(asv_subject_pass)

  # Step 6: filter ASVs that are present in >= n_subjects_thresh subjects
  pass_asvs <- subjects_per_asv >= n_subjects_thresh
  seqtab_filtered <- seqtab[, pass_asvs, drop = FALSE]
    
  message(sum(pass_asvs), " of ", ncol(seqtab), " ASVs retained (>= ", 
          relabund_thresh, " in >= ", n_subjects_thresh, " subjects)")
  return(seqtab_filtered)
}

unify_data <- function(seqtab_list, meta_list) {
  # get union of all ASV columns
  all_asvs <- unique(unlist(lapply(seqtab_list, colnames)))

  # function to fill missing ASVs with 0
  fill_missing_asvs <- function(seqtab_sub, all_asvs) {
    # ensure matrix
    seqtab_sub <- as.matrix(seqtab_sub)
    
    # add missing ASVs as zeros
    missing_asvs <- setdiff(all_asvs, colnames(seqtab_sub))
    if(length(missing_asvs) > 0){
      # create a zero matrix for missing ASVs
      zero_mat <- matrix(0, nrow = nrow(seqtab_sub), ncol = length(missing_asvs),
                        dimnames = list(rownames(seqtab_sub), missing_asvs))
      seqtab_sub <- cbind(seqtab_sub, zero_mat)
    }

    # reorder columns
    seqtab_sub <- seqtab_sub[, all_asvs, drop = FALSE]
    return(seqtab_sub)
  }

  # apply to each tissue
  seqtab_filled <- lapply(seqtab_list, fill_missing_asvs, all_asvs = all_asvs)

  # combine rows (samples)
  seqtab_combined <- do.call(rbind, seqtab_filled)

  # combine metadata
  meta_combined <- do.call(rbind, meta_list)
  rownames(meta_combined) <- meta_combined$MHMC.sampleID

  # check alignment
  seqtab_combined <- seqtab_combined[rownames(meta_combined), , drop = FALSE]

  if(!all(rownames(seqtab_combined) == rownames(meta_combined))){
    stop("Sample names in seqtab_combined and meta_combined do not match!")
  }
  return(list(seqtab = seqtab_combined, meta = meta_combined))
}

assign_asv_ids <- function(seqtab, taxa, tax_col = "label") {
  # seqtab: samples x ASVs matrix
  # taxa: ASV taxonomy table, rownames = ASV sequences
  # tax_col: column in taxa with human-readable taxonomy

  # Original ASV sequences
  asv_seqs <- colnames(seqtab)
  
  # Assign short IDs: ASV1, ASV2, ...
  asv_ids <- paste0("ASV", seq_along(asv_seqs))
  names(asv_ids) <- asv_seqs  # Map sequence -> ASV ID
  
  # Replace column names in seqtab
  colnames(seqtab) <- asv_ids[colnames(seqtab)]
  
  # Build mapping table
  asv_mapping <- data.frame(
    ASV_ID = asv_ids,
    Sequence = names(asv_ids),
    Taxonomy = taxa[names(asv_ids), tax_col],
    stringsAsFactors = FALSE
  )
  
  # Optional helper function to combine taxonomy + ASV ID
  get_asv_label <- function(sequence) {
    matched <- asv_mapping[asv_mapping$Sequence == sequence, ]
    if (nrow(matched) == 0) return(NA)
    paste(matched$Taxonomy, matched$ASV_ID, sep = " | ")
  }
  
  return(list(
    seqtab = seqtab,
    mapping = asv_mapping,
    get_asv_label = get_asv_label
  ))
}


plot_family_heatmap <- function(seqtab, taxa, top_n = 15, family_colors = NULL,
                                out_file = NULL, sample_order = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)

  # --- Prepare seqtab dataframe ---
  seqtab_df <- as.data.frame(t(seqtab))
  seqtab_df$ASV <- rownames(seqtab_df)

  # --- Map ASVs to family ---
  taxa$ASV <- rownames(taxa)
  seqtab_taxa <- seqtab_df %>%
    left_join(taxa[, c("ASV", "Family")], by = "ASV")

  # --- Sum counts per family ---
  family_counts <- seqtab_taxa %>%
    group_by(Family) %>%
    summarise(across(-ASV, sum), .groups = "drop")

  # --- Total reads per family & cumulative percentage ---
  family_totals <- family_counts %>%
    rowwise() %>%
    mutate(TotalReads = sum(c_across(-Family))) %>%
    ungroup() %>%
    arrange(desc(TotalReads))

  # --- Top N families ---
  top_families <- family_totals$Family[1:min(top_n, nrow(family_totals))]

  # Collapse remaining into "Other"
  family_counts_collapsed <- family_counts %>%
    mutate(Family = ifelse(Family %in% top_families, Family, "Other"))

  # Pivot to long format & compute relative abundance
  df_long <- family_counts_collapsed %>%
    pivot_longer(-Family, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Sample) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>%
    ungroup()

  # Compute mean relative abundance for labeling
  family_means <- df_long %>%
    filter(Family != "Other") %>%
    group_by(Family) %>%
    summarise(mean_rel = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_rel))

  legend_labels <- c(
    setNames(
      paste0(family_means$Family, " (", percent(family_means$mean_rel, accuracy = 0.1), ")"),
      family_means$Family
    ),
    Other = "Other"
  )
  df_long$Family <- factor(df_long$Family,
                           levels = c(family_means$Family, "Other"))

  # Default colors if not provided
  if (is.null(family_colors)) {
    n_fam <- length(levels(df_long$Family))
    family_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_fam),
      levels(df_long$Family)
    )
  }

  # --- Order samples ---
  if (!is.null(sample_order)) {
    # ensure factor levels follow sample_order
    df_long$Sample <- factor(df_long$Sample, levels = sample_order)
  } else {
    # Default: order by numeric prefix
    df_long <- df_long %>%
      mutate(Sample_num = as.numeric(str_extract(Sample, "^\\d+"))) %>%
      arrange(Sample_num) %>%
      mutate(Sample = factor(Sample, levels = unique(Sample)))
  }

  # --- Plot ---
  p <- ggplot(df_long, aes(x = Sample, y = RelAbund, fill = Family)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = family_colors, labels = legend_labels) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Sample", y = "Relative abundance", fill = "Family") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_blank()
    )

  # Save if out_file provided
  if(!is.null(out_file)) {
    pdf(out_file, width = 16, height = 8)
    print(p)
    dev.off()
  }

  return(list(plot = p, df_long = df_long))
}


plot_family_heatmap_v2 <- function(seqtab, taxa,
                                top_n = 15,
                                min_rel_abund = 0,
                                family_colors = NULL,
                                out_file = NULL,
                                sample_order = NULL) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)

  # -------------------------
  # 1. ASV → Family counts
  # -------------------------
  seqtab_df <- as.data.frame(t(seqtab))
  seqtab_df$ASV <- rownames(seqtab_df)

  taxa$ASV <- rownames(taxa)

  family_counts <- seqtab_df %>%
    left_join(taxa[, c("ASV","Family")], by = "ASV") %>%
    group_by(Family) %>%
    summarise(across(-ASV, sum), .groups="drop")

  # -------------------------
  # 2. Long format + rel abundance
  # -------------------------
  df_long <- family_counts %>%
    pivot_longer(-Family, names_to="Sample", values_to="Abundance") %>%
    group_by(Sample) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>%
    ungroup()

  # -------------------------
  # 3. Family statistics
  # -------------------------
  family_stats <- df_long %>%
    group_by(Family) %>%
    summarise(
      mean_rel = mean(RelAbund),
      total_reads = sum(Abundance),
      .groups="drop"
    ) %>%
    arrange(desc(mean_rel))

  keep_families <- family_stats %>%
    filter(mean_rel >= min_rel_abund) %>%
    slice_head(n = top_n) %>%
    pull(Family)

  # -------------------------
  # 4. Collapse rare families
  # -------------------------
  df_long <- df_long %>%
    mutate(Family = ifelse(Family %in% keep_families, Family, "Other"))

  # recompute means after collapse
  family_means <- df_long %>%
    filter(Family != "Other") %>%
    group_by(Family) %>%
    summarise(mean_rel = mean(RelAbund), .groups="drop") %>%
    arrange(desc(mean_rel))

  df_long$Family <- factor(df_long$Family,
                           levels=c(family_means$Family,"Other"))

  # -------------------------
  # 5. Legend labels
  # -------------------------
  legend_labels <- c(
    setNames(
      paste0(family_means$Family," (",
             percent(family_means$mean_rel, accuracy=0.1),")"),
      family_means$Family
    ),
    Other="Other"
  )

  # -------------------------
  # 6. Colors
  # -------------------------
  if(is.null(family_colors)) {
    n <- length(levels(df_long$Family))
    family_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(n),
      levels(df_long$Family)
    )
  }

  # -------------------------
  # 7. Sample order
  # -------------------------
  if(!is.null(sample_order)) {

    df_long$Sample <- factor(df_long$Sample, levels=sample_order)

  } else {

    df_long <- df_long %>%
      mutate(Sample_num = as.numeric(str_extract(Sample,"^\\d+"))) %>%
      arrange(Sample_num) %>%
      mutate(Sample=factor(Sample,levels=unique(Sample)))

  }

  # -------------------------
  # 8. Plot
  # -------------------------
  p <- ggplot(df_long,
              aes(x=Sample,y=RelAbund,fill=Family)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=family_colors,
                      labels=legend_labels) +
    scale_y_continuous(labels=percent_format(accuracy=1)) +
    labs(x="Sample",y="Relative abundance",fill="Family") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
      panel.grid.major.x=element_blank()
    )

  if(!is.null(out_file)){
    pdf(out_file,width=16,height=8)
    print(p)
    dev.off()
  }

  return(list(plot=p, df_long=df_long, family_means=family_means))
}





plot_beta_diversity <- function(ps,
                                color_var = "Sample_type",
                                dist_method = "bray",
                                ordinate_method = "PCoA",
                                label_var = NULL,
                                output_path = NULL,
                                title = NULL) {
  # Load required packages
  library(phyloseq)
  library(ggplot2)
  library(dplyr)

  # 1 Ordination
  ordination <- ordinate(ps, method = ordinate_method, distance = dist_method)
  eig_vals <- ordination$values$Eigenvalues
  var_explained <- eig_vals / sum(eig_vals)
  pc1 <- round(var_explained[1] * 100, 1)
  pc2 <- round(var_explained[2] * 100, 1)

  # 2 Extract ordination dataframe with metadata
  ord_df <- plot_ordination(ps, ordination, justDF = TRUE)

  # Ensure chosen variables exist
  if (!color_var %in% colnames(ord_df)) {
    stop(paste("Column", color_var, "not found in sample_data(ps)"))
  }
  if (!is.null(label_var) && !label_var %in% colnames(ord_df)) {
    stop(paste("Column", label_var, "not found in sample_data(ps)"))
  }

  # 3 Build ggplot
  p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = .data[[color_var]])) +
    geom_point(alpha = 0.8, size = 3) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.margin = margin(t = 10),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(color = guide_legend(nrow = 1)) +
    labs(
      title = title %||% paste0(ordinate_method, " of ", dist_method, " Distances"),
      x = paste0("Axis 1 (", pc1, "%)"),
      y = paste0("Axis 2 (", pc2, "%)"),
      color = color_var
    )

  # Optionally add labels
  if (!is.null(label_var)) {
    p <- p +
      geom_text(aes(label = .data[[label_var]]), vjust = -0.8, size = 3)
  }

  # 4 Save if output_path provided
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 6, height = 6)
  }

  # Return ggplot object
  return(p)
}
