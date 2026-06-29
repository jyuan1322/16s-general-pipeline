#!/usr/bin/env Rscript
# .libPaths("/data/bwh-comppath-seq/software/R-4.1.0//library/")
# run micromamba activate r-microbiome

# NOTE: compare to:
# https://github.com/gerberlab/cdiff_paper_analyses/blob/master/scripts/dada2/dada2.R
# Big-data adaptation reference: https://benjjneb.github.io/dada2/bigdata_paired.html

# Load required libraries
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)
library(pheatmap)
library(ggpubr)

library(ShortRead)
source("utils.R")

# Function to read config file
read_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  print(config_path)
  config <- read.table(config_path,
                      sep = "=",
                      header = FALSE,
                      row.names = 1,
                      strip.white = TRUE,
                      stringsAsFactors = FALSE,
                      quote = "\"",
                      comment.char = "#")

  # Convert to named list
  config_list <- as.list(config$V2)
  names(config_list) <- rownames(config)

  # Validate required fields
  required_fields <- c("data_path", "pattern_forward", "pattern_reverse",
                      "data_output_dir", "figure_output_dir",
                      "truncLen_f", "truncLen_r", "taxonomy_db", "taxonomy_db_species")

  missing_fields <- required_fields[!required_fields %in% names(config_list)]
  if (length(missing_fields) > 0) {
    stop("Missing required fields in config file: ",
         paste(missing_fields, collapse = ", "))
  }

  return(config_list)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript process_16s.R <path_to_config_file>")
}

# Read and validate config
config <- read_config(args[1])

# Create output directories if they don't exist
dir.create(config$data_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$figure_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$data_path, "filtered"), recursive = TRUE, showWarnings = FALSE)

print("Read config and paths created")

# List input files
fnFs <- sort(list.files(config$data_path,
                       pattern=config$pattern_forward,
                       full.names = TRUE))
fnRs <- sort(list.files(config$data_path,
                       pattern=config$pattern_reverse,
                       full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

print(fnFs)
# Plot quality profiles: one page per sample, F and R side-by-side
qual_pdf <- file.path(config$figure_output_dir, "quality_profiles_per_sample.pdf")
pdf(qual_pdf, width = 10, height = 5)

for (i in seq_along(sample.names)) {
  pF <- plotQualityProfile(fnFs[i]) +
    ggtitle(paste(sample.names[i], "- Forward"))
  pR <- plotQualityProfile(fnRs[i]) +
    ggtitle(paste(sample.names[i], "- Reverse"))

  print(ggarrange(pF, pR, ncol = 2))
}

dev.off()
cat("Saved per-sample quality profiles to:", qual_pdf, "\n")

# text summary of quality scores by bp
fastq_files <- c(fnFs, fnRs)
qtable <- summarize_quality_by_position(fastq_files, output_csv = file.path(config$data_output_dir, "read_quality_summary.csv"))


# Filter and trim
filtFs <- file.path(config$data_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(config$data_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print("Filtering and trimming reads")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    truncLen=c(as.numeric(config$truncLen_f), as.numeric(config$truncLen_r)),
                    maxN=as.numeric(config$maxN),
                    maxEE=c(as.numeric(config$maxEE), as.numeric(config$maxEE)),
                    truncQ=as.numeric(config$truncQ),
                    rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)
write.csv(out, file.path(config$data_output_dir, "filterAndTrim_summary.csv"))

# Learn error rates
# nbases is capped at 2e8: error learning does not require the full dataset,
# DADA2's own big-data workflow notes ~1e6-5e6 reads is generally sufficient,
# so this just bounds the cost without sacrificing model quality.
print("Learning error rates")
errF <- learnErrors(filtFs, nbases=2e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=2e8, multithread=TRUE)

# Save error plots
pdf(file.path(config$figure_output_dir, "error_plots.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Denoise + merge, streamed one sample at a time.
# This keeps memory flat regardless of sample count: each iteration loads only
# that sample's dereplicated reads, runs dada() on it, merges, then discards
# the heavy intermediate objects (derep + dada-class) before moving on.
# This produces identical results to passing all samples to derepFastq/dada at
# once, since dada()'s default pool=FALSE already processes samples
# independently -- the loop changes memory footprint, not the algorithm.
print("Denoising and merging (streamed per-sample)")

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

denoisedF_n <- numeric(length(sample.names))
denoisedR_n <- numeric(length(sample.names))
names(denoisedF_n) <- sample.names
names(denoisedR_n) <- sample.names

for (sam in sample.names) {
  cat("Processing sample:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  denoisedF_n[sam] <- sum(getUniques(ddF))

  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  denoisedR_n[sam] <- sum(getUniques(ddR))

  mergers[[sam]] <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)

  rm(derepF, derepR, ddF, ddR)
}

# Make sequence table and remove chimeras
# method="consensus" checks each sample independently then takes a majority
# vote across samples -- this scales better than "pooled" (which pools
# abundances across all samples before checking) and is less sensitive to any
# single noisy/high-chimera sample skewing the result.
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Calculate chimera removal statistics
chimera_stats <- sum(seqtab.nochim)/sum(seqtab)
write.table(chimera_stats,
            file=file.path(config$data_output_dir, "chimera_stats.txt"))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
              denoisedF_n[sample.names],
              denoisedR_n[sample.names],
              sapply(mergers, getN),
              rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- cbind(track,
               pct_retained = round(100 * track[, "nonchim"] / track[, "input"], 2))

# Save tracking statistics
write.table(track,
            file=file.path(config$data_output_dir, "dada2_stats.tsv"),
            sep="\t", quote=FALSE)

# Save sequence table
saveRDS(seqtab.nochim,
        file=file.path(config$data_output_dir, "seqtab_final.rds"))
# write.table(seqtab.nochim,
#         file=file.path(config$data_output_dir, "seqtab_final.tsv"))
# Replace long sequences with short ASV IDs before writing
# (you already load taxa for assignTaxonomy downstream, so do this after that step)
asv_result <- assign_asv_ids(seqtab.nochim, taxa, tax_col = "label")

# Write count table with short IDs as column names
write.table(asv_result$seqtab,
            file = file.path(config$data_output_dir, "seqtab_final.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Write the ID -> sequence -> taxonomy mapping separately
write.table(asv_result$mapping,
            file = file.path(config$data_output_dir, "asv_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


cat("Processing complete! Output files have been saved to:",
    config$data_output_dir, "and", config$figure_output_dir, "\n")


# Assign taxonomy
print("Assigning taxonomy")
taxa <- assignTaxonomy(seqtab.nochim, config$taxonomy_db, multithread=TRUE)
taxa <- addSpecies(taxa, config$taxonomy_db_species)
taxa <- data.frame(taxa)
taxa$label <- paste(taxa$Family, taxa$Genus, taxa$Species, sep = "_")
write.table(taxa, file = paste0(config$data_output_dir, "/taxa.csv"))

# evaluate taxa assignment rate
rate_table <- classification_rate_table(seqtab.nochim, taxa)
write.table(rate_table,
            file=file.path(config$data_output_dir, "taxa_assignment_rates.tsv"),
            sep="\t", quote=FALSE)

# evaluate percent composition of at genus level
genus_pct <- taxa_percent_table(seqtab.nochim, taxa, rank = "Genus")
family_pct <- taxa_percent_table(seqtab.nochim, taxa, rank = "Family")
write.table(genus_pct,
            file=file.path(config$data_output_dir, "genus_percent_table.tsv"),
            sep="\t", quote=FALSE)
write.table(family_pct,
            file=file.path(config$data_output_dir, "family_percent_table.tsv"),
            sep="\t", quote=FALSE)

# # make a barplot of composition
# # filter for prevalence
# genus_pct_filtered <- filter_taxa_by_prevalence(genus_pct, min_abundance=1, min_prevalence=0.1)
# # remove control samples (some are all zero and interfere with hclust)
# genus_pct_filtered <- genus_pct_filtered[grepl("^CP101", rownames(genus_pct_filtered)), ]
# gg <- plot_taxa_heatmap(genus_pct_filtered, tax_rank = "Genus", scale_rows=FALSE, cluster_rows=TRUE, cluster_cols=FALSE)
# ggsave(file.path(config$figure_output_dir, "genus_filtered_heatmap.pdf"), gg, width = 16, height = 6)

# # filter for prevalence
# family_pct_filtered <- filter_taxa_by_prevalence(family_pct, min_abundance=1, min_prevalence=0.1)
# # remove control samples (some are all zero and interfere with hclust)
# family_pct_filtered <- family_pct_filtered[grepl("^CP101", rownames(family_pct_filtered)), ]
# gg <- plot_taxa_heatmap(family_pct_filtered, tax_rank = "Family", scale_rows=FALSE, cluster_rows=TRUE, cluster_cols=FALSE)
# ggsave(file.path(config$figure_output_dir, "family_filtered_heatmap.pdf"), gg, width = 16, height = 6)