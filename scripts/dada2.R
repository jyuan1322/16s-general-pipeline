#!/usr/bin/env Rscript
# .libPaths("/data/bwh-comppath-seq/software/R-4.1.0//library/")
# run micromamba activate r-microbiome

# NOTE: compare to:
# https://github.com/gerberlab/cdiff_paper_analyses/blob/master/scripts/dada2/dada2.R

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
# Plot quality profiles
p1 <- plotQualityProfile(fnFs[1:4])
ggsave(file.path(config$figure_output_dir, "forward_read_quality.png"), p1)
p2 <- plotQualityProfile(fnRs[1:4])
ggsave(file.path(config$figure_output_dir, "reverse_read_quality.png"), p2)

# text summary of quality scores by bp
fastq_files <- c(fnFs, fnRs)
qtable <- summarize_quality_by_position(fastq_files, output_csv = file.path(config$data_output_dir, "read_quality_summary.csv"))


# Filter and trim
filtFs <- file.path(config$data_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(config$data_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

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
print("Learning error rates")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Save error plots
pdf(file.path(config$figure_output_dir, "error_plots.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Denoise
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Make sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)

# Calculate chimera removal statistics
chimera_stats <- sum(seqtab.nochim)/sum(seqtab)
write.table(chimera_stats,
            file=file.path(config$data_output_dir, "chimera_stats.txt"))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
              sapply(dadaFs, getN),
              sapply(dadaRs, getN),
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
write.table(seqtab.nochim,
        file=file.path(config$data_output_dir, "seqtab_final.tsv"))

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
