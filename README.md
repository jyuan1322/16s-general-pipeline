# 16S Amplicon Sequencing Pipeline

R pipeline for processing 16S rRNA amplicon sequencing data: quality filtering, denoising, and taxonomic assignment (DADA2), followed by downstream community analysis — composition heatmaps, alpha/beta diversity, and differential abundance testing. The pipeline is split into fixed processing code and per-study parameters, so it can be reused across sequencing runs and experiments without editing the scripts themselves.

## Pipeline overview

1. **`dada2.R`** (run first) — reads raw paired-end fastq files, runs quality filtering/trimming, learns and applies error models, denoises and merges reads, removes chimeras (DADA2), and assigns taxonomy against a reference database (e.g. SILVA). Outputs a sequence table, taxonomy table, and per-sample read-tracking statistics.
2. **`dada2_downstream_analysis.R`** (depends on step 1) — filters ASVs by prevalence and relative abundance, produces a family-level composition heatmap, runs alpha diversity (Shannon) and beta diversity (Bray-Curtis PCoA, PERMANOVA) analyses, and runs DESeq2-based differential abundance testing across a study-defined set of pairwise comparisons.
3. **`dada2_additional_figures.R`** (depends on step 1, independent of step 2) — produces supplementary QC plots: read counts retained at each processing step, and mean quality score across read cycles.

`utils.R` contains shared helper functions, including `read_config()`, and is sourced by all three scripts above.

## Configuration

Settings are split into two files so the same scripts can be reused across studies:

- **`config.ini`** — paths and DADA2 parameters that vary by sequencing run (raw data location, output directories, taxonomy database, quality/truncation thresholds). Also points to a `study_params` file.
- **`study_params.R`** — per-study analysis parameters: sample display order, treatment group labels, pairwise comparisons for differential abundance, samples to exclude, plot colors, and QC thresholds.

## Installation

Create and load micromamba environment for running R:
```
micromamba create -n r-microbiome -f environment.yml
micromamba activate r-microbiome
```

Inside R, install additional packages through BiocManager:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("dada2", "phyloseq", "ComplexHeatmap"))
```


### Setup

1. Copy `config_template.ini` to `config.ini` (or `config_<study>.ini`) and fill in real paths.
2. Copy `study_params_template.R` to `study_params_<study>.R` and fill in your study's sample order, treatment groups, comparisons, and colors.
3. Point to your study params file from the `study_params` field in your config.
4. Keep filled-in config and study params files out of version control if they contain sensitive or machine-specific paths (add them to `.gitignore`); only the templates need to be tracked.

## Usage

Run from the directory containing `utils.R` (scripts source it with a relative path):

```bash
Rscript dada2.R config.ini
Rscript dada2_downstream_analysis.R config.ini
Rscript dada2_additional_figures.R config.ini
```

## Repository structure

```
.
├── dada2.R                  # Step 1: DADA2 filtering, denoising, taxonomy assignment
├── dada2_downstream_analysis.R    # Step 2: composition heatmap, diversity, differential abundance
├── dada2_additional_figures.R     # Step 3: supplementary QC plots
├── utils.R                         # Shared helper functions, incl. read_config()
├── config_template.ini             # Path/parameter template (copy and fill in)
├── study_params_template.R         # Per-study parameter template (copy and fill in)
├── database/                       # Taxonomy reference database (e.g. SILVA)
├── metadata/                       # Sample metadata CSV(s)
├── data/
│   └── raw_fastq/                  # Raw paired-end fastq files
├── data_output/                    # Generated tables, stats, and CSVs
└── figure_output/                  # Generated plots
```

## Dependencies

R packages: `dada2`, `phyloseq`, `DESeq2`, `apeglm`, `vegan`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `reshape2`, `pheatmap`, `ggpubr`, `ggrepel`, `scales`, `forcats`, `rstatix`, `viridis`, `ShortRead`