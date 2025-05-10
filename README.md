# F. prausnitzii Treatment RNA-Seq Analysis

## Project Overview
This repository contains the analysis of Bulk RNA-Seq data from alpha-synuclein overexpressing (ASO) mice, their wild-type (WT) counterparts, and ASO mice treated with F. prausnitzii. The analysis aims to understand the transcriptional changes associated with F. prausnitzii treatment in the context of alpha-synuclein overexpression.

## Repository Structure

### Data Organization
- `data/input/`: Contains raw input data files
  - Raw count matrices
  - Sample metadata
- `data/interim/`: Intermediate processed data files
- `data/results/`: Final processed data tables and results

### Analysis Code
- `notebook/bulk-rnaseq-analysis.qmd`: Primary Quarto notebook containing the complete analysis workflow
- `notebook/R_scripts/`: Additional R scripts for specific analyses
  - `quick_enrich.R`: Script for quick enrichment analysis

### Results and Visualizations
- `figures/`: Contains all generated figures
  - `DESeq2_VolcanoPlots/`: Differential expression volcano plots
  - `GO_enrichment/`: Gene Ontology enrichment analysis plots
  - PCA plots and other visualizations

## Analysis Overview
The analysis pipeline includes:
1. Quality control and preprocessing of raw RNA-Seq data
2. Differential expression analysis using DESeq2
3. Gene Ontology enrichment analysis
4. Comparative analysis between treatment conditions
5. Visualization of results through various plots

## Getting Started
To reproduce the analysis:
1. Place raw data files in `data/input/`
2. Run the main analysis notebook: `notebook/bulk-rnaseq-analysis.qmd`
3. Generated results will be saved in `data/results/` and figures in `figures/`
