# MALAT1 Correlation with a Gene Target in Ovarian Cancer

This repository contains an R script developed as part of a research project the author did under the **Taiwan Experience Education Program (TEEP)**. The objective of this project is to identify and analyze potential interactions or correlations between the long non-coding RNA **MALAT1** and an **undisclosed gene** in **ovarian cancer**, using publicly available data from **The Cancer Genome Atlas (TCGA)**.

## 📌 Project Overview

- **Objective**: Investigate the expression and potential correlation between MALAT1 and a confidential target gene in ovarian cancer.
- **Dataset**: TCGA RNA-seq data (STAR-Counts format).
- **Programming Language**: R
- **Tools & Packages**: `TCGAbiolinks`, `DESeq2`, `ggplot2`, `dplyr`, and base R functions.

## 🧪 Workflow Summary

1. **Data Retrieval**
   - Download TCGA-OV RNA-seq data using `TCGAbiolinks`.

2. **Preprocessing**
   - Normalize and filter the expression data using `DESeq2`.

3. **Analysis**
   - Extract expression levels of MALAT1, target gene, and miRNAs.
   - Perform correlation analysis.
   - Visualize expression patterns and correlation results.

4. **Interpretation**
   - Explore potential interaction based on co-expression data.

> ⚠️ **Note**: The identity of the gene is currently undisclosed due to research confidentiality agreements.

## 📁 Repository Contents

- `malat-1-correlation-analysis.R`: Main R script containing the full analysis pipeline.
- `README.md`: Description of the project and usage instructions.

## 🚀 Getting Started

### Prerequisites

Install the necessary R packages before running the script:

```R
install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2"))
install.packages("ggplot2")
install.packages("dplyr")
