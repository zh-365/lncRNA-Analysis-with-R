library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(maftools)
library(dplyr)

## Correlation between MALAT-1 and the Target Gene
# Set the project to "TCGA-OV" for ovarian cancer
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

# Download and prepare the data
GDCdownload(query)
ov_data <- GDCprepare(query)

# Extract expression matrix
count_data <- assay(ov_data)

# Strip version numbers from Ensembl IDs
rownames(count_data) <- sub("\\.\\d+$", "", rownames(count_data))

# Add gene symbols as a column to your data
count_data <- as.data.frame(count_data)
count_data$ensembl_gene_id <- rownames(count_data)
count_data <- merge(count_data, gene_map, by = "ensembl_gene_id")
rownames(count_data) <- count_data$hgnc_symbol

# Filter for MALAT1(ENSG00000251562) and Target Gene
genes_of_interest <- c("ENSG00000251562", "**ENSEMBL ID of target gene**")
filtered_data <- count_data[rownames(count_data) %in% genes_of_interest, ]

# Convert data to a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = filtered_data,
                              colData = colData(ov_data),
                              design = ~ 1)

# Normalize the data
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Extract normalized counts for MALAT1 and Target Gene
malat1_expr <- normalized_counts["ENSG00000251562", ]
target_gene_expr <- normalized_counts["**ENSEMBL ID of target gene**", ]

# Perform correlation test
cor_test <- cor.test(malat1_expr, target_gene_expr, method = "pearson") # Or spearman

# Output correlation coefficient and p-value
cat("Correlation coefficient:", cor_test$estimate, "\n")
cat("P-value:", cor_test$p.value, "\n")

## Data Visualization for Overall Correlation between MALAT-1 and Target Gene
#Prepare data for plotting
plot_data <- data.frame(MALAT1 = malat1_expr, target_gene = target_gene_expr)

# Plot using ggplot2
ggplot(plot_data, aes(x = MALAT1, y = Target_Gene)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Correlation between MALAT-1 and Gene Y in Ovarian Cancer",
       x = "MALAT1 Expression (Normalized)",
       y = "Gene Y Expression (Normalized)") +
  theme_minimal()

## Correlation between MALAT-1 and the Target Gene with miRNAs
# Query and download TCGA RNA-seq data for ovarian cancer (miRNA and gene expression)
query_miRNA <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification"
)

GDCdownload(query_miRNA)
miRNA_data <- GDCprepare(query_miRNA)

query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
ov_data <- GDCprepare(query)

# Prepare miRNA expression matrix
rownames(miRNA_data) <- miRNA_data[, 1]
miRNA_data <- miRNA_data[, -1]

columns_to_remove <- grepl("cross-mapped_", colnames(miRNA_data))
miRNA_data <- miRNA_data[, !columns_to_remove]

columns_to_remove2 <- grepl("reads_per_million", colnames(miRNA_data))
miRNA_data <- miRNA_data[, !columns_to_remove2]

colnames(miRNA_data) <- gsub("read_count_", "", colnames(miRNA_data))

# Prepare gene expression matrix
rownames(gene_matrix) <- make.unique(gene_matrix[, 1])
gene_matrix <- gene_matrix[, -1]

# Match sample IDs between the two datasets
common_samples <- intersect(colnames(miRNA_data), colnames(gene_matrix))
miRNA_matrix <- miRNA_data[, common_samples]
gene_matrix <- gene_matrix[, common_samples]

# Extract Target Gene and MALAT1 expressions
target_gene_expr <- as.numeric(gene_matrix["**ENSEMBL ID of target gene**", ])
MALAT1_expr <- as.numeric(gene_matrix["ENSG00000251562", ])

# Compute correlations
compute_correlation <- function(miRNA_expr, gene_expr) {
  cor.test(miRNA_expr, gene_expr, method = "spearman")
}

# Correlations with Target Gene
correlations <- apply(miRNA_matrix, 1, function(x) compute_correlation(x, target_gene_expr))
results_protein <- data.frame(
  miRNA = rownames(miRNA_matrix),
  correlation = sapply(correlations, function(res) res$estimate),
  p_value = sapply(correlations, function(res) res$p.value)
)

# Correlations with MALAT1
correlations_MALAT1 <- apply(miRNA_matrix, 1, function(x) compute_correlation(x, MALAT1_expr))
results_MALAT1 <- data.frame(
  miRNA = rownames(miRNA_matrix),
  correlation = sapply(correlations_MALAT1, function(res) res$estimate),
  p_value = sapply(correlations_MALAT1, function(res) res$p.value)
)

# Set thresholds
cor_threshold <- 0.3
p_threshold <- 0.001

# Filter highly correlated miRNAs
high_target_gene <- results_target_gene %>%
  filter(abs(correlation) > cor_threshold & p_value < p_threshold)

high_MALAT1 <- results_MALAT1 %>%
  filter(abs(correlation) > cor_threshold & p_value < p_threshold)

# Find overlapping miRNAs
overlap_miRNAs <- intersect(high_target_gene$miRNA, high_MALAT1$miRNA)
overlap_data <- data.frame(
  miRNA = overlap_miRNAs,
  correlation_target_gene = high_target_gene$correlation[match(overlap_miRNAs, high_target_gene$miRNA)],
  correlation_MALAT1 = high_MALAT1$correlation[match(overlap_miRNAs, high_MALAT1$miRNA)]
)

# Scatter plot of correlations
ggplot(overlap_data, aes(x = correlation_target_gene, y = correlation_MALAT1)) +
  geom_point(color = "blue", size = 3) +
  geom_text(aes(label = miRNA), vjust = -0.5, size = 3, color = "red") +
  theme_minimal() +
  labs(
    title = "Overlapping miRNAs Correlated with Target Gene and MALAT1",
    x = "Correlation with Target Gene",
    y = "Correlation with MALAT1"
  )


