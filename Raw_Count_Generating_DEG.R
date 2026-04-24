library(DESeq2)
library(readxl)
library(writexl)
library(dplyr)

# Read raw count file
counts <- read.delim("gene_count.xls", row.names = 1)

# Extract gene id explicitly
counts$gene_id <- rownames(counts)

# Extract no-condition lists
expr <- counts %>%
  select(WT21, WT31, Het11, Het21, Het31) %>%
  as.matrix()

rownames(expr) <- counts$gene_id

# Construct group information
phenotypes <- c("WT", "WT", "Het", "Het", "Het")
colData <- data.frame(phenotypes = factor(phenotypes))
rownames(colData) <- colnames(expr)

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = expr,
  colData = colData,
  design = ~ phenotypes
)

# run DESeq2
dds <- DESeq(dds)

# Prepare normalized counts for Heatmap
norm_counts <- counts(dds, normalized=TRUE)

# Extract DESeq results of GA vs GG
res <- results(dds, contrast = c("phenotypes", "Het", "WT"))

# Add back gene id
res$gene_id <- rownames(res)

# Add back gene_names (symbol)
res_df <- as.data.frame(res) %>%
  left_join(counts %>% select(gene_id, gene_name), by = "gene_id")

# Save the result
write_xlsx(as.data.frame(res_df), "Rodent_DEG_Het_vs_WT.xlsx")
write.csv(norm_counts, "normalized_counts.csv")
