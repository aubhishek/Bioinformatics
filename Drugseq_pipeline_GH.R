# Remove all variables
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
###Alternative
#this_file <- commandArgs(trailingOnly = FALSE)
#script_path <- normalizePath(sub("--file=", "", this_file[grep("--file=", this_file)]))
#setwd(dirname(script_path))

library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(Rtsne)
library(umap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(msigdbr)
library(enrichplot)
library(fgsea)
library(reshape2)


counts_df_gene_name <- read.csv("merged_counts_df_gene_name_filtered.csv")
counts_df_gene_name <- counts_df_gene_name[, -1]
counts_df_gene_name$gene_name <- as.character(counts_df_gene_name$gene_name)
counts_df_gene_name$gene_name <- make.unique(counts_df_gene_name$gene_name)
rownames(counts_df_gene_name) <- counts_df_gene_name$gene_name
#counts_df_gene_name <- counts_df_gene_name[, -1]
#counts_df_gene_name <- counts_df_gene_name[, !colnames(counts_df_gene_name) %in% "gene_name", drop = FALSE]
View(counts_df_gene_name)




#################################
##DESEQ2 Analysis

# Ensure merged_counts is loaded
# If not already loaded, uncomment:
# merged_counts <- fread("merged_counts_gene_names.csv", header = TRUE, row.names = 1)
# Prepare count matrix for DESeq2 (remove gene_name column)
counts_df_gene_name_wo_gn <- counts_df_gene_name[, !colnames(counts_df_gene_name) %in% c("gene_name"), drop = FALSE]
merged_counts <- counts_df_gene_name_wo_gn

# Extract column names (no gene_name column in merged_counts)
sample_cols <- colnames(merged_counts)
cat("First few sample_cols:\n")
print(head(sample_cols))  # Debug: Confirm column names

# Filter for Cell_line_Genotype_1 samples only (exclude Cell_line_1A4_Genotype_1 and Empty)
Genotype_1_cols <- grep("Samples[1-4]_[A-Z][0-9]{2}__Cell_line_Genotype_1_", sample_cols, value = TRUE)
Genotype_1_cols <- Genotype_1_cols[!grepl("Empty", Genotype_1_cols)]  # Exclude Empty wells if present
cat("First few Genotype_1_cols:\n")
print(head(Genotype_1_cols))  # Debug: Confirm Genotype_1 columns

# Extract conditions by removing SamplesX_[A-Z][0-9]{2}__Cell_line_Genotype_1_
conditions <- sub("Samples[1-4]_[A-Z][0-9]{2}__Cell_line_Genotype_1_", "", Genotype_1_cols)
cat("After first sub (remove prefix):\n")
print(head(conditions))  # Debug: Should see "1_24hr", "Compound-1_100nM_1_24hr", etc.

# Assign "NT" to specific no-treatment conditions
nt_conditions <- c("1_24hr", "2_24hr", "3_24hr", "4_24hr", "5_24hr", "6_24hr")
conditions <- ifelse(conditions %in% nt_conditions, "NT", conditions)
cat("After NT assignment:\n")
print(head(conditions))  # Debug: Should see "NT" for no-treatment samples

# Remove replicate and time suffixes from remaining conditions
conditions <- sub("_[1-6]_24hr$", "", conditions)
cat("After removing suffixes:\n")
print(head(conditions))  # Debug: Should see "NT", "DMSO", "Compound-1_100nM", etc.

# Get unique conditions, excluding "Empty" if present
unique_conditions <- unique(conditions)
unique_conditions <- unique_conditions[!unique_conditions %in% c("Empty")]
cat("Unique conditions:\n")
print(unique_conditions)  # Should show "NT", "DMSO", "Compound-1_100nM", "Compound-2_100nM", etc.

# List to store dataframes for Cell_line_Genotype_1 only
condition_dfs <- list()
prefix <- "Samples[1-4]_[A-Z][0-9]{2}__Cell_line_Genotype_1_"
for (cond in unique_conditions) {
  # Find columns matching this condition in Cell_line_Genotype_1 only
  if (cond == "NT") {
    # For NT, match columns with just a number before _24hr
    cond_cols <- grep(paste0(prefix, "[1-6]_24hr$"), Genotype_1_cols, value = TRUE)
  } else {
    # For other conditions, match the specific treatment
    cond_cols <- grep(paste0(prefix, cond, "_[1-6]_24hr$"), Genotype_1_cols, value = TRUE)
  }
  if (length(cond_cols) > 0) {
    # Create a data.frame with gene_name column
    cond_df <- data.frame(
      gene_name = rownames(merged_counts),
      merged_counts[, cond_cols]
    )
    condition_dfs[[cond]] <- cond_df
    cat("Created dataframe for", cond, "with", length(cond_cols), "replicates\n")
  } else {
    cat("No columns found for", cond, "\n")
  }
}

# Assign dataframes to global environment
list2env(condition_dfs, envir = .GlobalEnv)

# save the dataframes
for (cond in unique_conditions) {
  if (cond %in% names(condition_dfs)) {
    write.csv(condition_dfs[[cond]], file = paste0(cond, "_dataframe.csv"), row.names = FALSE)
  } else {
    cat("Warning: No dataframe for", cond, "- Skipping save\n")
  }
}


# Check condition dataframes for confirmation
head(DMSO)  # Should have gene_name + DMSO replicate columns
head(NT)    # Should have gene_name + NT replicate columns (e.g., Samples1_A01__Cell_line_Genotype_1_1_24hr, Samples2_A01__Cell_line_Genotype_1_1_24hr)

# Function to perform pairwise DESeq2 analysis
run_deseq_pairwise <- function(count_data, col_data, treatment_condition, control_condition = "DMSO", output_dir) {
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = col_data, 
                                design = ~ condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results (treatment vs control)
  res <- results(dds, contrast = c("condition", treatment_condition, control_condition))
  
  # Save to CSV
  write.csv(as.data.frame(res), 
            file = file.path(output_dir, paste0(treatment_condition, "_vs_", control_condition, "_deseq_results.csv")),
            row.names = TRUE)
  cat("Completed DESeq2 for", treatment_condition, "vs", control_condition, "\n")
  
  return(as.data.frame(res))
}

# Create output directory
output_dir <- "DESeq_Results_Cell_line_Genotype_1_DMSO"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Prepare and run DESeq2 for each condition vs DMSO in Cell_line_Genotype_1
deseq_results_list <- list()

for (cond in unique_conditions) {
  if (cond == "DMSO") next
  
  # Get control and treatment dataframes
  control_df <- condition_dfs[["DMSO"]]
  treatment_df <- condition_dfs[[cond]]
  
  # Handle duplicate gene names by renaming (e.g., gene.1, gene.2)
  control_df$gene_name <- make.unique(control_df$gene_name)
  treatment_df$gene_name <- make.unique(treatment_df$gene_name)
  
  # Merge control and treatment
  count_data_df <- merge(control_df, treatment_df, by = "gene_name")
  
  # Convert to matrix with proper row names
  gene_names <- count_data_df$gene_name
  count_data <- as.matrix(count_data_df[, -1])  # Remove gene_name
  rownames(count_data) <- gene_names
  
  # Create colData with exact column names
  n_control <- ncol(control_df) - 1  # Number of DMSO replicates
  n_treatment <- ncol(treatment_df) - 1  # Number of treatment replicates
  col_data <- data.frame(
    condition = factor(c(rep("DMSO", n_control), rep(cond, n_treatment)),
                       levels = c("DMSO", cond)),
    row.names = colnames(count_data)
  )
  
  # Run DESeq2
  result <- run_deseq_pairwise(count_data, col_data, treatment_condition = cond, output_dir = output_dir)
  deseq_results_list[[cond]] <- result
}

# Save results list for later use
saveRDS(deseq_results_list, file = file.path(output_dir, "deseq_results_list.rds"))
cat("Saved DESeq2 results list to", file = file.path(output_dir, "deseq_results_list.rds"), "\n")





##############################
##Volcano plot

# Create a folder for volcano plots
output_dir <- "Volcano_Plots_Cell_line_Genotype_1_DMSO"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define comparison titles for Cell_line_Genotype_1 conditions
comparison_titles <- list(
  "NT" = "No Treatment vs DMSO",
  "Compound-4_100nM" = "Compound-4 (100nM) vs DMSO",
  "Compound-6_1uM" = "Compound-6 (1uM) vs DMSO",
  "Compound-5_1uM" = "Compound-5 (1uM) vs DMSO",
  "Compound-7_100nM" = "Compound-7 (100nM) vs DMSO",
  "Compound-3_100nM" = "Compound-3 (100nM) vs DMSO",
  "Compound-1_100nM" = "Compound-1 (100nM) vs DMSO",
  "Compound-2_100nM" = "Compound-2 (100nM) vs DMSO"
)

# Get list of DESeq2 result files from Cell_line_Genotype_1_DMSO
result_dir <- "DESeq_Results_Cell_line_Genotype_1_DMSO"
deseq_result_files <- list.files(result_dir, pattern = "_vs_DMSO_deseq_results.csv", full.names = TRUE)
conditions <- gsub("_vs_DMSO_deseq_results.csv", "", basename(deseq_result_files))

# Verify all expected conditions have result files
expected_conditions <- c("NT", "Compound-4_100nM", "Compound-6_1uM", "Compound-5_1uM", 
                         "Compound-7_100nM", "Compound-3_100nM", "Compound-1_100nM", "Compound-2_100nM")
missing_conditions <- setdiff(expected_conditions, conditions)
if (length(missing_conditions) > 0) {
  cat("Warning: Missing result files for conditions:", paste(missing_conditions, collapse = ", "), "\n")
}

# Loop over each DESeq2 result file
for (i in seq_along(conditions)) {
  condition <- conditions[i]
  file_path <- deseq_result_files[i]
  
  # Load and filter DESeq2 results
  deseq_res <- read.csv(file_path, row.names = 1)
  deseq_res <- deseq_res[!is.na(deseq_res$log2FoldChange) & !is.na(deseq_res$padj), ]
  
  # Create data frame for plotting
  plot_data <- data.frame(
    log2FoldChange = deseq_res$log2FoldChange,
    negLog10Pval = -log10(deseq_res$padj),
    padj = deseq_res$padj,
    row.names = rownames(deseq_res)
  )
  
  # Save the plot_data as CSV
  write.csv(plot_data, 
            file = file.path(output_dir, paste0(condition, "_vs_DMSO_plot_data.csv")), 
            row.names = TRUE)
  
  # Filter out Inf and extreme values
  plot_data <- plot_data[!is.infinite(plot_data$negLog10Pval), ]  # Remove Inf
  plot_data <- plot_data[plot_data$negLog10Pval <= 310, ]  # Cap at 310
  
  # Add color column based on thresholds
  plot_data$color <- "grey"
  plot_data$color[plot_data$log2FoldChange >= 1 & plot_data$padj < 0.05] <- "red"  # Upregulated
  plot_data$color[plot_data$log2FoldChange <= -1 & plot_data$padj < 0.05] <- "green"  # Downregulated
  
  # Generate volcano plot
  volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = negLog10Pval, color = color)) +
    geom_point(size = 2) + scale_color_identity() + theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(), 
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 30)
    ) +
    labs(title = paste("Volcano Plot:", comparison_titles[[condition]]), x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
    xlim(-10, 10) + ylim(0, 100)
  
  # Save the plot
  ggsave(filename = file.path(output_dir, paste0(condition, "_vs_DMSO_volcano_plot.png")),
         plot = volcano_plot, width = 10, height = 8, bg = "white")
  
  # Print progress
  cat("Generated volcano plot for", condition, "vs DMSO - Saved to", output_dir, "\n")
}



########################
#Dimensional reduction
library(Rtsne)
library(umap)


# Create output directory
output_dir <- "Dim_Reduction_Plots_Cell_line_Genotype_1_DMSO"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Step 1: Prepare count data and colData
# Ensure merged_counts is loaded
# If not loaded, uncomment:
# merged_counts <- fread("merged_counts_gene_names.csv", header = TRUE, row.names = 1)

# Create count_data with gene_name column
count_data <- data.frame(
  gene_name = rownames(merged_counts),
  merged_counts
)

# Rename duplicate gene names
count_data$gene_name <- make.unique(count_data$gene_name)

# Set row names and remove gene_name column
rownames(count_data) <- count_data$gene_name
count_data <- count_data[, -1]

# Extract conditions from column names
sample_cols <- colnames(count_data)
conditions <- sub("Samples[1-4]_[A-Z][0-9]{2}__Cell_line_Genotype_1_", "", sample_cols)
# Assign "NT" to no-treatment samples
nt_conditions <- paste0(1:6, "_24hr")
conditions <- ifelse(conditions %in% nt_conditions, "NT", conditions)
# Remove replicate and time suffixes
conditions <- sub("_[1-6]_24hr$", "", conditions)

# Create colData with DMSO as reference level
col_data <- data.frame(
  condition = factor(conditions, levels = c("DMSO", setdiff(unique(conditions), "DMSO"))),
  row.names = sample_cols
)

# Exclude Empty samples and non-Cell_line_Genotype_1 samples
keep_samples <- col_data$condition != "Empty" & grepl("Cell_line_Genotype_1", rownames(col_data))
count_data <- count_data[, keep_samples]
col_data <- col_data[keep_samples, , drop = FALSE]
cat("Number of samples after filtering:", ncol(count_data), "\n")

# Step 2: Filter low-count genes and normalize
# Keep genes with non-zero counts in at least 10% of samples
min_nonzero_samples <- ceiling(ncol(count_data) * 0.1)
keep_genes <- rowSums(count_data > 0) >= min_nonzero_samples
count_data_filtered <- count_data[keep_genes, ]
cat("Number of genes after non-zero filter:", nrow(count_data_filtered), "\n")

# Additional filter: Ensure some genes have non-zero counts in all samples
keep_genes_nonzero_all <- rowSums(count_data_filtered > 0) == ncol(count_data_filtered)
if (sum(keep_genes_nonzero_all) == 0) {
  cat("Warning: No genes with non-zero counts across all samples. Relaxing filter.\n")
  # Fall back to genes with mean count > 1
  keep_genes <- rowMeans(count_data_filtered) > 1
  count_data_filtered <- count_data_filtered[keep_genes, ]
} else {
  cat("Found", sum(keep_genes_nonzero_all), "genes with non-zero counts across all samples\n")
}

# Check for all-zero genes
any_all_zeros <- any(rowSums(count_data_filtered) == 0)
cat("Any genes with all zeros after filtering:", any_all_zeros, "\n")
if (any_all_zeros) {
  keep_genes <- rowSums(count_data_filtered) > 0
  count_data_filtered <- count_data_filtered[keep_genes, ]
  cat("Removed", sum(!keep_genes), "all-zero genes\n")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data_filtered, 
                              colData = col_data, 
                              design = ~ condition)

# Use vst() with fallback to manual log transformation
tryCatch({
  vsd <- vst(dds, blind = TRUE)
  vsd_mat <- assay(vsd)
  cat("vst transformation successful\n")
}, error = function(e) {
  cat("vst failed:", e$message, "\nFalling back to manual log transformation with pseudo-counts\n")
  sample_sums <- colSums(count_data_filtered)
  cat("Sample sums (first 5):", sample_sums[1:5], "...\n")
  zero_samples <- sum(sample_sums == 0)
  cat("Number of all-zero samples:", zero_samples, "\n")
  if (zero_samples > 0) {
    keep_samples <- sample_sums > 0
    count_data_filtered <<- count_data_filtered[, keep_samples]
    col_data <<- col_data[keep_samples, , drop = FALSE]
    cat("Removed", zero_samples, "all-zero samples\n")
    dds <<- DESeqDataSetFromMatrix(countData = count_data_filtered, colData = col_data, design = ~ condition)
  }
  size_factors <- colSums(count_data_filtered) / median(colSums(count_data_filtered))
  cat("Size factors (first 5):", size_factors[1:5], "...\n")
  norm_counts <- t(t(count_data_filtered) / size_factors)
  vsd_mat <<- log2(norm_counts + 1)
  cat("Manual log transformation completed\n")
})

# Step 3: Filter out zero-variance genes and samples
gene_vars <- apply(vsd_mat, 1, var)
cat("Number of genes before variance filter:", nrow(vsd_mat), "\n")
cat("Number of genes with NA variance:", sum(is.na(gene_vars)), "\n")
cat("Number of genes with zero variance:", sum(gene_vars == 0, na.rm = TRUE), "\n")
vsd_mat_filtered <- vsd_mat[!is.na(gene_vars) & gene_vars > 0, ]
cat("Number of genes after variance filter:", nrow(vsd_mat_filtered), "\n")

# Check sample variance
col_vars <- apply(vsd_mat_filtered, 2, var)
cat("Number of samples with zero variance:", sum(col_vars == 0), "\n")
if (any(col_vars == 0)) {
  keep_samples <- col_vars > 0
  vsd_mat_filtered <<- vsd_mat_filtered[, keep_samples]
  col_data <<- col_data[keep_samples, , drop = FALSE]
  cat("Removed", sum(col_vars == 0), "samples with zero variance\n")
}
cat("Final dimensions of vsd_mat_filtered:", dim(vsd_mat_filtered), "\n")

# Check if matrix is usable
if (nrow(vsd_mat_filtered) < 2) {
  stop("Error: Too few genes with variance after filtering. Cannot perform PCA, t-SNE, or UMAP.")
}
if (ncol(vsd_mat_filtered) < 2) {
  stop("Error: Too few samples with variance after filtering. Cannot perform PCA, t-SNE, or UMAP.")
}

# Step 4: PCA
pca_res <- prcomp(t(vsd_mat_filtered), scale. = TRUE)
cat("PCA results dimensions:", dim(pca_res$x), "\n")
cat("First few PC1 values:", head(pca_res$x[, 1]), "\n")
cat("First few PC2 values:", head(pca_res$x[, 2]), "\n")

# Create PCA data frame
pca_data <- as.data.frame(pca_res$x[, 1:2])
cat("PCA data dimensions:", dim(pca_data), "\n")
cat("Head of PCA data:\n")
print(head(pca_data))

# Ensure condition column matches the samples
cat("Number of samples in col_data:", nrow(col_data), "\n")
pca_data$condition <- col_data$condition
cat("Head of PCA data with condition:\n")
print(head(pca_data))

# Ensure condition is a factor
pca_data$condition <- factor(pca_data$condition, levels = unique(col_data$condition))
cat("Levels of condition in pca_data:", levels(pca_data$condition), "\n")

# Calculate the percentage of variance explained by each principal component
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

# Define color palette
unique_conditions <- unique(col_data$condition)
colors <- rainbow(length(unique_conditions))
names(colors) <- unique_conditions
cat("Color palette:\n")
print(colors)

# Check plot data
cat("Summary of pca_data before plotting:\n")
print(summary(pca_data))

# Generate PCA plot with explicit color mapping
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of Cell_line_Genotype_1 Samples (No Empty)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "right")
print(pca_plot)
ggsave(file.path(output_dir, "pca_Cell_line_Genotype_1_dmso.png"), pca_plot, width = 8, height = 6, bg = "white")
cat("PCA plot saved to", file.path(output_dir, "pca_Cell_line_Genotype_1_dmso.png"), "\n")

# Generate PCA plot without legend
pca_plot_no_legend <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of Cell_line_Genotype_1 Samples (No Empty)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "pca_Cell_line_Genotype_1_dmso_no_legend.png"), pca_plot_no_legend, width = 10, height = 8, bg = "white")
cat("PCA plot (no legend) saved to", file.path(output_dir, "pca_Cell_line_Genotype_1_dmso_no_legend.png"), "\n")


# Create PCA plot with conditions
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors, name = "Condition") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "PCA of KO Plate Samples (No Empty, No Human RNA)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )

# Save the PCA plot
ggsave(
  filename = file.path(output_dir, "pca_with_variance_explained.png"),
  plot = pca_plot, width = 10, height = 8, bg = "white"
)



# Step 5: t-SNE
set.seed(42)
tsne_res <- Rtsne(t(vsd_mat_filtered), dims = 2, perplexity = min(5, (ncol(vsd_mat_filtered) - 1) / 3), max_iter = 1000)
tsne_data <- as.data.frame(tsne_res$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$condition <- factor(col_data$condition, levels = unique(col_data$condition))
cat("t-SNE data dimensions:", dim(tsne_data), "\n")
cat("Head of t-SNE data:\n")
print(head(tsne_data))

# t-SNE plot with legend
tsne_plot_with_legend <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of Cell_line_Genotype_1 Samples (No Empty)", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "right")
ggsave(file.path(output_dir, "tsne_Cell_line_Genotype_1_dmso_with_legend.png"), tsne_plot_with_legend, width = 12, height = 8, bg = "white")
cat("t-SNE plot (with legend) saved to", file.path(output_dir, "tsne_Cell_line_Genotype_1_dmso_with_legend.png"), "\n")

# t-SNE plot without legend
tsne_plot_no_legend <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of Cell_line_Genotype_1 Samples (No Empty)", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "tsne_Cell_line_Genotype_1_dmso_no_legend.png"), tsne_plot_no_legend, width = 10, height = 8, bg = "white")
cat("t-SNE plot (no legend) saved to", file.path(output_dir, "tsne_Cell_line_Genotype_1_dmso_no_legend.png"), "\n")


# Step 6: UMAP
set.seed(42)
umap_res <- umap(t(vsd_mat_filtered), n_components = 2)
umap_data <- as.data.frame(umap_res$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_data$condition <- factor(col_data$condition, levels = unique(col_data$condition))
cat("UMAP data dimensions:", dim(umap_data), "\n")
cat("Head of UMAP data:\n")
print(head(umap_data))

# UMAP plot with legend
umap_plot_with_legend <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "UMAP of Cell_line_Genotype_1 Samples (No Empty)", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "right")
ggsave(file.path(output_dir, "umap_Cell_line_Genotype_1_dmso_with_legend.png"), umap_plot_with_legend, width = 12, height = 8, bg = "white")
cat("UMAP plot (with legend) saved to", file.path(output_dir, "umap_Cell_line_Genotype_1_dmso_with_legend.png"), "\n")

# UMAP plot without legend
umap_plot_no_legend <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "UMAP of Cell_line_Genotype_1 Samples (No Empty)", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "umap_Cell_line_Genotype_1_dmso_no_legend.png"), umap_plot_no_legend, width = 10, height = 8, bg = "white")
cat("UMAP plot (no legend) saved to", file.path(output_dir, "umap_Cell_line_Genotype_1_dmso_no_legend.png"), "\n")

# Print completion
cat("Generated PCA, t-SNE, and UMAP plots - Saved to", output_dir, "\n")




###########################
# Loading DESEQs

library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(data.table)

# Define function to load DESeq2 results
load_deseq_results <- function(file_path, gene_names) {
  res <- read.csv(file_path, row.names = 1)
  cat("Loaded DESeq2 results from", file_path, "- Rows:", nrow(res), "\n")
  
  # Check for gene name consistency
  if (!all(rownames(res) %in% gene_names)) {
    cat("Warning: Some gene names in DESeq2 results do not match provided gene_names.\n")
    cat("Missing genes (first few):", head(setdiff(rownames(res), gene_names)), "\n")
  }
  
  # Allow for duplicate gene names (e.g., gene.1, gene.2)
  if (nrow(res) < length(gene_names)) {
    cat("Note: DESeq2 results have fewer rows (", nrow(res), 
        ") than gene names (", length(gene_names), ") due to filtering or duplicates.\n")
  }
  
  return(res)
}

# Define directory and load results
result_dir <- "DESeq_Results_Cell_line_Genotype_1_DMSO"
deseq_result_files <- list.files(result_dir, pattern = "_vs_DMSO_deseq_results.csv", full.names = TRUE)
cat("Found", length(deseq_result_files), "DESeq2 result files in", result_dir, "\n")

# Verify expected conditions
expected_conditions <- c("NT", "Compound-4_100nM", "Compound-6_1uM", "Compound-5_1uM", 
                         "Compound-7_100nM", "Compound-3_100nM", "Compound-1_100nM", "Compound-2_100nM")
conditions <- gsub("_vs_DMSO_deseq_results.csv", "", basename(deseq_result_files))
missing_conditions <- setdiff(expected_conditions, conditions)
if (length(missing_conditions) > 0) {
  cat("Warning: Missing result files for conditions:", paste(missing_conditions, collapse = ", "), "\n")
}

# Get gene names from merged_counts
# Ensure merged_counts is loaded
# If not loaded, uncomment:
# merged_counts <- fread("merged_counts_gene_names.csv", header = TRUE, row.names = 1)
gene_names <- rownames(merged_counts)
cat("Number of gene names:", length(gene_names), "\n")

# Load DESeq2 results with gene names
deseq_results <- lapply(deseq_result_files, load_deseq_results, gene_names = gene_names)
names(deseq_results) <- conditions
cat("Loaded DESeq2 results for conditions:", names(deseq_results), "\n")










#######################
##Heatmap
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(data.table)

# Define the heatmap function (unchanged)
generate_heatmap <- function(counts_data, col_data, deseq_results_file, comparison_name, cond_name,
                             padj_threshold = 0.05, lfc_threshold = 1, 
                             output_dir = "Heatmaps_Cell_line_Genotype_1_DMSO") {
  
  # Rename duplicate gene names
  counts_data$gene_name <- make.unique(counts_data$gene_name)
  cat("Dimensions of counts_data after renaming duplicates:", dim(counts_data), "\n")
  
  # Set row names and remove gene_name column explicitly
  rownames(counts_data) <- counts_data$gene_name
  counts_data <- counts_data[, colnames(counts_data) != "gene_name"]
  cat("Dimensions of counts_data after removing gene_name:", dim(counts_data), "\n")
  
  # Verify column types
  col_types <- sapply(counts_data, class)
  non_numeric_cols <- names(col_types)[!col_types %in% c("numeric", "integer")]
  if (length(non_numeric_cols) > 0) {
    cat("Error: Non-numeric columns detected:", non_numeric_cols, "\n")
    stop("Counts data contains non-numeric columns. Please check input data.")
  }
  
  # Check for all-zero samples
  sample_sums <- colSums(counts_data)
  cat("Sample sums (first 5):", sample_sums[1:5], "...\n")
  zero_samples <- sum(sample_sums == 0)
  if (zero_samples > 0) {
    keep_samples <- sample_sums > 0
    counts_data <- counts_data[, keep_samples]
    col_data <- col_data[keep_samples, , drop = FALSE]
    cat("Removed", zero_samples, "all-zero samples\n")
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                                colData = col_data, 
                                design = ~ condition)
  
  # Run DESeq2 and apply variance-stabilizing transformation
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  vsd_counts <- assay(vsd)
  cat("Dimensions of VST counts:", dim(vsd_counts), "\n")
  
  # Load DESeq2 results
  deseq_results <- read.csv(deseq_results_file, row.names = 1)
  cat("Dimensions of DESeq2 results:", dim(deseq_results), "\n")
  
  # Convert columns to numeric
  deseq_results$log2FoldChange <- as.numeric(as.character(deseq_results$log2FoldChange))
  deseq_results$padj <- as.numeric(as.character(deseq_results$padj))
  
  # Filter for significant genes
  sig_genes <- deseq_results[
    !is.na(deseq_results$padj) & 
      deseq_results$padj < padj_threshold & 
      abs(deseq_results$log2FoldChange) > lfc_threshold, 
  ]
  sig_gene_names <- rownames(sig_genes)
  cat("Number of significant genes (padj <", padj_threshold, ", |log2FC| >", lfc_threshold, "):", length(sig_gene_names), "\n")
  
  # Ensure sig_gene_names exist in vsd_counts
  sig_gene_names <- sig_gene_names[sig_gene_names %in% rownames(vsd_counts)]
  cat("Number of significant genes after matching with VST counts:", length(sig_gene_names), "\n")
  
  # Check if there are any significant genes
  if (length(sig_gene_names) < 1) {
    cat("No significant genes for", comparison_name, "- Skipping heatmap\n")
    return(NULL)
  }
  
  # Subset VST counts to significant genes
  vsd_sig <- vsd_counts[sig_gene_names, ]
  cat("Dimensions of VST counts for significant genes:", dim(vsd_sig), "\n")
  
  # Check for zero-variance rows in vsd_sig
  row_vars <- apply(vsd_sig, 1, var)
  cat("Number of zero-variance genes in vsd_sig:", sum(row_vars == 0), "\n")
  if (any(row_vars == 0)) {
    keep_genes <- row_vars > 0
    vsd_sig <- vsd_sig[keep_genes, ]
    cat("Removed", sum(row_vars == 0), "zero-variance genes from vsd_sig\n")
    if (nrow(vsd_sig) < 1) {
      cat("After removing zero-variance genes, no genes remain for", comparison_name, "- Skipping heatmap\n")
      return(NULL)
    }
  }
  
  # Check for NA values in vsd_sig
  if (any(is.na(vsd_sig))) {
    cat("Warning: NA values found in vsd_sig - Removing rows with NAs\n")
    vsd_sig <- vsd_sig[complete.cases(vsd_sig), ]
    cat("Dimensions of vsd_sig after removing NAs:", dim(vsd_sig), "\n")
    if (nrow(vsd_sig) < 1) {
      cat("After removing NAs, no genes remain for", comparison_name, "- Skipping heatmap\n")
      return(NULL)
    }
  }
  
  # Annotation for columns
  annotation_col <- data.frame(Condition = col_data$condition)
  rownames(annotation_col) <- colnames(vsd_sig)
  
  # Define color palette
  colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  
  # Generate heatmap
  heatmap_plot <- pheatmap(vsd_sig,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           annotation_col = annotation_col,
                           color = colors,
                           scale = "row",
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           main = paste("Heatmap of Significant DE Genes:", comparison_name),
                           fontsize = 10)
  
  # Create output directory if it doesn’t exist
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Save the plot with a unique filename based on cond_name
  heatmap_filename <- paste0("heatmap_", cond_name, "_vs_DMSO.png")
  png(filename = file.path(output_dir, heatmap_filename),
      width = 10, height = 8, units = "in", res = 300)
  print(heatmap_plot)
  dev.off()
  cat("Saved heatmap to:", file.path(output_dir, heatmap_filename), "\n")
  
  return(heatmap_plot)
}

# Define comparisons (using DMSO as control, updated to use dots)
comparisons <- list(
  NT = c(control = "DMSO", treatment = "NT"),
  Compound_4_100nM = c(control = "DMSO", treatment = "Compound.4_100nM"),
  Compound_6_1uM = c(control = "DMSO", treatment = "Compound.6_1uM"),
  Compound_5_1uM = c(control = "DMSO", treatment = "Compound.5_1uM"),
  Compound_7_100nM = c(control = "DMSO", treatment = "Compound.7_100nM"),
  Compound_3_100nM = c(control = "DMSO", treatment = "Compound.3_100nM"),
  Compound_1_100nM = c(control = "DMSO", treatment = "Compound.1_100nM"),
  Compound_2_100nM = c(control = "DMSO", treatment = "Compound.2_100nM")
)

# Define comparison titles for better labeling (updated to use dots)
comparison_titles <- list(
  "NT" = "NT vs DMSO",
  "Compound.4_100nM" = "Compound.4 (100nM) vs DMSO",
  "Compound.6_1uM" = "Compound.6 (1uM) vs DMSO",
  "Compound.5_1uM" = "Compound.5 (1uM) vs DMSO",
  "Compound.7_100nM" = "Compound.7 (100nM) vs DMSO",
  "Compound.3_100nM" = "Compound.3 (100nM) vs DMSO",
  "Compound.1_100nM" = "Compound.1 (100nM) vs DMSO",
  "Compound.2_100nM" = "Compound.2 (100nM) vs DMSO"
)

# Loop over comparisons to generate heatmaps
output_dir <- "Heatmaps_Cell_line_Genotype_1_DMSO"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Ensure merged_counts is loaded
# If not loaded, uncomment:
# merged_counts <- fread("merged_counts_gene_names.csv", header = TRUE, row.names = 1)

for (cond_name in names(comparisons)) {
  # Get control and treatment dataframes from condition_dfs
  control_df <- condition_dfs[["DMSO"]]
  treatment_df <- condition_dfs[[comparisons[[cond_name]]["treatment"]]]
  
  # Check if dataframes exist and have data
  if (is.null(control_df) || is.null(treatment_df) || ncol(control_df) == 0 || ncol(treatment_df) == 0) {
    cat("Error: Invalid or empty dataframe for", cond_name, "- Skipping\n")
    cat("Control dimensions:", if (is.null(control_df)) "NULL" else dim(control_df), "\n")
    cat("Treatment dimensions:", if (is.null(treatment_df)) "NULL" else dim(treatment_df), "\n")
    next
  }
  
  # Debug: Check column names and dimensions
  cat("Control columns:", colnames(control_df), "\n")
  cat("Treatment columns:", colnames(treatment_df), "\n")
  cat("Control dimensions before selection:", dim(control_df), "\n")
  cat("Treatment dimensions before selection:", dim(treatment_df), "\n")
  
  # Exclude known non-numeric columns
  control_df <- control_df[, !colnames(control_df) %in% c("gene_name", "gene_name.1"), drop = FALSE]
  treatment_df <- treatment_df[, !colnames(treatment_df) %in% c("gene_name", "gene_name.1"), drop = FALSE]
  
  # Check if dataframes have valid sample columns
  if (ncol(control_df) == 0 || ncol(treatment_df) == 0) {
    cat("Error: No valid sample columns for", cond_name, "- Skipping\n")
    cat("Control columns after selection:", colnames(control_df), "\n")
    cat("Treatment columns after selection:", colnames(treatment_df), "\n")
    next
  }
  
  # Debug: Check dimensions after selection
  cat("Control dimensions after selection:", dim(control_df), "\n")
  cat("Treatment dimensions after selection:", dim(treatment_df), "\n")
  
  # Add gene_name column
  control_df <- data.frame(gene_name = rownames(control_df), control_df)
  treatment_df <- data.frame(gene_name = rownames(treatment_df), treatment_df)
  
  # Verify column types before merge
  control_cols <- sapply(control_df[, -1], class)
  treatment_cols <- sapply(treatment_df[, -1], class)
  if (any(!control_cols %in% c("numeric", "integer")) || any(!treatment_cols %in% c("numeric", "integer"))) {
    cat("Error: Non-numeric columns in condition_dfs for", cond_name, "\n")
    cat("Control non-numeric:", names(control_cols)[!control_cols %in% c("numeric", "integer")], "\n")
    cat("Treatment non-numeric:", names(treatment_cols)[!treatment_cols %in% c("numeric", "integer")], "\n")
    next
  }
  
  # Merge count data
  counts_data <- merge(control_df, treatment_df, by = "gene_name")
  cat("Merged counts data for", cond_name, "vs DMSO - Dimensions:", dim(counts_data), "\n")
  
  # Create colData
  n_control <- ncol(control_df) - 1  # Exclude gene_name
  n_treatment <- ncol(treatment_df) - 1
  col_data <- data.frame(
    condition = factor(c(rep("DMSO", n_control), rep(comparisons[[cond_name]]["treatment"], n_treatment)),
                       levels = c("DMSO", comparisons[[cond_name]]["treatment"])),
    row.names = colnames(counts_data)[-1]  # Exclude gene_name
  )
  
  # Clean condition names for DESeq2
  levels(col_data$condition) <- gsub("\\.", "_", levels(col_data$condition))
  
  cat("col_data for", cond_name, "vs DMSO - Dimensions:", dim(col_data), "\n")
  
  # DESeq2 results file path (updated to use dots)
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(comparisons[[cond_name]]["treatment"], "_vs_DMSO_deseq_results.csv"))
  cat("DESeq2 results file for", cond_name, ":", deseq_results_file, "\n")
  
  # Check if the DESeq2 results file exists
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping heatmap for", cond_name, "\n")
    next
  }
  
  # Generate heatmap
  heatmap_plot <- generate_heatmap(counts_data = counts_data,
                                   col_data = col_data,
                                   deseq_results_file = deseq_results_file,
                                   comparison_name = comparison_titles[[cond_name]],
                                   cond_name = cond_name,
                                   padj_threshold = 0.05,
                                   lfc_threshold = 1,
                                   output_dir = output_dir)
  
  cat("Generated heatmap for", cond_name, "vs DMSO - Saved to", output_dir, "\n")
}

# Save relevant data
save(merged_counts, condition_dfs, file = "heatmap_data_Cell_line_Genotype_1_dmso.rds")
cat("Saved heatmap data to heatmap_data_Cell_line_Genotype_1_dmso.rds\n")



#######################
#Pathway Analysis
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(enrichplot)
library(fgsea)
library(reshape2)

# Ensure output directory exists
top_output_dir <- "Top_Gene_and_Pathway_Results"
if (!dir.exists(top_output_dir)) {
  dir.create(top_output_dir, recursive = TRUE)
}

# Define comparisons (using DMSO as control, with dots)
comparisons <- list(
  NT = c(control = "DMSO", treatment = "NT"),
  Compound_4_100nM = c(control = "DMSO", treatment = "Compound.4_100nM"),
  Compound_6_1uM = c(control = "DMSO", treatment = "Compound.6_1uM"),
  Compound_5_1uM = c(control = "DMSO", treatment = "Compound.5_1uM"),
  Compound_7_100nM = c(control = "DMSO", treatment = "Compound.7_100nM"),
  Compound_3_100nM = c(control = "DMSO", treatment = "Compound.3_100nM"),
  Compound_1_100nM = c(control = "DMSO", treatment = "Compound.1_100nM"),
  Compound_2_100nM = c(control = "DMSO", treatment = "Compound.2_100nM")
)

# Initialize a data frame to store the counts of up/down-regulated genes
significant_counts <- data.frame(Comparison = character(), Upregulated = integer(), Downregulated = integer())

# Perform pairwise comparisons and count upregulated/downregulated genes
for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep="_")
  cat("Processing comparison:", contrast_name, "\n")
  
  # Load DESeq2 results
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(comparison["treatment"], "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, row.names = 1)
  
  # Convert columns to numeric
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  
  # Count significantly upregulated genes (log2FoldChange > 0, padj < 0.05)
  upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05 & !is.na(res$padj), na.rm = TRUE)
  
  # Count significantly downregulated genes (log2FoldChange < 0, padj < 0.05)
  downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05 & !is.na(res$padj), na.rm = TRUE)
  
  # Append the counts to the data frame
  significant_counts <- rbind(significant_counts, data.frame(Comparison = contrast_name, Upregulated = upregulated, Downregulated = downregulated))
}

# Print the counts to verify
print(significant_counts)
# Save the counts to a CSV file
write.csv(significant_counts, file = "significant_gene_counts.csv", row.names = FALSE)
# Reshape the data to long format for ggplot2
long_data <- melt(significant_counts, id.vars = "Comparison", variable.name = "Regulation", value.name = "Count")

# Generate the stacked bar plot
ggplot(long_data, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Significantly Upregulated and Downregulated Genes",
       x = "Comparison", y = "Number of Genes") +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "firebrick")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("significant_gene_counts_stacked_barplot.png", width = 8, height = 6)

# Create output directory for pathway results
output_dir <- "Pathway_Analysis_Cell_line_Genotype_1_DMSO"
if (!dir.exists(output_dir)) dir.create(output_dir)



# Initialize matrices for top genes (rows: comparisons, columns: top 10 genes)
top_up_genes_matrix <- matrix(NA, nrow = length(comparisons), ncol = 10)
top_down_genes_matrix <- matrix(NA, nrow = length(comparisons), ncol = 10)
rownames(top_up_genes_matrix) <- names(comparisons)
rownames(top_down_genes_matrix) <- names(comparisons)
colnames(top_up_genes_matrix) <- paste0("Gene_", 1:10)
colnames(top_down_genes_matrix) <- paste0("Gene_", 1:10)

for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep = "_")
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(contrast_name, "_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", contrast_name, "- Skipping\n")
    next
  }
  
  # Load DESeq2 results
  res <- read.csv(deseq_results_file, row.names = 1)
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- sub("\\.\\d+$", "", rownames(res_df))
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  # Filter and sort upregulated genes
  up_genes_df <- res_df[res_df$log2FoldChange > 1 & res_df$padj < 0.05, ]
  if (nrow(up_genes_df) > 0) {
    up_genes_df <- up_genes_df[order(up_genes_df$log2FoldChange, decreasing = TRUE), ]
    top_up_genes <- head(up_genes_df$Symbol, 10)
    top_up_genes_matrix[cond_name, 1:min(10, length(top_up_genes))] <- top_up_genes
  }
  
  # Filter and sort downregulated genes
  down_genes_df <- res_df[res_df$log2FoldChange < -1 & res_df$padj < 0.05, ]
  if (nrow(down_genes_df) > 0) {
    down_genes_df <- down_genes_df[order(down_genes_df$log2FoldChange, decreasing = FALSE), ]
    top_down_genes <- head(down_genes_df$Symbol, 10)
    top_down_genes_matrix[cond_name, 1:min(10, length(top_down_genes))] <- top_down_genes
  }
  
  cat("Extracted top genes for", contrast_name, "\n")
}

# Save top genes to CSV
write.csv(top_up_genes_matrix,
          file = file.path(top_output_dir, "top10_upregulated_genes.csv"),
          quote = TRUE, na = "")
write.csv(top_down_genes_matrix,
          file = file.path(top_output_dir, "top10_downregulated_genes.csv"),
          quote = TRUE, na = "")

# Step 2: Extract top 10 pathways for each comparison
top_pathways_matrix <- matrix(NA, nrow = length(comparisons), ncol = 10)
rownames(top_pathways_matrix) <- names(comparisons)
colnames(top_pathways_matrix) <- paste0("Pathway_", 1:10)

for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep = "_")
  deseq_results_file <- file.path("DESeq_Results", paste0(contrast_name, "_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", contrast_name, "- Skipping\n")
    next
  }
  
  # Load DESeq2 results
  res <- read.csv(deseq_results_file, row.names = 1)
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- sub("\\.\\d+$", "", rownames(res_df))
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  # Split into up- and downregulated genes
  up_genes <- res_df$Symbol[res_df$log2FoldChange > 0.5 & res_df$padj < 0.1]
  down_genes <- res_df$Symbol[res_df$log2FoldChange < -0.5 & res_df$padj < 0.1]
  
  # Map gene symbols to Entrez IDs
  entrez_ids_up <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  entrez_ids_down <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  # Rank genes for Hallmark GSEA
  ranked_genes <- res_df$log2FoldChange
  names(ranked_genes) <- res_df$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Collect pathways
  condition_pathways <- data.frame(pathway = character(), padj = numeric(), stringsAsFactors = FALSE)
  
  # Hallmark GSEA
  msigdb <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
  all_genes_in_sets <- unique(gene_sets_hallmark$gene)
  
  ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
  if (length(ranked_genes_filtered) > 0) {
    gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                         stats = ranked_genes_filtered,
                         minSize = 15,
                         maxSize = 500,
                         nPermSimple = 1000)
    
    up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < 0.1, ]
    down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < 0.1, ]
    
    if (nrow(up_enriched_hallmark) > 0) {
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("Hallmark_", up_enriched_hallmark$pathway),
        padj = up_enriched_hallmark$padj
      ))
    }
    if (nrow(down_enriched_hallmark) > 0) {
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("Hallmark_", down_enriched_hallmark$pathway),
        padj = down_enriched_hallmark$padj
      ))
    }
  }
  
  # KEGG Enrichment for upregulated genes
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_ids_up$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      kegg_up_df <- as.data.frame(kegg_up)
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("KEGG_", kegg_up_df$Description),
        padj = kegg_up_df$p.adjust
      ))
    }
  }
  
  # KEGG Enrichment for downregulated genes
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_ids_down$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      kegg_down_df <- as.data.frame(kegg_down)
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("KEGG_", kegg_down_df$Description),
        padj = kegg_down_df$p.adjust
      ))
    }
  }
  
  # GO Enrichment (Biological Process) for upregulated genes
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    go_up <- tryCatch({
      enrichGO(gene = entrez_ids_up$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = 0.1,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      go_up_df <- as.data.frame(go_up)
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("GO_BP_", go_up_df$Description),
        padj = go_up_df$p.adjust
      ))
    }
  }
  
  # GO Enrichment (Biological Process) for downregulated genes
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    go_down <- tryCatch({
      enrichGO(gene = entez_ids_down$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = 0.1,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      go_down_df <- as.data.frame(go_down)
      condition_pathways <- rbind(condition_pathways, data.frame(
        pathway = paste0("GO_BP_", go_down_df$Description),
        padj = go_down_df$p.adjust
      ))
    }
  }
  
  # Select top 10 pathways by padj
  if (nrow(condition_pathways) > 0) {
    condition_pathways <- condition_pathways[order(condition_pathways$padj), ]
    top_pathways <- head(condition_pathways$pathway, 10)
    top_pathways_matrix[cond_name, 1:min(10, length(top_pathways))] <- top_pathways
  }
  
  cat("Extracted top pathways for", contrast_name, "\n")
}

# Save top pathways to CSV
write.csv(top_pathways_matrix,
          file = file.path(top_output_dir, "top10_pathways.csv"),
          quote = TRUE, na = "")

cat("Top genes and pathways CSV files saved in", top_output_dir, "\n")


# Separate analysis for upregulated and downregulated genes
for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep="_")
  cat("Performing pathway analysis for up/downregulated genes for:", contrast_name, "\n")
  
  # Load DESeq2 results
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(comparison["treatment"], "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, row.names = 1)
  
  # Convert columns to numeric
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  
  # Add gene symbols
  res_annot <- as.data.frame(res)
  res_annot$Symbol <- rownames(res_annot)
  
  # Clean gene symbols
  res_annot <- res_annot[!is.na(res_annot$Symbol) & res_annot$Symbol != "", ]
  res_annot$Symbol <- trimws(res_annot$Symbol)
  res_annot$Symbol <- gsub("[^[:alnum:]_-]", "", res_annot$Symbol)
  
  # Check for duplicates
  if (any(duplicated(res_annot$Symbol))) {
    cat("Warning: Duplicate gene symbols found for", contrast_name, "- Making unique\n")
    res_annot$Symbol <- make.unique(res_annot$Symbol)
  }
  
  # Upregulated genes
  upregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange > 0 & res_annot$padj < 0.05 & !is.na(res_annot$padj)]
  
  # Check if there are valid gene symbols to proceed
  if (length(upregulated_genes) > 0) {
    up_entrez_ids <- tryCatch({
      bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("Error mapping upregulated gene symbols for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(up_entrez_ids) && nrow(up_entrez_ids) > 0) {
      # Perform GO and KEGG enrichment analysis for upregulated genes
      go_up_enrich <- enrichGO(gene = up_entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
      kegg_up_enrich <- enrichKEGG(gene = up_entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      kegg_up_enrich <- setReadable(kegg_up_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    } else {
      cat("No upregulated genes mapped to Entrez IDs for", contrast_name, "\n")
      go_up_enrich <- NULL
      kegg_up_enrich <- NULL
    }
  } else {
    cat("No upregulated genes found for:", contrast_name, "\n")
    go_up_enrich <- NULL
    kegg_up_enrich <- NULL
  }
  
  # Downregulated genes
  downregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange < 0 & res_annot$padj < 0.05 & !is.na(res_annot$padj)]
  
  # Check if there are valid gene symbols to proceed
  if (length(downregulated_genes) > 0) {
    down_entrez_ids <- tryCatch({
      bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("Error mapping downregulated gene symbols for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(down_entrez_ids) && nrow(down_entrez_ids) > 0) {
      # Perform GO and KEGG enrichment analysis for downregulated genes
      go_down_enrich <- enrichGO(gene = down_entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
      kegg_down_enrich <- enrichKEGG(gene = down_entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      kegg_down_enrich <- setReadable(kegg_down_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    } else {
      cat("No downregulated genes mapped to Entrez IDs for", contrast_name, "\n")
      go_down_enrich <- NULL
      kegg_down_enrich <- NULL
    }
  } else {
    cat("No downregulated genes found for:", contrast_name, "\n")
    go_down_enrich <- NULL
    kegg_down_enrich <- NULL
  }
  
  # Save upregulated and downregulated pathway results
  if (!is.null(go_up_enrich)) {
    write.table(as.data.frame(go_up_enrich), file = file.path(output_dir, paste0("GO_upregulated_", contrast_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(go_down_enrich)) {
    write.table(as.data.frame(go_down_enrich), file = file.path(output_dir, paste0("GO_downregulated_", contrast_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(kegg_up_enrich)) {
    write.table(as.data.frame(kegg_up_enrich), file = file.path(output_dir, paste0("KEGG_upregulated_", contrast_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(kegg_down_enrich)) {
    write.table(as.data.frame(kegg_down_enrich), file = file.path(output_dir, paste0("KEGG_downregulated_", contrast_name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  # Plot top GO and KEGG pathways for both upregulated and downregulated genes
  if (!is.null(go_up_enrich) && nrow(as.data.frame(go_up_enrich)) > 0) {
    go_up_plot <- barplot(go_up_enrich, showCategory = 10, title = paste0("Top GO Terms (Upregulated): ", contrast_name))
    print(go_up_plot)
    ggsave(file.path(output_dir, paste0("go_upregulated_", contrast_name, ".png")), go_up_plot, width = 8, height = 6)
  }
  
  if (!is.null(go_down_enrich) && nrow(as.data.frame(go_down_enrich)) > 0) {
    go_down_plot <- barplot(go_down_enrich, showCategory = 10, title = paste0("Top GO Terms (Downregulated): ", contrast_name))
    print(go_down_plot)
    ggsave(file.path(output_dir, paste0("go_downregulated_", contrast_name, ".png")), go_down_plot, width = 8, height = 6)
  }
  
  if (!is.null(kegg_up_enrich) && nrow(as.data.frame(kegg_up_enrich)) > 0) {
    kegg_up_plot <- barplot(kegg_up_enrich, showCategory = 10, title = paste0("Top KEGG Pathways (Upregulated): ", contrast_name))
    print(kegg_up_plot)
    ggsave(file.path(output_dir, paste0("kegg_upregulated_", contrast_name, ".png")), kegg_up_plot, width = 8, height = 6)
  }
  
  if (!is.null(kegg_down_enrich) && nrow(as.data.frame(kegg_down_enrich)) > 0) {
    kegg_down_plot <- barplot(kegg_down_enrich, showCategory = 10, title = paste0("Top KEGG Pathways (Downregulated): ", contrast_name))
    print(kegg_down_plot)
    ggsave(file.path(output_dir, paste0("kegg_downregulated_", contrast_name, ".png")), kegg_down_plot, width = 8, height = 6)
  }
}

# Cnet plots for each comparison (GO and KEGG for upregulated and downregulated genes)
for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep="_")
  cat("Generating cnet plots for:", contrast_name, "\n")
  
  # Load DESeq2 results
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(comparison["treatment"], "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, row.names = 1)
  
  # Convert columns to numeric
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  
  # Add gene symbols
  res_annot <- as.data.frame(res)
  res_annot$Symbol <- rownames(res_annot)
  
  # Clean gene symbols
  res_annot <- res_annot[!is.na(res_annot$Symbol) & res_annot$Symbol != "", ]
  res_annot$Symbol <- trimws(res_annot$Symbol)
  res_annot$Symbol <- gsub("[^[:alnum:]_-]", "", res_annot$Symbol)
  
  # Check for duplicates
  if (any(duplicated(res_annot$Symbol))) {
    cat("Warning: Duplicate gene symbols found for", contrast_name, "- Making unique\n")
    res_annot$Symbol <- make.unique(res_annot$Symbol)
  }
  
  # Create named log2FoldChange vector for cnet plots
  fold_changes <- res_annot$log2FoldChange
  names(fold_changes) <- res_annot$Symbol
  
  # Upregulated genes
  upregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange > 0 & res_annot$padj < 0.05 & !is.na(res_annot$padj)]
  if (length(upregulated_genes) > 0) {
    up_entrez_ids <- tryCatch({
      bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("Error mapping upregulated gene symbols for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(up_entrez_ids) && nrow(up_entrez_ids) > 0) {
      go_up_enrich <- enrichGO(gene = up_entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
      kegg_up_enrich <- enrichKEGG(gene = up_entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      kegg_up_enrich <- setReadable(kegg_up_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      # Cnet plots for upregulated
      if (!is.null(go_up_enrich) && nrow(as.data.frame(go_up_enrich)) > 0) {
        cnet_go_up <- cnetplot(go_up_enrich, categorySize = "pvalue", foldChange = fold_changes, showCategory = 5)
        ggsave(file.path(output_dir, paste0("cnet_go_upregulated_", contrast_name, ".png")), cnet_go_up, width = 10, height = 8)
      }
      if (!is.null(kegg_up_enrich) && nrow(as.data.frame(kegg_up_enrich)) > 0) {
        cnet_kegg_up <- cnetplot(kegg_up_enrich, categorySize = "pvalue", foldChange = fold_changes, showCategory = 5)
        ggsave(file.path(output_dir, paste0("cnet_kegg_upregulated_", contrast_name, ".png")), cnet_kegg_up, width = 10, height = 8)
      }
    }
  }
  
  # Downregulated genes
  downregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange < 0 & res_annot$padj < 0.05 & !is.na(res_annot$padj)]
  if (length(downregulated_genes) > 0) {
    down_entrez_ids <- tryCatch({
      bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("Error mapping downregulated gene symbols for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(down_entrez_ids) && nrow(down_entrez_ids) > 0) {
      go_down_enrich <- enrichGO(gene = down_entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
      kegg_down_enrich <- enrichKEGG(gene = down_entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      kegg_down_enrich <- setReadable(kegg_down_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      # Cnet plots for downregulated
      if (!is.null(go_down_enrich) && nrow(as.data.frame(go_down_enrich)) > 0) {
        cnet_go_down <- cnetplot(go_down_enrich, categorySize = "pvalue", foldChange = fold_changes, showCategory = 5)
        ggsave(file.path(output_dir, paste0("cnet_go_downregulated_", contrast_name, ".png")), cnet_go_down, width = 10, height = 8)
      }
      if (!is.null(kegg_down_enrich) && nrow(as.data.frame(kegg_down_enrich)) > 0) {
        cnet_kegg_down <- cnetplot(kegg_down_enrich, categorySize = "pvalue", foldChange = fold_changes, showCategory = 5)
        ggsave(file.path(output_dir, paste0("cnet_kegg_downregulated_", contrast_name, ".png")), cnet_kegg_down, width = 10, height = 8)
      }
    }
  }
}







#Up and down separately
# Load required libraries
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

# Create output directory
top_output_dir <- "Top_Pathways"
if (!dir.exists(top_output_dir)) dir.create(top_output_dir)

#Initialize matrices for top 10 pathways
top_up_pathways_matrix <- matrix(NA, nrow = length(comparisons), ncol = 10)
rownames(top_up_pathways_matrix) <- names(comparisons)
colnames(top_up_pathways_matrix) <- paste0("Pathway_", 1:10)

top_down_pathways_matrix <- matrix(NA, nrow = length(comparisons), ncol = 10)
rownames(top_down_pathways_matrix) <- names(comparisons)
colnames(top_down_pathways_matrix) <- paste0("Pathway_", 1:10)

for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep = "_")
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(contrast_name, "_deseq_results.csv"))
  
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", contrast_name, "- Skipping\n")
    next
  }
  
  # Load DESeq2 results
  res <- read.csv(deseq_results_file, row.names = 1)
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$padj <- as.numeric(as.character(res$padj))
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- trimws(toupper(sub("\\.\\d+$", "", rownames(res_df))))
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  
  # Split into up- and downregulated genes
  up_genes <- res_df$Symbol[res_df$log2FoldChange > 0.5 & res_df$padj < 0.1]
  down_genes <- res_df$Symbol[res_df$log2FoldChange < -0.5 & res_df$padj < 0.1]
  
  # Map gene symbols to Entrez IDs
  entrez_ids_up <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  entrez_ids_down <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  # Rank genes for Hallmark GSEA
  ranked_genes <- res_df$log2FoldChange
  names(ranked_genes) <- res_df$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Initialize pathway dataframes
  up_condition_pathways <- data.frame(pathway = character(), padj = numeric(), stringsAsFactors = FALSE)
  down_condition_pathways <- data.frame(pathway = character(), padj = numeric(), stringsAsFactors = FALSE)
  
  # Hallmark GSEA
  msigdb <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
  all_genes_in_sets <- unique(gene_sets_hallmark$gene)
  
  ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
  if (length(ranked_genes_filtered) > 0) {
    gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                         stats = ranked_genes_filtered,
                         minSize = 15,
                         maxSize = 500,
                         nPermSimple = 1000)
    
    up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < 0.1, ]
    down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < 0.1, ]
    
    if (nrow(up_enriched_hallmark) > 0) {
      up_condition_pathways <- rbind(up_condition_pathways, data.frame(
        pathway = paste0("Hallmark_", up_enriched_hallmark$pathway),
        padj = up_enriched_hallmark$padj
      ))
    }
    if (nrow(down_enriched_hallmark) > 0) {
      down_condition_pathways <- rbind(down_condition_pathways, data.frame(
        pathway = paste0("Hallmark_", down_enriched_hallmark$pathway),
        padj = down_enriched_hallmark$padj
      ))
    }
  }
  
  # KEGG Enrichment for upregulated genes
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_ids_up$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("Error in KEGG enrichment (up) for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      kegg_up_df <- as.data.frame(kegg_up)
      up_condition_pathways <- rbind(up_condition_pathways, data.frame(
        pathway = paste0("KEGG_", kegg_up_df$Description),
        padj = kegg_up_df$p.adjust
      ))
    }
  }
  
  # KEGG Enrichment for downregulated genes
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_ids_down$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("Error in KEGG enrichment (down) for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      kegg_down_df <- as.data.frame(kegg_down)
      down_condition_pathways <- rbind(down_condition_pathways, data.frame(
        pathway = paste0("KEGG_", kegg_down_df$Description),
        padj = kegg_down_df$p.adjust
      ))
    }
  }
  
  # GO Enrichment (Biological Process) for upregulated genes
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    go_up <- tryCatch({
      enrichGO(gene = entrez_ids_up$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = 0.1,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("Error in GO enrichment (up) for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      go_up_df <- as.data.frame(go_up)
      up_condition_pathways <- rbind(up_condition_pathways, data.frame(
        pathway = paste0("GO_BP_", go_up_df$Description),
        padj = go_up_df$p.adjust
      ))
    }
  }
  
  # GO Enrichment (Biological Process) for downregulated genes
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    go_down <- tryCatch({
      enrichGO(gene = entrez_ids_down$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = 0.1,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("Error in GO enrichment (down) for", contrast_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      go_down_df <- as.data.frame(go_down)
      down_condition_pathways <- rbind(down_condition_pathways, data.frame(
        pathway = paste0("GO_BP_", go_down_df$Description),
        padj = go_down_df$p.adjust
      ))
    }
  }
  
  # Select top 10 pathways for upregulated genes
  if (nrow(up_condition_pathways) > 0) {
    up_condition_pathways <- up_condition_pathways[order(up_condition_pathways$padj), ]
    top_up_pathways <- head(up_condition_pathways$pathway, 10)
    top_up_pathways_matrix[cond_name, 1:min(10, length(top_up_pathways))] <- top_up_pathways
  }
  
  # Select top 10 pathways for downregulated genes
  if (nrow(down_condition_pathways) > 0) {
    down_condition_pathways <- down_condition_pathways[order(down_condition_pathways$padj), ]
    top_down_pathways <- head(down_condition_pathways$pathway, 10)
    top_down_pathways_matrix[cond_name, 1:min(10, length(top_down_pathways))] <- top_down_pathways
  }
  
  cat("Extracted top pathways for", contrast_name, "\n")
}

# Save top pathways to CSV
write.csv(top_up_pathways_matrix,
          file = file.path(top_output_dir, "top10_upregulated_pathways.csv"),
          quote = TRUE, na = "")
write.csv(top_down_pathways_matrix,
          file = file.path(top_output_dir, "top10_downregulated_pathways.csv"),
          quote = TRUE, na = "")

cat("Top pathways CSV files saved in", top_output_dir, "\n")




##########################
#Network
library(msigdbr)
library(enrichplot)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)

# Create output directory
output_dir <- "Cnet_Pathway_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)


# Lists to store enrichment results
all_enrich_results <- list()

# Enrichment for up/downregulated genes
for (cond_name in names(comparisons)) {
  comparison <- comparisons[[cond_name]]
  contrast_name <- paste(comparison["treatment"], "vs", comparison["control"], sep = "_")
  cat("Performing pathway analysis for up/downregulated genes for:", contrast_name, "\n")
  
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(contrast_name, "_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, row.names = 1)
  res$log2FoldChange <- as.numeric(as.character(res$log2FoldChange))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$padj <- as.numeric(as.character(res$padj))
  
  res_annot <- as.data.frame(res)
  res_annot$Symbol <- sub("\\.\\d+$", "", rownames(res_annot))
  
  # Use relaxed thresholds
  upregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange > 0.5 & res_annot$padj < 0.1]
  downregulated_genes <- res_annot$Symbol[res_annot$log2FoldChange < -0.5 & res_annot$padj < 0.1]
  
  # Map symbols to Entrez IDs for GO and KEGG
  up_entrez_ids <- tryCatch({
    bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  down_entrez_ids <- tryCatch({
    bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", contrast_name, ":", e$message, "\n")
    NULL
  })
  
  # Rank genes by log2FoldChange for MsigDB (Hallmark) GSEA
  ranked_genes <- res_annot$log2FoldChange
  names(ranked_genes) <- res_annot$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Upregulated genes enrichment
  if (length(upregulated_genes) > 0 && !is.null(up_entrez_ids)) {
    tryCatch({
      # GO enrichment
      go_up_enrich <- enrichGO(gene = up_entrez_ids$ENTREZID, 
                               OrgDb = org.Hs.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "fdr", 
                               pvalueCutoff = 0.05, 
                               readable = TRUE)
      
      if (!is.null(go_up_enrich) && nrow(as.data.frame(go_up_enrich)) > 0) {
        go_up_df <- as.data.frame(go_up_enrich)
        go_up_df_standard <- data.frame(
          Description = go_up_df$Description,
          geneID = go_up_df$geneID,
          pvalue = go_up_df$pvalue,
          Condition = cond_name,
          GeneSetType = "GO",
          Direction = "Up",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("GO_Up_", cond_name)]] <- go_up_df_standard
      }
      
      # KEGG enrichment
      kegg_up_enrich <- enrichKEGG(gene = up_entrez_ids$ENTREZID, 
                                   organism = 'hsa', 
                                   pvalueCutoff = 0.05)
      kegg_up_enrich <- setReadable(kegg_up_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      if (!is.null(kegg_up_enrich) && nrow(as.data.frame(kegg_up_enrich)) > 0) {
        kegg_up_df <- as.data.frame(kegg_up_enrich)
        kegg_up_df_standard <- data.frame(
          Description = kegg_up_df$Description,
          geneID = kegg_up_df$geneID,
          pvalue = kegg_up_df$pvalue,
          Condition = cond_name,
          GeneSetType = "KEGG",
          Direction = "Up",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("KEGG_Up_", cond_name)]] <- kegg_up_df_standard
      }
    }, error = function(e) {
      cat("Error in upregulated enrichment for", contrast_name, ":", e$message, "\n")
    })
  } else {
    cat("No upregulated genes found for:", contrast_name, "\n")
  }
  
  # Downregulated genes enrichment
  if (length(downregulated_genes) > 0 && !is.null(down_entrez_ids)) {
    tryCatch({
      # GO enrichment
      go_down_enrich <- enrichGO(gene = down_entrez_ids$ENTREZID, 
                                 OrgDb = org.Hs.eg.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "fdr", 
                                 pvalueCutoff = 0.05, 
                                 readable = TRUE)
      
      if (!is.null(go_down_enrich) && nrow(as.data.frame(go_down_enrich)) > 0) {
        go_down_df <- as.data.frame(go_down_enrich)
        go_down_df_standard <- data.frame(
          Description = go_down_df$Description,
          geneID = go_down_df$geneID,
          pvalue = go_down_df$pvalue,
          Condition = cond_name,
          GeneSetType = "GO",
          Direction = "Down",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("GO_Down_", cond_name)]] <- go_down_df_standard
      }
      
      # KEGG enrichment
      kegg_down_enrich <- enrichKEGG(gene = down_entrez_ids$ENTREZID, 
                                     organism = 'hsa', 
                                     pvalueCutoff = 0.05)
      kegg_up_enrich <- setReadable(kegg_up_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      if (!is.null(kegg_down_enrich) && nrow(as.data.frame(kegg_down_enrich)) > 0) {
        kegg_down_df <- as.data.frame(kegg_down_enrich)
        kegg_down_df_standard <- data.frame(
          Description = kegg_down_df$Description,
          geneID = kegg_down_df$geneID,
          pvalue = kegg_down_df$pvalue,
          Condition = cond_name,
          GeneSetType = "KEGG",
          Direction = "Down",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("KEGG_Down_", cond_name)]] <- kegg_down_df_standard
      }
    }, error = function(e) {
      cat("Error in downregulated enrichment for", contrast_name, ":", e$message, "\n")
    })
  } else {
    cat("No downregulated genes found for:", contrast_name, "\n")
  }
  
  # MsigDB (Hallmark) GSEA
  tryCatch({
    msigdb <- msigdbr(species = "Homo sapiens", category = "H")
    gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
    all_genes_in_sets <- unique(gene_sets_hallmark$gene)
    
    ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
    if (length(ranked_genes_filtered) > 0) {
      gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                           stats = ranked_genes_filtered,
                           minSize = 15,
                           maxSize = 500,
                           nPermSimple = 1000)
      
      up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < 0.1, ]
      down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < 0.1, ]
      
      if (nrow(up_enriched_hallmark) > 0) {
        hallmark_up_df <- as.data.frame(up_enriched_hallmark)
        hallmark_up_df$Description <- hallmark_up_df$pathway
        hallmark_up_df$geneID <- sapply(hallmark_up_df$leadingEdge, function(x) paste(x, collapse = "/"))
        hallmark_up_df_standard <- data.frame(
          Description = hallmark_up_df$Description,
          geneID = hallmark_up_df$geneID,
          pvalue = hallmark_up_df$pval,
          Condition = cond_name,
          GeneSetType = "MsigDB",
          Direction = "Up",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("MsigDB_Up_", cond_name)]] <- hallmark_up_df_standard
      }
      
      if (nrow(down_enriched_hallmark) > 0) {
        hallmark_down_df <- as.data.frame(down_enriched_hallmark)
        hallmark_down_df$Description <- hallmark_down_df$pathway
        hallmark_down_df$geneID <- sapply(hallmark_down_df$leadingEdge, function(x) paste(x, collapse = "/"))
        hallmark_down_df_standard <- data.frame(
          Description = hallmark_down_df$Description,
          geneID = hallmark_down_df$geneID,
          pvalue = hallmark_down_df$pval,
          Condition = cond_name,
          GeneSetType = "MsigDB",
          Direction = "Down",
          stringsAsFactors = FALSE
        )
        all_enrich_results[[paste0("MsigDB_Down_", cond_name)]] <- hallmark_down_df_standard
      }
    }
  }, error = function(e) {
    cat("Error in MsigDB GSEA for", contrast_name, ":", e$message, "\n")
  })
  
  # Combine all enrichment results for this comparison
  if (length(all_enrich_results) == 0) {
    cat("No enrichment results available for", contrast_name, ". Skipping cnet plot.\n")
    next
  }
  
  combined_results <- do.call(rbind, all_enrich_results[grepl(cond_name, names(all_enrich_results))])
  
  # Merge pathways that appear in multiple conditions (for this comparison)
  pathway_groups <- split(combined_results, combined_results$Description)
  
  merged_results <- lapply(pathway_groups, function(group) {
    conditions <- unique(group$Condition)
    merged_pathway_id <- paste0(group$Description[1], "_", paste(conditions, collapse = "_"))
    
    all_genes <- unique(unlist(strsplit(paste(group$geneID, collapse = "/"), "/")))
    
    min_pval <- min(group$pvalue)
    
    gene_set_type <- group$GeneSetType[1]
    
    directions <- unique(group$Direction)
    direction <- if (length(directions) == 1) directions else "Mixed"
    
    data.frame(
      Description = group$Description[1],
      geneID = paste(all_genes, collapse = "/"),
      pvalue = min_pval,
      Condition = paste(conditions, collapse = "_"),
      GeneSetType = gene_set_type,
      Direction = direction,
      PathwayID = merged_pathway_id,
      ConditionCount = length(conditions),
      stringsAsFactors = FALSE
    )
  })
  
  merged_results <- do.call(rbind, merged_results)
  
  # Prepare gene-to-pathway mappings
  gene_pathway_list <- list()
  for (i in 1:nrow(merged_results)) {
    genes <- unlist(strsplit(merged_results$geneID[i], "/"))
    for (gene in genes) {
      if (!(gene %in% names(gene_pathway_list))) {
        gene_pathway_list[[gene]] <- c()
      }
      gene_pathway_list[[gene]] <- c(gene_pathway_list[[gene]], merged_results$PathwayID[i])
    }
  }
  
  # Create a pathway-to-metadata mapping
  pathway_metadata <- data.frame(
    PathwayID = merged_results$PathwayID,
    GeneSetType = merged_results$GeneSetType,
    Direction = merged_results$Direction,
    pvalue = merged_results$pvalue,
    ConditionCount = merged_results$ConditionCount,
    stringsAsFactors = FALSE
  )
  
  # Create a gene-to-direction mapping
  gene_direction <- list()
  for (i in 1:nrow(merged_results)) {
    genes <- unlist(strsplit(merged_results$geneID[i], "/"))
    direction <- merged_results$Direction[i]
    for (gene in genes) {
      if (gene %in% names(gene_direction)) {
        existing_pval <- gene_direction[[gene]]$pvalue
        if (merged_results$pvalue[i] < existing_pval) {
          gene_direction[[gene]] <- list(direction = direction, pvalue = merged_results$pvalue[i])
        }
      } else {
        gene_direction[[gene]] <- list(direction = direction, pvalue = merged_results$pvalue[i])
      }
    }
  }
  
  # Debug: Check gene_pathway_list
  cat("Number of genes in gene_pathway_list:", length(gene_pathway_list), "\n")
  cat("Sample gene_pathway_list entries:", head(names(gene_pathway_list), 5), "\n")
  
  # Check if gene_pathway_list is empty
  if (length(gene_pathway_list) == 0) {
    cat("Error: gene_pathway_list is empty for", contrast_name, ". Skipping cnet plot.\n")
    next
  }
  
  # Filter gene_pathway_list to include only the top 50 pathways
  all_pathways <- unique(unlist(gene_pathway_list))
  top_pathways <- head(all_pathways, 50)
  cat("Number of pathways selected for", contrast_name, ":", length(top_pathways), "\n")
  
  # Subset gene_pathway_list to include only genes associated with the top pathways
  filtered_gene_pathway_list <- list()
  for (gene in names(gene_pathway_list)) {
    pathways <- gene_pathway_list[[gene]]
    pathways_to_keep <- pathways[pathways %in% top_pathways]
    if (length(pathways_to_keep) > 0) {
      filtered_gene_pathway_list[[gene]] <- pathways_to_keep
    }
  }
  
  # Debug: Check filtered_gene_pathway_list
  cat("Number of genes in filtered_gene_pathway_list:", length(filtered_gene_pathway_list), "\n")
  cat("Sample filtered_gene_pathway_list entries:", head(names(filtered_gene_pathway_list), 5), "\n")
  
  # Check if filtered_gene_pathway_list is empty
  if (length(filtered_gene_pathway_list) == 0) {
    cat("Error: filtered_gene_pathway_list is empty after filtering pathways for", contrast_name, ". Skipping cnet plot.\n")
    next
  }
  
  # Generate the base cnet plot with filtered pathways
  cnet <- cnetplot(
    x = filtered_gene_pathway_list,
    showCategory = length(top_pathways),
    categorySize = "pvalue",
    color.params = list(foldChange = NULL)
  )
  
  # Extract node data
  cnet_data <- ggplot_build(cnet)$data[[2]]
  
  # Debug: Check contents of cnet_data
  cat("Number of rows in cnet_data:", nrow(cnet_data), "\n")
  
  # Determine unique nodes in cnet_data
  if ("group" %in% colnames(cnet_data)) {
    unique_nodes <- unique(cnet_data$group)
    cat("Number of unique nodes (based on group):", length(unique_nodes), "\n")
    cat("Sample group values in cnet_data:", head(unique_nodes, 5), "\n")
  } else {
    all_genes <- names(filtered_gene_pathway_list)
    all_pathways <- unique(unlist(filtered_gene_pathway_list))
    expected_nodes <- c(all_pathways, all_genes)
    unique_nodes <- expected_nodes
    cat("Number of unique nodes (based on expected):", length(unique_nodes), "\n")
  }
  
  # Reconstruct node names if necessary
  if (all(is.na(cnet_data$name)) || all(cnet_data$name == "")) {
    cat("Warning: cnet_data$name is empty. Manually reconstructing node names...\n")
    
    all_genes <- names(filtered_gene_pathway_list)
    all_pathways <- unique(unlist(filtered_gene_pathway_list))
    expected_nodes <- c(all_pathways, all_genes)
    
    if (length(unique_nodes) != length(expected_nodes)) {
      cat("Warning: Number of unique nodes (", length(unique_nodes), ") does not match expected nodes (", length(expected_nodes), "). Using expected nodes.\n")
      unique_nodes <- expected_nodes
    }
    
    node_name_mapping <- rep(NA, length(unique_nodes))
    node_name_mapping[1:length(all_pathways)] <- all_pathways
    node_name_mapping[(length(all_pathways) + 1):length(unique_nodes)] <- all_genes
    
    if ("group" %in% colnames(cnet_data)) {
      group_counts <- table(cnet_data$group)
      cat("Group distribution in cnet_data:\n")
      print(group_counts)
      
      num_pathways <- length(all_pathways)
      num_genes <- nrow(cnet_data) - num_pathways
      indices <- rep(NA, nrow(cnet_data))
      indices[cnet_data$group == 1] <- 1:num_pathways
      indices[cnet_data$group == 2] <- (num_pathways + 1):(num_pathways + num_genes)
      
      cnet_data$name <- node_name_mapping[indices]
    } else {
      indices <- 1:nrow(cnet_data)
      cnet_data$name <- node_name_mapping[indices]
    }
  }
  
  # Debug: Check contents of pathway_metadata and cnet_data
  cat("Number of pathways in pathway_metadata:", length(pathway_metadata$PathwayID), "\n")
  cat("Sample PathwayIDs:", head(pathway_metadata$PathwayID, 5), "\n")
  cat("Sample node names in cnet_data:", head(cnet_data$name, 5), "\n")
  cat("Number of matches:", sum(cnet_data$name %in% pathway_metadata$PathwayID, na.rm = TRUE), "\n")
  
  # Reconstruct edge data
  node_coords <- data.frame(
    name = cnet_data$name,
    x = cnet_data$x,
    y = cnet_data$y,
    stringsAsFactors = FALSE
  )
  node_coords <- node_coords[!duplicated(node_coords$name), ]
  
  edges_list <- list()
  edge_idx <- 1
  for (gene in names(filtered_gene_pathway_list)) {
    pathways <- filtered_gene_pathway_list[[gene]]
    gene_coords <- node_coords[node_coords$name == gene, c("x", "y")]
    if (nrow(gene_coords) == 0) {
      cat("Warning: No coordinates found for gene:", gene, "\n")
      next
    }
    
    for (pathway in pathways) {
      pathway_coords <- node_coords[node_coords$name == pathway, c("x", "y")]
      if (nrow(pathway_coords) == 0) {
        cat("Warning: No coordinates found for pathway:", pathway, "\n")
        next
      }
      
      edges_list[[edge_idx]] <- data.frame(
        x = gene_coords$x,
        y = gene_coords$y,
        xend = pathway_coords$x,
        yend = pathway_coords$y,
        stringsAsFactors = FALSE
      )
      edge_idx <- edge_idx + 1
    }
  }
  
  cnet_edges <- do.call(rbind, edges_list)
  
  # Debug: Check reconstructed edge data
  cat("Number of edges reconstructed:", if (is.null(cnet_edges)) 0 else nrow(cnet_edges), "\n")
  if (!is.null(cnet_edges)) {
    cat("Sample edge data:\n")
    print(head(cnet_edges, 5))
  }
  
  # Ensure PathwayID and cnet_data$name are in the same format
  pathway_metadata$PathwayID <- as.character(trimws(pathway_metadata$PathwayID))
  cnet_data$name <- as.character(trimws(cnet_data$name))
  
  # Create a mapping of node names to shapes
  node_shapes <- rep("circle", nrow(cnet_data))
  is_pathway <- cnet_data$name %in% pathway_metadata$PathwayID
  matched_indices <- match(cnet_data$name[is_pathway], pathway_metadata$PathwayID)
  matched_gene_set_types <- pathway_metadata$GeneSetType[matched_indices]
  node_shapes[is_pathway] <- ifelse(matched_gene_set_types == "GO", "diamond",
                                    ifelse(matched_gene_set_types == "KEGG", "triangle", "square"))
  
  cnet_data$shape <- node_shapes
  
  # Assign colors based on Direction
  cnet_data$color <- ifelse(cnet_data$name %in% names(gene_direction),
                            ifelse(sapply(gene_direction[cnet_data$name], function(x) x$direction) == "Up", "red",
                                   ifelse(sapply(gene_direction[cnet_data$name], function(x) x$direction) == "Down", "green", "purple")),
                            "grey")
  
  # Adjust node sizes based on ConditionCount
  cnet_data$size <- ifelse(cnet_data$name %in% pathway_metadata$PathwayID,
                           3 + pathway_metadata$ConditionCount[match(cnet_data$name, pathway_metadata$PathwayID)] * 2,
                           3)
  
  # Rebuild the plot with reconstructed edges (with labels)
  cnet_custom <- ggplot()
  if (!is.null(cnet_edges) && nrow(cnet_edges) > 0) {
    cnet_custom <- cnet_custom +
      geom_segment(data = cnet_edges, aes(x = x, y = y, xend = xend, yend = yend), color = "grey", alpha = 0.5)
  } else {
    cat("Warning: No edges reconstructed for", contrast_name, ". Plotting nodes only.\n")
  }
  cnet_custom <- cnet_custom +
    geom_point(data = cnet_data, aes(x = x, y = y, size = size, shape = shape, color = color)) +
    scale_shape_manual(values = c("circle" = 16, "diamond" = 18, "triangle" = 17, "square" = 15)) +
    scale_color_identity() +
    scale_size_continuous(range = c(3, 15)) +
    geom_text_repel(data = cnet_data, 
                    aes(x = x, y = y, label = name),
                    color = "black",
                    size = 3,
                    box.padding = 0.5,
                    max.overlaps = 20,
                    segment.color = "grey50",
                    segment.size = 0.2) +
    theme_void() +
    ggtitle(paste("Cnet Plot: Pathways for", contrast_name)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save the cnet plot (with labels)
  ggsave(file.path(output_dir, paste0("cnet_plot_", contrast_name, ".png")), cnet_custom, width = 20, height = 15, dpi = 300)
  
  cat("Cnet plot saved as", file.path(output_dir, paste0("cnet_plot_", contrast_name, ".png")), "\n")
  
  # Rebuild the cnet plot without labels
  cnet_custom_no_labels <- ggplot()
  if (!is.null(cnet_edges) && nrow(cnet_edges) > 0) {
    cnet_custom_no_labels <- cnet_custom_no_labels +
      geom_segment(data = cnet_edges, aes(x = x, y = y, xend = xend, yend = yend), color = "grey", alpha = 0.5)
  } else {
    cat("Warning: No edges reconstructed for", contrast_name, ". Plotting nodes only.\n")
  }
  cnet_custom_no_labels <- cnet_custom_no_labels +
    geom_point(data = cnet_data, aes(x = x, y = y, size = size, shape = shape, color = color)) +
    scale_shape_manual(values = c("circle" = 16, "diamond" = 18, "triangle" = 17, "square" = 15)) +
    scale_color_identity() +
    scale_size_continuous(range = c(3, 15)) +
    theme_void() +
    ggtitle(paste("Cnet Plot: Pathways for", contrast_name)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save the cnet plot (without labels)
  ggsave(file.path(output_dir, paste0("cnet_plot_no_labels_", contrast_name, ".png")), cnet_custom_no_labels, width = 20, height = 15, dpi = 300)
  
  cat("Cnet plot without labels saved as", file.path(output_dir, paste0("cnet_plot_no_labels_", contrast_name, ".png")), "\n")
  
  # Save the merged enrichment results as a CSV
  write.csv(merged_results,
            file = file.path(output_dir, paste0("cnet_merged_pathways_", contrast_name, ".csv")),
            quote = FALSE,
            row.names = FALSE)
  
  cat("Merged enrichment results saved as", file.path(output_dir, paste0("cnet_merged_pathways_", contrast_name, ".csv")), "\n")
}







##################################
###Epistasis

'
# Install and load required libraries
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  BiocManager::install("msigdbr")
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
'


library(VennDiagram)
library(ggVennDiagram)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)

# Configuration parameters
config <- list(
  padj_threshold = 0.05,
  lfc_threshold = 1,
  pvalue_cutoff = 0.1
)

# Define comparisons for treatment vs DMSO
comparisons <- list(
  NT = "NT",
  Compound_1_100nM = "Compound-1_100nM",
  Compound_3_100nM = "Compound-3_100nM",
  Compound_4_100nM = "Compound-4_100nM",
  Compound_7_100nM = "Compound-7_100nM",
  Compound_2_100nM = "Compound-2_100nM",
  Compound_5_1uM = "Compound-5_1uM",
  Compound_6_1uM = "Compound-6_1uM"
)

# Create output directories
venn_output_dir <- "Venn_Diagrams_Treatment_vs_DMSO"
epistasis_output_dir <- "Epistasis_Results_Treatment_vs_DMSO"
boxplot_output_dir <- "Boxplot_Analysis_Treatment_vs_DMSO"
if (!dir.exists(venn_output_dir)) dir.create(venn_output_dir)
if (!dir.exists(epistasis_output_dir)) dir.create(epistasis_output_dir)
if (!dir.exists(boxplot_output_dir)) dir.create(boxplot_output_dir)

# Step 1: Extract up- and downregulated genes for each condition
gene_lists_up <- list()
gene_lists_down <- list()

for (cond_name in names(comparisons)) {
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond_name, "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("Error: DESeq2 results file does not exist:", deseq_results_file, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, stringsAsFactors = FALSE)
  if (!"gene_name" %in% colnames(res)) {
    cat("Error: 'gene_name' column missing in", deseq_results_file, "- Skipping\n")
    next
  }
  
  # Identify problematic gene names (e.g., '1-Mar', '2-Mar') and duplicates
  problematic_genes <- res$gene_name[res$gene_name %in% c("1-Mar", "2-Mar")]
  duplicate_genes <- res$gene_name[duplicated(res$gene_name) | duplicated(res$gene_name, fromLast = TRUE)]
  genes_to_remove <- unique(c(problematic_genes, duplicate_genes))
  
  if (length(genes_to_remove) > 0) {
    cat("Removing", length(genes_to_remove), "problematic/duplicate genes from", deseq_results_file, ":\n")
    print(genes_to_remove)
    res <- res[!res$gene_name %in% genes_to_remove, ]
  }
  
  # Ensure gene_name is unique
  if (any(duplicated(res$gene_name))) {
    cat("Warning: Additional duplicates found in", deseq_results_file, "after initial filtering. Skipping condition.\n")
    next
  }
  
  rownames(res) <- res$gene_name
  res$log2FoldChange <- as.numeric(res$log2FoldChange)
  res$padj <- as.numeric(res$padj)
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- res_df$gene_name
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  # Extract up- and downregulated genes
  upregulated_genes <- res_df$Symbol[res_df$log2FoldChange > config$lfc_threshold & res_df$padj < config$padj_threshold]
  downregulated_genes <- res_df$Symbol[res_df$log2FoldChange < -config$lfc_threshold & res_df$padj < config$padj_threshold]
  
  gene_lists_up[[cond_name]] <- unique(upregulated_genes)
  gene_lists_down[[cond_name]] <- unique(downregulated_genes)
  
  cat("Condition:", cond_name, "- Upregulated genes:", length(upregulated_genes), 
      "Downregulated genes:", length(downregulated_genes), "\n")
}

# Step 2: Generate Venn diagrams for pairwise comparisons using ggVennDiagram
conditions <- names(comparisons)

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    
    # Upregulated genes Venn diagram
    up_genes1 <- gene_lists_up[[cond1]]
    up_genes2 <- gene_lists_up[[cond2]]
    
    if (length(up_genes1) > 0 && length(up_genes2) > 0) {
      tryCatch({
        venn_data <- list(up_genes1, up_genes2)
        names(venn_data) <- c(cond1, cond2)
        
        venn_plot_up <- ggVennDiagram(venn_data, 
                                      label_alpha = 0.5, 
                                      set_color = c("#FF9999", "#99CCFF"),
                                      set_size = 5) +
          ggtitle(paste("Upregulated Genes:", cond1, "vs", cond2)) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        
        ggsave(file.path(venn_output_dir, paste0("venn_upregulated_", cond1, "_vs_", cond2, ".png")), 
               venn_plot_up, width = 6, height = 6, bg = "white")
        
        common_up_genes <- intersect(up_genes1, up_genes2)
        cat("Common upregulated genes between", cond1, "and", cond2, ":", length(common_up_genes), "\n")
        
        if (length(common_up_genes) > 0) {
          write.table(data.frame(Gene = common_up_genes),
                      file = file.path(epistasis_output_dir, paste0("common_upregulated_", cond1, "_vs_", cond2, ".tsv")),
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }, error = function(e) {
        cat("Error generating upregulated Venn diagram for", cond1, "vs", cond2, ":", e$message, "\n")
      })
    }
    
    # Downregulated genes Venn diagram
    down_genes1 <- gene_lists_down[[cond1]]
    down_genes2 <- gene_lists_down[[cond2]]
    
    if (length(down_genes1) > 0 && length(down_genes2) > 0) {
      tryCatch({
        venn_data <- list(down_genes1, down_genes2)
        names(venn_data) <- c(cond1, cond2)
        
        venn_plot_down <- ggVennDiagram(venn_data, 
                                        label_alpha = 0.5, 
                                        set_color = c("#FF9999", "#99CCFF"),
                                        set_size = 5) +
          ggtitle(paste("Downregulated Genes:", cond1, "vs", cond2)) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        
        ggsave(file.path(venn_output_dir, paste0("venn_downregulated_", cond1, "_vs_", cond2, ".png")), 
               venn_plot_down, width = 6, height = 6, bg = "white")
        
        common_down_genes <- intersect(down_genes1, down_genes2)
        cat("Common downregulated genes between", cond1, "and", cond2, ":", length(common_down_genes), "\n")
        
        if (length(common_down_genes) > 0) {
          write.table(data.frame(Gene = common_down_genes),
                      file = file.path(epistasis_output_dir, paste0("common_downregulated_", cond1, "_vs_", cond2, ".tsv")),
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }, error = function(e) {
        cat("Error generating downregulated Venn diagram for", cond1, "vs", cond2, ":", e$message, "\n")
      })
    }
  }
}

# Step 3: Heatmap of overlap matrix
n_conditions <- length(conditions)
up_overlap_matrix <- matrix(0, nrow = n_conditions, ncol = n_conditions)
down_overlap_matrix <- matrix(0, nrow = n_conditions, ncol = n_conditions)
rownames(up_overlap_matrix) <- conditions
colnames(up_overlap_matrix) <- conditions
rownames(down_overlap_matrix) <- conditions
colnames(down_overlap_matrix) <- conditions

for (i in 1:(n_conditions - 1)) {
  for (j in (i + 1):n_conditions) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    
    # Upregulated genes overlap
    up_genes1 <- gene_lists_up[[cond1]]
    up_genes2 <- gene_lists_up[[cond2]]
    
    if (length(up_genes1) > 0 && length(up_genes2) > 0) {
      common_up_genes <- intersect(up_genes1, up_genes2)
      total_unique_up <- length(unique(c(up_genes1, up_genes2)))
      if (total_unique_up > 0) {
        percent_overlap_up <- (length(common_up_genes) / total_unique_up) * 100
        up_overlap_matrix[cond1, cond2] <- percent_overlap_up
        up_overlap_matrix[cond2, cond1] <- percent_overlap_up
      }
      
      cat("Percentage overlap (upregulated) between", cond1, "and", cond2, ":", percent_overlap_up, "%\n")
    }
    
    # Downregulated genes overlap
    down_genes1 <- gene_lists_down[[cond1]]
    down_genes2 <- gene_lists_down[[cond2]]
    
    if (length(down_genes1) > 0 && length(down_genes2) > 0) {
      common_down_genes <- intersect(down_genes1, down_genes2)
      total_unique_down <- length(unique(c(down_genes1, down_genes2)))
      if (total_unique_down > 0) {
        percent_overlap_down <- (length(common_down_genes) / total_unique_down) * 100
        down_overlap_matrix[cond1, cond2] <- percent_overlap_down
        down_overlap_matrix[cond2, cond1] <- percent_overlap_down
      }
      
      cat("Percentage overlap (downregulated) between", cond1, "and", cond2, ":", percent_overlap_down, "%\n")
    }
  }
}

# Generate heatmap for upregulated genes overlap
heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
pheatmap(up_overlap_matrix,
         color = heatmap_colors,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Percentage Overlap of Upregulated Genes (Treatment vs DMSO)",
         fontsize = 8,
         display_numbers = TRUE,
         number_format = "%.1f",
         filename = file.path(epistasis_output_dir, "overlap_heatmap_upregulated_treatment_vs_dmso.png"),
         width = 10,
         height = 8)

write.csv(up_overlap_matrix,
          file = file.path(epistasis_output_dir, "overlap_matrix_upregulated_treatment_vs_dmso.csv"),
          quote = FALSE)


#####Special!!!
down_overlap_matrix_jitter <- jitter(as.matrix(down_overlap_matrix), amount = 0.01)
pheatmap(down_overlap_matrix_jitter,
         color = heatmap_colors,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Percentage Overlap of Downregulated Genes (Treatment vs DMSO)",
         fontsize = 8,
         display_numbers = TRUE,
         number_format = "%.1f",
         filename = file.path(epistasis_output_dir, "overlap_heatmap_downregulated_treatment_vs_dmso.png"),
         width = 10,
         height = 8)
'
# Generate heatmap for downregulated genes overlap
pheatmap(down_overlap_matrix,
         color = heatmap_colors,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Percentage Overlap of Downregulated Genes (Treatment vs DMSO)",
         fontsize = 8,
         display_numbers = TRUE,
         number_format = "%.1f",
         filename = file.path(epistasis_output_dir, "overlap_heatmap_downregulated_treatment_vs_dmso.png"),
         width = 10,
         height = 8)
'
write.csv(down_overlap_matrix,
          file = file.path(epistasis_output_dir, "overlap_matrix_downregulated_treatment_vs_dmso.csv"),
          quote = FALSE)

cat("Overlap heatmaps and data saved in", epistasis_output_dir, "\n")

# Step 4: Measure epistasis for common genes and visualize as heatmaps
up_epistasis_data <- list()
down_epistasis_data <- list()
all_common_up_genes <- c()
all_common_down_genes <- c()
condition_pairs <- c()

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    
    deseq_file1 <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond1, "_vs_DMSO_deseq_results.csv"))
    deseq_file2 <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond2, "_vs_DMSO_deseq_results.csv"))
    
    if (!file.exists(deseq_file1) || !file.exists(deseq_file2)) {
      cat("Skipping epistasis analysis for", cond1, "vs", cond2, "- DESeq2 files missing\n")
      next
    }
    
    res1 <- read.csv(deseq_file1, stringsAsFactors = FALSE)
    if (!"gene_name" %in% colnames(res1)) stop("Error: 'gene_name' column missing in ", deseq_file1)
    
    # Remove problematic genes
    problematic_genes1 <- res1$gene_name[res1$gene_name %in% c("1-Mar", "2-Mar")]
    duplicate_genes1 <- res1$gene_name[duplicated(res1$gene_name) | duplicated(res1$gene_name, fromLast = TRUE)]
    genes_to_remove1 <- unique(c(problematic_genes1, duplicate_genes1))
    if (length(genes_to_remove1) > 0) {
      cat("Removing", length(genes_to_remove1), "problematic/duplicate genes from", deseq_file1, ":\n")
      print(genes_to_remove1)
      res1 <- res1[!res1$gene_name %in% genes_to_remove1, ]
    }
    
    rownames(res1) <- res1$gene_name
    res1$log2FoldChange <- as.numeric(res1$log2FoldChange)
    res1$Symbol <- res1$gene_name
    
    res2 <- read.csv(deseq_file2, stringsAsFactors = FALSE)
    if (!"gene_name" %in% colnames(res2)) stop("Error: 'gene_name' column missing in ", deseq_file2)
    
    # Remove problematic genes
    problematic_genes2 <- res2$gene_name[res2$gene_name %in% c("1-Mar", "2-Mar")]
    duplicate_genes2 <- res2$gene_name[duplicated(res2$gene_name) | duplicated(res2$gene_name, fromLast = TRUE)]
    genes_to_remove2 <- unique(c(problematic_genes2, duplicate_genes2))
    if (length(genes_to_remove2) > 0) {
      cat("Removing", length(genes_to_remove2), "problematic/duplicate genes from", deseq_file2, ":\n")
      print(genes_to_remove2)
      res2 <- res2[!res2$gene_name %in% genes_to_remove2, ]
    }
    
    rownames(res2) <- res2$gene_name
    res2$log2FoldChange <- as.numeric(res2$log2FoldChange)
    res2$Symbol <- res2$gene_name
    
    common_up_genes <- intersect(gene_lists_up[[cond1]], gene_lists_up[[cond2]])
    common_down_genes <- intersect(gene_lists_down[[cond1]], gene_lists_down[[cond2]])
    
    pair_label <- paste(cond1, "_vs_", cond2, sep = "")
    condition_pairs <- c(condition_pairs, pair_label)
    
    if (length(common_up_genes) > 0) {
      lfc1_up <- res1$log2FoldChange[match(common_up_genes, res1$Symbol)]
      lfc2_up <- res2$log2FoldChange[match(common_up_genes, res2$Symbol)]
      
      valid_idx_up <- !is.na(lfc1_up) & !is.na(lfc2_up)
      if (sum(valid_idx_up) == 0) {
        cat("No valid LFC values for upregulated genes between", cond1, "and", cond2, "- Skipping\n")
        next
      }
      
      common_up_genes <- common_up_genes[valid_idx_up]
      lfc1_up <- lfc1_up[valid_idx_up]
      lfc2_up <- lfc2_up[valid_idx_up]
      
      expected_lfc_up <- lfc1_up + lfc2_up
      epistasis_score_up <- abs(lfc1_up) + abs(lfc2_up) - abs(expected_lfc_up)
      
      epistasis_up <- data.frame(
        Gene = common_up_genes,
        LFC_Cond1 = lfc1_up,
        LFC_Cond2 = lfc2_up,
        Expected_LFC = expected_lfc_up,
        Epistasis_Score = epistasis_score_up
      )
      
      up_epistasis_data[[pair_label]] <- setNames(epistasis_score_up, common_up_genes)
      all_common_up_genes <- unique(c(all_common_up_genes, common_up_genes))
      
      tryCatch({
        write.table(epistasis_up,
                    file = file.path(epistasis_output_dir, paste0("epistasis_upregulated_", cond1, "_vs_", cond2, ".tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
      }, error = function(e) {
        cat("Error saving epistasis results for upregulated genes", cond1, "vs", cond2, ":", e$message, "\n")
      })
    }
    
    if (length(common_down_genes) > 0) {
      lfc1_down <- res1$log2FoldChange[match(common_down_genes, res1$Symbol)]
      lfc2_down <- res2$log2FoldChange[match(common_down_genes, res2$Symbol)]
      
      valid_idx_down <- !is.na(lfc1_down) & !is.na(lfc2_down)
      if (sum(valid_idx_down) == 0) {
        cat("No valid LFC values for downregulated genes between", cond1, "and", cond2, "- Skipping\n")
        next
      }
      
      common_down_genes <- common_down_genes[valid_idx_down]
      lfc1_down <- lfc1_down[valid_idx_down]
      lfc2_down <- lfc2_down[valid_idx_down]
      
      expected_lfc_down <- lfc1_down + lfc2_down
      epistasis_score_down <- abs(lfc1_down) + abs(lfc2_down) - abs(expected_lfc_down)
      
      epistasis_down <- data.frame(
        Gene = common_down_genes,
        LFC_Cond1 = lfc1_down,
        LFC_Cond2 = lfc2_down,
        Expected_LFC = expected_lfc_down,
        Epistasis_Score = epistasis_score_down
      )
      
      down_epistasis_data[[pair_label]] <- setNames(epistasis_score_down, common_down_genes)
      all_common_down_genes <- unique(c(all_common_down_genes, common_down_genes))
      
      tryCatch({
        write.table(epistasis_down,
                    file = file.path(epistasis_output_dir, paste0("epistasis_downregulated_", cond1, "_vs_", cond2, ".tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
      }, error = function(e) {
        cat("Error saving epistasis results for downregulated genes", cond1, "vs", cond2, ":", e$message, "\n")
      })
    }
  }
}

# Generate heatmap for upregulated genes
if (length(all_common_up_genes) > 0 && length(condition_pairs) > 0) {
  up_epistasis_matrix <- matrix(NA, nrow = length(all_common_up_genes), ncol = length(condition_pairs))
  rownames(up_epistasis_matrix) <- all_common_up_genes
  colnames(up_epistasis_matrix) <- condition_pairs
  
  for (pair in condition_pairs) {
    scores <- up_epistasis_data[[pair]]
    if (!is.null(scores)) {
      up_epistasis_matrix[names(scores), pair] <- scores
    }
  }
  
  up_epistasis_matrix[is.na(up_epistasis_matrix)] <- 0
  
  heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  score_range <- range(up_epistasis_matrix, na.rm = TRUE)
  if (diff(score_range) == 0) {
    cat("Skipping heatmap for upregulated genes: No variation in epistasis scores (all values are", score_range[1], ")\n")
  } else {
    breaks <- seq(score_range[1], score_range[2], length.out = 101)
    tryCatch({
      pheatmap(up_epistasis_matrix,
               color = heatmap_colors,
               breaks = breaks,
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = length(all_common_up_genes) <= 50,
               show_colnames = TRUE,
               main = "Epistasis Scores for Upregulated Genes (Treatment vs DMSO)",
               fontsize = 10,
               filename = file.path(epistasis_output_dir, "epistasis_heatmap_upregulated_treatment_vs_dmso.png"),
               width = 10,
               height = max(8, length(all_common_up_genes) * 0.1))
    }, error = function(e) {
      cat("Error in clustering for upregulated genes heatmap:", e$message, "\n")
      cat("Retrying without clustering...\n")
      pheatmap(up_epistasis_matrix,
               color = heatmap_colors,
               breaks = breaks,
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               show_rownames = length(all_common_up_genes) <= 50,
               show_colnames = TRUE,
               main = "Epistasis Scores for Upregulated Genes (Treatment vs DMSO, No Clustering)",
               fontsize = 10,
               filename = file.path(epistasis_output_dir, "epistasis_heatmap_upregulated_treatment_vs_dmso.png"),
               width = 10,
               height = max(8, length(all_common_up_genes) * 0.1))
    })
  }
}

# Generate heatmap for downregulated genes
if (length(all_common_down_genes) > 0 && length(condition_pairs) > 0) {
  down_epistasis_matrix <- matrix(NA, nrow = length(all_common_down_genes), ncol = length(condition_pairs))
  rownames(down_epistasis_matrix) <- all_common_down_genes
  colnames(down_epistasis_matrix) <- condition_pairs
  
  for (pair in condition_pairs) {
    scores <- down_epistasis_data[[pair]]
    if (!is.null(scores)) {
      down_epistasis_matrix[names(scores), pair] <- scores
    }
  }
  
  down_epistasis_matrix[is.na(down_epistasis_matrix)] <- 0
  
  heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  score_range <- range(down_epistasis_matrix, na.rm = TRUE)
  if (diff(score_range) == 0) {
    cat("Skipping heatmap for downregulated genes: No variation in epistasis scores (all values are", score_range[1], ")\n")
  } else {
    breaks <- seq(score_range[1], score_range[2], length.out = 101)
    tryCatch({
      pheatmap(down_epistasis_matrix,
               color = heatmap_colors,
               breaks = breaks,
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = length(all_common_down_genes) <= 50,
               show_colnames = TRUE,
               main = "Epistasis Scores for Downregulated Genes (Treatment vs DMSO)",
               fontsize = 10,
               filename = file.path(epistasis_output_dir, "epistasis_heatmap_downregulated_treatment_vs_dmso.png"),
               width = 10,
               height = max(8, length(all_common_down_genes) * 0.1))
    }, error = function(e) {
      cat("Error in clustering for downregulated genes heatmap:", e$message, "\n")
      cat("Retrying without clustering...\n")
      pheatmap(down_epistasis_matrix,
               color = heatmap_colors,
               breaks = breaks,
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               show_rownames = length(all_common_down_genes) <= 50,
               show_colnames = TRUE,
               main = "Epistasis Scores for Downregulated Genes (Treatment vs DMSO, No Clustering)",
               fontsize = 10,
               filename = file.path(epistasis_output_dir, "epistasis_heatmap_downregulated_treatment_vs_dmso.png"),
               width = 10,
               height = max(8, length(all_common_down_genes) * 0.1))
    })
  }
}

# Step 5: Pathway epistasis - Combined pathway enrichment
all_gene_sets <- c()
enrichment_results <- list()

for (cond_name in conditions) {
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond_name, "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", cond_name, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, stringsAsFactors = FALSE)
  if (!"gene_name" %in% colnames(res)) {
    cat("Error: 'gene_name' column missing in", deseq_results_file, "- Skipping\n")
    next
  }
  
  # Remove problematic genes
  problematic_genes <- res$gene_name[res$gene_name %in% c("1-Mar", "2-Mar")]
  duplicate_genes <- res$gene_name[duplicated(res$gene_name) | duplicated(res$gene_name, fromLast = TRUE)]
  genes_to_remove <- unique(c(problematic_genes, duplicate_genes))
  if (length(genes_to_remove) > 0) {
    cat("Removing", length(genes_to_remove), "problematic/duplicate genes from", deseq_results_file, ":\n")
    print(genes_to_remove)
    res <- res[!res$gene_name %in% genes_to_remove, ]
  }
  
  rownames(res) <- res$gene_name
  res$log2FoldChange <- as.numeric(res$log2FoldChange)
  res$padj <- as.numeric(res$padj)
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- res_df$gene_name
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  up_genes <- res_df$Symbol[res_df$log2FoldChange > config$lfc_threshold & res_df$padj < config$padj_threshold]
  down_genes <- res_df$Symbol[res_df$log2FoldChange < -config$lfc_threshold & res_df$padj < config$padj_threshold]
  
  entrez_ids_up <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  entrez_ids_down <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  ranked_genes <- res_df$log2FoldChange
  names(ranked_genes) <- res_df$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  msigdb <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
  all_genes_in_sets <- unique(gene_sets_hallmark$gene)
  
  ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
  if (length(ranked_genes_filtered) > 0) {
    gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                         stats = ranked_genes_filtered,
                         minSize = 15,
                         maxSize = 500,
                         nPermSimple = 1000)
    
    up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < config$padj_threshold, ]
    down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < config$padj_threshold, ]
    
    if (nrow(up_enriched_hallmark) > 0) {
      up_enriched_hallmark$pathway <- paste0("Hallmark_", up_enriched_hallmark$pathway)
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        up_enriched_hallmark[, c("pathway", "padj")]
      )
      all_gene_sets <- unique(c(all_gene_sets, up_enriched_hallmark$pathway))
    }
    if (nrow(down_enriched_hallmark) > 0) {
      down_enriched_hallmark$pathway <- paste0("Hallmark_", down_enriched_hallmark$pathway)
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        down_enriched_hallmark[, c("pathway", "padj")]
      )
      all_gene_sets <- unique(c(all_gene_sets, down_enriched_hallmark$pathway))
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_ids_up$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      kegg_up_df <- as.data.frame(kegg_up)
      kegg_up_df$pathway <- paste0("KEGG_", kegg_up_df$Description)
      kegg_up_df <- kegg_up_df[, c("pathway", "p.adjust")]
      colnames(kegg_up_df) <- c("pathway", "padj")
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        kegg_up_df
      )
      all_gene_sets <- unique(c(all_gene_sets, kegg_up_df$pathway))
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_ids_down$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      kegg_down_df <- as.data.frame(kegg_down)
      kegg_down_df$pathway <- paste0("KEGG_", kegg_down_df$Description)
      kegg_down_df <- kegg_down_df[, c("pathway", "p.adjust")]
      colnames(kegg_down_df) <- c("pathway", "padj")
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        kegg_down_df
      )
      all_gene_sets <- unique(c(all_gene_sets, kegg_down_df$pathway))
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    go_up <- tryCatch({
      enrichGO(gene = entrez_ids_up$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      go_up_df <- as.data.frame(go_up)
      go_up_df$pathway <- paste0("GO_BP_", go_up_df$Description)
      go_up_df <- go_up_df[, c("pathway", "p.adjust")]
      colnames(go_up_df) <- c("pathway", "padj")
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        go_up_df
      )
      all_gene_sets <- unique(c(all_gene_sets, go_up_df$pathway))
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    go_down <- tryCatch({
      enrichGO(gene = entrez_ids_down$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      go_down_df <- as.data.frame(go_down)
      go_down_df$pathway <- paste0("GO_BP_", go_down_df$Description)
      go_down_df <- go_down_df[, c("pathway", "p.adjust")]
      colnames(go_down_df) <- c("pathway", "padj")
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        go_down_df
      )
      all_gene_sets <- unique(c(all_gene_sets, go_down_df$pathway))
    }
  }
  
  cat("Processed", cond_name, "- Total unique gene sets:", length(all_gene_sets), "\n")
}

# Construct heatmap matrix
if (length(all_gene_sets) == 0) {
  stop("No significant gene sets found across all conditions.")
}

condition_dirs <- names(enrichment_results)
heatmap_matrix <- matrix(NA, nrow = length(condition_dirs), ncol = length(all_gene_sets))
rownames(heatmap_matrix) <- condition_dirs
colnames(heatmap_matrix) <- all_gene_sets

for (cond_dir in condition_dirs) {
  enriched <- enrichment_results[[cond_dir]]
  if (!is.null(enriched)) {
    for (i in 1:nrow(enriched)) {
      pathway <- enriched$pathway[i]
      padj <- enriched$padj[i]
      heatmap_matrix[cond_dir, pathway] <- -log10(padj)
    }
  }
}

heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

tryCatch({
  pheatmap(heatmap_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = length(all_gene_sets) <= 50,
           main = "Combined Pathway Enrichment Across Treatment vs DMSO Conditions (-log10(padj))",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "combined_pathway_enrichment_heatmap_treatment_vs_dmso.png"),
           width = max(10, length(all_gene_sets) * 0.2),
           height = max(8, length(condition_dirs) * 0.2))
}, error = function(e) {
  cat("Error in clustering for pathway enrichment heatmap:", e$message, "\n")
  cat("Retrying without clustering...\n")
  pheatmap(heatmap_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = length(all_gene_sets) <= 50,
           main = "Combined Pathway Enrichment Across Treatment vs DMSO Conditions (-log10(padj), No Clustering)",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "combined_pathway_enrichment_heatmap_treatment_vs_dmso.png"),
           width = max(10, length(all_gene_sets) * 0.2),
           height = max(8, length(condition_dirs) * 0.2))
})

write.csv(heatmap_matrix,
          file = file.path(epistasis_output_dir, "combined_pathway_enrichment_heatmap_data_treatment_vs_dmso.csv"),
          quote = FALSE)

cat("Combined pathway enrichment heatmap and data saved in", epistasis_output_dir, "\n")

# Step 6: Pathway epistasis - Top 20 gene sets
enrichment_results <- list()
top_gene_sets <- list(Hallmark = c(), KEGG = c(), GO_BP = c())

for (cond_name in conditions) {
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond_name, "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", cond_name, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, stringsAsFactors = FALSE)
  if (!"gene_name" %in% colnames(res)) {
    cat("Error: 'gene_name' column missing in", deseq_results_file, "- Skipping\n")
    next
  }
  
  # Remove problematic genes
  problematic_genes <- res$gene_name[res$gene_name %in% c("1-Mar", "2-Mar")]
  duplicate_genes <- res$gene_name[duplicated(res$gene_name) | duplicated(res$gene_name, fromLast = TRUE)]
  genes_to_remove <- unique(c(problematic_genes, duplicate_genes))
  if (length(genes_to_remove) > 0) {
    cat("Removing", length(genes_to_remove), "problematic/duplicate genes from", deseq_results_file, ":\n")
    print(genes_to_remove)
    res <- res[!res$gene_name %in% genes_to_remove, ]
  }
  
  rownames(res) <- res$gene_name
  res$log2FoldChange <- as.numeric(res$log2FoldChange)
  res$padj <- as.numeric(res$padj)
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- res_df$gene_name
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  up_genes <- res_df$Symbol[res_df$log2FoldChange > config$lfc_threshold & res_df$padj < config$padj_threshold]
  down_genes <- res_df$Symbol[res_df$log2FoldChange < -config$lfc_threshold & res_df$padj < config$padj_threshold]
  
  entrez_ids_up <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  entrez_ids_down <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  ranked_genes <- res_df$log2FoldChange
  names(ranked_genes) <- res_df$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  msigdb <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
  all_genes_in_sets <- unique(gene_sets_hallmark$gene)
  
  ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
  if (length(ranked_genes_filtered) > 0) {
    gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                         stats = ranked_genes_filtered,
                         minSize = 15,
                         maxSize = 500,
                         nPermSimple = 1000)
    
    up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < config$padj_threshold, ]
    down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < config$padj_threshold, ]
    
    if (nrow(up_enriched_hallmark) > 0) {
      up_enriched_hallmark$pathway <- paste0("Hallmark_", up_enriched_hallmark$pathway)
      up_enriched_hallmark <- up_enriched_hallmark[order(up_enriched_hallmark$padj), ]
      top_gene_sets$Hallmark <- unique(c(top_gene_sets$Hallmark, head(up_enriched_hallmark$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        up_enriched_hallmark[, c("pathway", "padj")]
      )
    }
    if (nrow(down_enriched_hallmark) > 0) {
      down_enriched_hallmark$pathway <- paste0("Hallmark_", down_enriched_hallmark$pathway)
      down_enriched_hallmark <- down_enriched_hallmark[order(down_enriched_hallmark$padj), ]
      top_gene_sets$Hallmark <- unique(c(top_gene_sets$Hallmark, head(down_enriched_hallmark$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        down_enriched_hallmark[, c("pathway", "padj")]
      )
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_ids_up$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      kegg_up_df <- as.data.frame(kegg_up)
      kegg_up_df$pathway <- paste0("KEGG_", kegg_up_df$Description)
      kegg_up_df <- kegg_up_df[, c("pathway", "p.adjust")]
      colnames(kegg_up_df) <- c("pathway", "padj")
      kegg_up_df <- kegg_up_df[order(kegg_up_df$padj), ]
      top_gene_sets$KEGG <- unique(c(top_gene_sets$KEGG, head(kegg_up_df$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        kegg_up_df
      )
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_ids_down$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      kegg_down_df <- as.data.frame(kegg_down)
      kegg_down_df$pathway <- paste0("KEGG_", kegg_down_df$Description)
      kegg_down_df <- kegg_down_df[, c("pathway", "p.adjust")]
      colnames(kegg_down_df) <- c("pathway", "padj")
      kegg_down_df <- kegg_down_df[order(kegg_down_df$padj), ]
      top_gene_sets$KEGG <- unique(c(top_gene_sets$KEGG, head(kegg_down_df$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        kegg_down_df
      )
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    go_up <- tryCatch({
      enrichGO(gene = entrez_ids_up$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      go_up_df <- as.data.frame(go_up)
      go_up_df$pathway <- paste0("GO_BP_", go_up_df$Description)
      go_up_df <- go_up_df[, c("pathway", "p.adjust")]
      colnames(go_up_df) <- c("pathway", "padj")
      go_up_df <- go_up_df[order(go_up_df$padj), ]
      top_gene_sets$GO_BP <- unique(c(top_gene_sets$GO_BP, head(go_up_df$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_up")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_up")]],
        go_up_df
      )
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    go_down <- tryCatch({
      enrichGO(gene = entrez_ids_down$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      go_down_df <- as.data.frame(go_down)
      go_down_df$pathway <- paste0("GO_BP_", go_down_df$Description)
      go_down_df <- go_down_df[, c("pathway", "p.adjust")]
      colnames(go_down_df) <- c("pathway", "padj")
      go_down_df <- go_down_df[order(go_down_df$padj), ]
      top_gene_sets$GO_BP <- unique(c(top_gene_sets$GO_BP, head(go_down_df$pathway, 20)))
      enrichment_results[[paste0(cond_name, "_down")]] <- rbind(
        enrichment_results[[paste0(cond_name, "_down")]],
        go_down_df
      )
    }
  }
  
  cat("Processed", cond_name, "- Total top gene sets (Hallmark:", length(top_gene_sets$Hallmark), 
      "KEGG:", length(top_gene_sets$KEGG), "GO_BP:", length(top_gene_sets$GO_BP), ")\n")
}

# Combine top gene sets
all_top_gene_sets <- unique(unlist(top_gene_sets))

# Construct heatmap matrix with top gene sets
if (length(all_top_gene_sets) == 0) {
  stop("No significant top gene sets found across all conditions.")
}

condition_dirs <- names(enrichment_results)
heatmap_matrix <- matrix(NA, nrow = length(condition_dirs), ncol = length(all_top_gene_sets))
rownames(heatmap_matrix) <- condition_dirs
colnames(heatmap_matrix) <- all_top_gene_sets

for (cond_dir in condition_dirs) {
  enriched <- enrichment_results[[cond_dir]]
  if (!is.null(enriched)) {
    for (i in 1:nrow(enriched)) {
      pathway <- enriched$pathway[i]
      if (pathway %in% all_top_gene_sets) {
        padj <- enriched$padj[i]
        heatmap_matrix[cond_dir, pathway] <- -log10(padj)
      }
    }
  }
}

heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

tryCatch({
  pheatmap(heatmap_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Top Pathway Enrichment Across Treatment vs DMSO Conditions (-log10(padj))",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "top20_pathway_enrichment_heatmap_treatment_vs_dmso.png"),
           width = max(10, length(all_top_gene_sets) * 0.2),
           height = max(8, length(condition_dirs) * 0.2))
}, error = function(e) {
  cat("Error in clustering for top pathway enrichment heatmap:", e$message, "\n")
  cat("Retrying without clustering...\n")
  pheatmap(heatmap_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Top Pathway Enrichment Across Treatment vs DMSO Conditions (-log10(padj), No Clustering)",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "top20_pathway_enrichment_heatmap_treatment_vs_dmso.png"),
           width = max(10, length(all_top_gene_sets) * 0.2),
           height = max(8, length(condition_dirs) * 0.2))
})

write.csv(heatmap_matrix,
          file = file.path(epistasis_output_dir, "top20_pathway_enrichment_heatmap_data_treatment_vs_dmso.csv"),
          quote = FALSE)

cat("Top 20 pathway enrichment heatmap and data saved in", epistasis_output_dir, "\n")

# Step 7: Pathway epistasis correlation plot/hypergeometric p-values
pathway_sets <- list()

for (cond_name in conditions) {
  deseq_results_file <- file.path("DESeq_Results_Cell_line_Genotype_1_DMSO", paste0(cond_name, "_vs_DMSO_deseq_results.csv"))
  if (!file.exists(deseq_results_file)) {
    cat("DESeq2 results file does not exist for", cond_name, "- Skipping\n")
    next
  }
  
  res <- read.csv(deseq_results_file, stringsAsFactors = FALSE)
  if (!"gene_name" %in% colnames(res)) {
    cat("Error: 'gene_name' column missing in", deseq_results_file, "- Skipping\n")
    next
  }
  
  # Remove problematic genes
  problematic_genes <- res$gene_name[res$gene_name %in% c("1-Mar", "2-Mar")]
  duplicate_genes <- res$gene_name[duplicated(res$gene_name) | duplicated(res$gene_name, fromLast = TRUE)]
  genes_to_remove <- unique(c(problematic_genes, duplicate_genes))
  if (length(genes_to_remove) > 0) {
    cat("Removing", length(genes_to_remove), "problematic/duplicate genes from", deseq_results_file, ":\n")
    print(genes_to_remove)
    res <- res[!res$gene_name %in% genes_to_remove, ]
  }
  
  rownames(res) <- res$gene_name
  res$log2FoldChange <- as.numeric(res$log2FoldChange)
  res$padj <- as.numeric(res$padj)
  res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  res_df <- as.data.frame(res)
  res_df$Symbol <- res_df$gene_name
  res_df <- res_df[!is.na(res_df$Symbol) & res_df$Symbol != "", ]
  res_df$Symbol <- trimws(res_df$Symbol)
  
  up_genes <- res_df$Symbol[res_df$log2FoldChange > config$lfc_threshold & res_df$padj < config$padj_threshold]
  down_genes <- res_df$Symbol[res_df$log2FoldChange < -config$lfc_threshold & res_df$padj < config$padj_threshold]
  
  entrez_ids_up <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping upregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  entrez_ids_down <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error mapping downregulated genes for", cond_name, ":", e$message, "\n")
    NULL
  })
  
  ranked_genes <- res_df$log2FoldChange
  names(ranked_genes) <- res_df$Symbol
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  condition_pathways <- c()
  
  msigdb <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets_hallmark <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol, stringsAsFactors = FALSE)
  all_genes_in_sets <- unique(gene_sets_hallmark$gene)
  
  ranked_genes_filtered <- ranked_genes[names(ranked_genes) %in% all_genes_in_sets]
  if (length(ranked_genes_filtered) > 0) {
    gsea_result <- fgsea(pathways = split(gene_sets_hallmark$gene, gene_sets_hallmark$term),
                         stats = ranked_genes_filtered,
                         minSize = 15,
                         maxSize = 500,
                         nPermSimple = 1000)
    
    up_enriched_hallmark <- gsea_result[gsea_result$NES > 0 & gsea_result$padj < config$padj_threshold, ]
    down_enriched_hallmark <- gsea_result[gsea_result$NES < 0 & gsea_result$padj < config$padj_threshold, ]
    
    if (nrow(up_enriched_hallmark) > 0) {
      condition_pathways <- c(condition_pathways, paste0("Hallmark_", up_enriched_hallmark$pathway))
    }
    if (nrow(down_enriched_hallmark) > 0) {
      condition_pathways <- c(condition_pathways, paste0("Hallmark_", down_enriched_hallmark$pathway))
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_ids_up$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      kegg_up_df <- as.data.frame(kegg_up)
      condition_pathways <- c(condition_pathways, paste0("KEGG_", kegg_up_df$Description))
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_ids_down$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = config$pvalue_cutoff,
                 pAdjustMethod = "fdr")
    }, error = function(e) {
      cat("KEGG enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      kegg_down_df <- as.data.frame(kegg_down)
      condition_pathways <- c(condition_pathways, paste0("KEGG_", kegg_down_df$Description))
    }
  }
  
  if (!is.null(entrez_ids_up) && nrow(entrez_ids_up) > 0) {
    go_up <- tryCatch({
      enrichGO(gene = entrez_ids_up$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for upregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
      go_up_df <- as.data.frame(go_up)
      condition_pathways <- c(condition_pathways, paste0("GO_BP_", go_up_df$Description))
    }
  }
  
  if (!is.null(entrez_ids_down) && nrow(entrez_ids_down) > 0) {
    go_down <- tryCatch({
      enrichGO(gene = entrez_ids_down$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pvalueCutoff = config$pvalue_cutoff,
               pAdjustMethod = "fdr",
               readable = FALSE)
    }, error = function(e) {
      cat("GO enrichment failed for downregulated genes in", cond_name, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
      go_down_df <- as.data.frame(go_down)
      condition_pathways <- c(condition_pathways, paste0("GO_BP_", go_down_df$Description))
    }
  }
  
  pathway_sets[[cond_name]] <- unique(condition_pathways)
  cat("Collected", length(condition_pathways), "pathways for", cond_name, "\n")
}

# Compute hypergeometric p-values for pathway overlap
valid_conditions <- names(pathway_sets)
n_conditions <- length(valid_conditions)
overlap_pval_matrix <- matrix(NA, nrow = n_conditions, ncol = n_conditions)
rownames(overlap_pval_matrix) <- valid_conditions
colnames(overlap_pval_matrix) <- valid_conditions

all_pathways <- unique(unlist(pathway_sets))
total_pathways <- length(all_pathways)

for (i in 1:n_conditions) {
  for (j in 1:n_conditions) {
    if (i == j) {
      overlap_pval_matrix[i, j] <- 0
      next
    }
    
    pathways1 <- pathway_sets[[valid_conditions[i]]]
    pathways2 <- pathway_sets[[valid_conditions[j]]]
    
    common_pathways <- length(intersect(pathways1, pathways2))
    size1 <- length(pathways1)
    size2 <- length(pathways2)
    
    pval <- phyper(common_pathways - 1, size2, total_pathways - size2, size1, lower.tail = FALSE)
    overlap_pval_matrix[i, j] <- -log10(pval)
    overlap_pval_matrix[j, i] <- -log10(pval)
  }
}

heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
overlap_pval_matrix[is.na(overlap_pval_matrix)] <- 0

tryCatch({
  pheatmap(overlap_pval_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Pathway Overlap Between Treatment vs DMSO Conditions (-log10(p-value))",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "pathway_overlap_heatmap_treatment_vs_dmso.png"),
           width = max(8, n_conditions * 0.5),
           height = max(8, n_conditions * 0.5))
}, error = function(e) {
  cat("Error in clustering for pathway overlap heatmap:", e$message, "\n")
  cat("Retrying without clustering...\n")
  pheatmap(overlap_pval_matrix,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Pathway Overlap Between Treatment vs DMSO Conditions (-log10(p-value), No Clustering)",
           fontsize = 8,
           display_numbers = FALSE,
           na_col = "grey",
           filename = file.path(epistasis_output_dir, "pathway_overlap_heatmap_treatment_vs_dmso.png"),
           width = max(8, n_conditions * 0.5),
           height = max(8, n_conditions * 0.5))
})

write.csv(overlap_pval_matrix,
          file = file.path(epistasis_output_dir, "pathway_overlap_heatmap_data_treatment_vs_dmso.csv"),
          quote = FALSE)

cat("Pathway overlap heatmap and data saved in", epistasis_output_dir, "\n")



#############################
#Boxplot for gene of interest
# Load count data
counts_df_gene_name <- read.csv("merged_counts_df_gene_name_filtered.csv")
if (!"gene_name" %in% colnames(counts_df_gene_name)) {
  stop("Error: 'gene_name' column missing in merged_counts_df_gene_name_filtered.csv")
}

# Prepare count matrix for DESeq2
counts_df_gene_name_wo_gn <- counts_df_gene_name[, !colnames(counts_df_gene_name) %in% c("gene_name"), drop = FALSE]

# Remove samples with all-zero counts
nonzero_samples <- colSums(counts_df_gene_name_wo_gn) > 0
counts_df_gene_name_wo_gn <- counts_df_gene_name_wo_gn[, nonzero_samples, drop = FALSE]

# Create colData
sample_cols <- colnames(counts_df_gene_name_wo_gn)
conditions <- rep(NA, length(sample_cols))

# Assign conditions
conditions[grepl("DMSO", sample_cols)] <- "DMSO"
drug_conditions <- c("Compound-1_100nM", "Compound-3_100nM", "Compound-4_100nM", "Compound-7_100nM", 
                     "Compound-2_100nM", "Compound-5_1uM", "Compound-6_1uM")
for (drug in drug_conditions) {
  conditions[grepl(gsub("-", "\\.", drug), sample_cols) & is.na(conditions)] <- drug
}
conditions[is.na(conditions)] <- "NT"

# Create colData
col_data <- data.frame(
  condition = factor(conditions, levels = c("DMSO", "NT", drug_conditions)),
  row.names = sample_cols,
  stringsAsFactors = FALSE
)

# Filter samples
col_data <- col_data[!is.na(col_data$condition), ]
counts_df_gene_name_wo_gn <- counts_df_gene_name_wo_gn[, rownames(col_data), drop = FALSE]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_df_gene_name_wo_gn),
                              colData = col_data,
                              design = ~ condition)

# Use poscounts for size factor estimation
dds <- estimateSizeFactors(dds, type = "poscounts")

# Apply VST normalization
vst_counts <- assay(vst(dds, blind = TRUE))

# Save VST-normalized counts
write.csv(vst_counts, file = file.path(boxplot_output_dir, "VST_normalized_counts.csv"), quote = FALSE)

# Define gene of interest
gene_of_interest <- "KEAP1"

# Validate gene
if (!(gene_of_interest %in% rownames(vst_counts))) {
  stop("Gene ", gene_of_interest, " not found in the count data.")
}

# Get VST-normalized expression values
gene_expression <- vst_counts[gene_of_interest, ]
expression_df <- data.frame(
  Sample = colnames(vst_counts),
  Expression = as.numeric(gene_expression),
  Condition = col_data$condition,
  stringsAsFactors = FALSE
)

# Clean sample names for plotting
expression_df$Sample_Clean <- sapply(expression_df$Sample, function(s) {
  condition <- "NT"
  if (grepl("DMSO", s)) {
    condition <- "DMSO"
  } else {
    for (drug in drug_conditions) {
      if (grepl(gsub("-", "\\.", drug), s)) {
        condition <- drug
        break
      }
    }
  }
  replicate <- sub(".*_([0-9])_24hr$", "\\1", s)
  paste(condition, replicate, sep = "_")
})

# Create boxplot
p <- ggplot(expression_df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.size = 1) +
  geom_point(aes(x = Condition, y = Expression),
             position = position_jitter(width = 0.15),
             size = 2, alpha = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "top") +
  labs(title = paste("Expression of", gene_of_interest, "(Treatment vs DMSO)"),
       x = "Treatment Condition",
       y = "VST-Normalized Expression") +
  scale_fill_brewer(palette = "Set3", name = "Condition")

# Save the plot
ggsave(file.path(boxplot_output_dir, paste0(tolower(gene_of_interest), "_expression_boxplot_treatment_vs_dmso.png")),
       plot = p, width = 12, height = 8, bg = "white")

cat("Boxplot for", gene_of_interest, "saved in", boxplot_output_dir, "\n")


