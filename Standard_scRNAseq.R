
setwd("/Users/az/Desktop/YouTubeTutorials-main/Outside/Own_scRNAseq")

# load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pheatmap)
#BiocManager::install("DESeq2")
library(DESeq2)
library(SingleR)

# Load the NSCLC dataset
hcc364.sparse <- Read10X_h5(filename = 'raw_feature_bc_matrix.h5',
                            use.names = TRUE,
                            unique.features = TRUE)
#str(hcc364.sparse)
hcc364.seurat_hdf5 <- CreateSeuratObject(counts = hcc364.sparse)
#saveRDS(hcc364.seurat_hdf5, '../hcc364.seurat_hdf5.rds')

hcc364.seurat_hdf5
str(hcc364.seurat_hdf5)
View(hcc364.seurat_hdf5@meta.data)

#some basic operations
Cells(hcc364.seurat_hdf5)
length (Cells(hcc364.seurat_hdf5))
row.names(hcc364.seurat_hdf5)
length((row.names(hcc364.seurat_hdf5)))
colnames(hcc364.seurat_hdf5)

# % MT reads
hcc364.seurat_hdf5$percent.mt <- PercentageFeatureSet(hcc364.seurat_hdf5, pattern = '^MT-')
View(hcc364.seurat_hdf5@meta.data)

#nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
#View(nsclc.seurat.obj@meta.data)

#DUsing ggplot
#ggplot(nsclc.seurat.obj@meta.data, aes(x= nsclc.seurat.obj@meta.data$nCount_RNA, y= nsclc.seurat.obj@meta.data$nFeature_RNA)) + geom_point()

#metadata_slot1 <- nsclc.seurat.obj@meta.data$nFeature_RNA
#metadata_slot2 <- nsclc.seurat.obj@meta.data$nCount_RNA
#scatterplot <-ggplot(data = pbmc.seurat@meta.data, aes(x = metadata_slot1, y = metadata_slot2)) +
  #geom_point()
## Label the x-axis
#scatter_plot <- scatter_plot + xlab("X-axis label")
#scatter_plot <- scatter_plot + ylab("Y-axis label")

#Replace the missing values
hcc364.seurat_hdf5@meta.data$nCount_RNA[is.na(nsclc.seurat.obj@meta.data$nCount_RNA)] <- 0
hcc364.seurat_hdf5@meta.data$nFeature_RNA[is.na(nsclc.seurat.obj@meta.data$nFeature_RNA)] <- 0
hcc364.seurat_hdf5@meta.data$percent.mt[is.na(nsclc.seurat.obj@meta.data$percent.mt)] <- 0
hcc364.seurat_hdf5@meta.data[is.na(hcc364.seurat_hdf5@meta.data)] <- 0


head(hcc364.seurat_hdf5@meta.data)
View(hcc364.seurat_hdf5@meta.data)

#Violinplot using ggplot
#ggplot(data = pbmc.seurat@meta.data, aes(x = metadata_feature, y = 1)) +
#geom_violin() +
  #ylab("") +  # Remove y-axis label
  #xlab("Metadata Feature")  # Set x-axis label


VlnPlot(hcc364.seurat_hdf5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(hcc364.seurat_hdf5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
hcc364.seurat_hdf5 <- subset(hcc364.seurat_hdf5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                          percent.mt < 5)

# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
hcc364.seurat_hdf5 <- NormalizeData(hcc364.seurat_hdf5)
str(hcc364.seurat_hdf5)

#after and before norm run this
# set seed and put two plots in one figure
#set.seed(123)
#par(mfrow=c(1,2))
# original expression distribution
raw_geneExp = as.vector(hcc364.seurat_hdf5[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
#
par(mar=c(5,6,4,1)+.1)
hist(raw_geneExp, breaks =100, xlab = "raw_geneExp", ylab = "frequency")


# 4. Identify highly variable features --------------
hcc364.seurat_hdf5 <- FindVariableFeatures(hcc364.seurat_hdf5)
hcc364.seurat_hdf5
#nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
VariableFeaturePlot(hcc364.seurat_hdf5)
VariableFeatures(hcc364.seurat_hdf5)
#head(VariableFeatures(hcc364.seurat_hdf5[['RNA']]))
head(VariableFeatures(hcc364.seurat_hdf5))
#tail(VariableFeatures(hcc364.seurat_hdf5[['RNA']]))
tail(VariableFeatures(hcc364.seurat_hdf5))
top10 <- head(VariableFeatures(hcc364.seurat_hdf5), 10)
top10
LabelPoints(plot = VariableFeaturePlot(hcc364.seurat_hdf5), points = top10, repel = TRUE)






# plot variable features with and without labels

#Does not work
#plot(top10)
#average_expression <- AverageExpression(object = nsclc.seurat.obj)
#DOES NOT WORK!!std_deviation <- apply(nsclc.seurat.obj@assays$RNA, 1, sd)
#install.packages("matrixStats")
#library(matrixStats)
#standard_deviation = rowVars(nsclc.seurat.obj@assays$RNA) #Error in rowVars(nsclc.seurat.obj@assays$RNA) : 
#Argument 'x' must be a matrix or a vector.


#count_matrix <- GetAssayData(object = nsclc.seurat.obj, assay = "RNA")
#std_deviation <- apply(count_matrix, 1, sd)


#std_deviation <- rowVars(assay(nsclc.seurat.obj, slot = "counts"))

#ggplot (top10, aes= (x= average_expression y= std_deviation)) + geom_point()


#plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
#LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Scaling -------------

#nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj)
#all.genes <- rownames(hcc364.seurat_hdf5)
#head(all.genes, 10)
#hcc364.seurat_hdf5 <- ScaleData(hcc364.seurat_hdf5, features = all.genes)
hcc364.seurat_hdf5 <- ScaleData(hcc364.seurat_hdf5)
str(hcc364.seurat_hdf5)

# 6. Perform Linear dimensionality reduction --------------
hcc364.seurat_hdf5 <- RunPCA(hcc364.seurat_hdf5, features = VariableFeatures(object = hcc364.seurat_hdf5))
str(hcc364.seurat_hdf5)
print(hcc364.seurat_hdf5[['pca']])
DimHeatmap(hcc364.seurat_hdf5, dims = 1:10, cells= 500, balanced= TRUE, fast = FALSE)
VizDimLoadings(hcc364.seurat_hdf5, dims = 1:2, reduction = "pca")

# visualize PCA results
#print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(nsclc.seurat.obj, dims = 1:5, cells = 500, balanced = TRUE)

Seurat::DimHeatmap(
  hcc364.seurat_hdf5,
  dims = 1:10,
  cells = 860,
  balanced = TRUE,
  fast = FALSE
) +
  ggplot2::scale_fill_gradientn(colors = c("blue", "yellow"))


# determine dimensionality of the data
ElbowPlot(hcc364.seurat_hdf5)


# 7. Clustering ------------
hcc364.seurat_hdf5 <- FindNeighbors(hcc364.seurat_hdf5, dims = 1:15)
hcc364.seurat_hdf5 <- FindClusters(hcc364.seurat_hdf5, resolution = c(0.1, 0.2, 0,.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
View(hcc364.seurat_hdf5@meta.data)
DimPlot(hcc364.seurat_hdf5, group.by = "RNA_snn_res.0.2", label = TRUE)
# understanding resolution
#nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
#View(nsclc.seurat.obj@meta.data)
#DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(hcc364.seurat_hdf5)
Idents(hcc364.seurat_hdf5) <- "RNA_snn_res.0.2"
Idents(hcc364.seurat_hdf5)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
hcc364.seurat_hdf5 <- RunUMAP(hcc364.seurat_hdf5, dims = 1:15)
hcc364.seurat_hdf5 <- RunTSNE(hcc364.seurat_hdf5, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
#DimPlot(hcc364.seurat_hd5, reduction = "umap")
#DimPlot(hcc364.seurat_hd5, reduction = "pca")
#DimPlot(hcc364.seurat_hd5, reduction = "tsne")


# individual clusters
DimPlot(hcc364.seurat_hdf5, reduction = "umap")
DimPlot(hcc364.seurat_hdf5, reduction = "pca")
DimPlot(hcc364.seurat_hdf5, reduction = "tsne")

#checking cluster id
cluster_ids <- hcc364.seurat_hdf5@meta.data$seurat_clusters
summary(cluster_ids)
cluster_freq <- (summary(cluster_ids))
type(cluster_freq)
cluster_freq[2,3]

cluster_freq_df <- data.frame(cluster_freq)
typeof(cluster_freq_df)
cluster_freq_df
cluster_freq_df <- rownames_to_column(cluster_freq_df, var = "cluster_id")
#colnames(cluster_freq_df)[colnames(cluster_freq_df) == ""] <- "cluster_id"
cluster_freq_df

barplot(summary(cluster_ids), xlab = "custer_id", ylab= "frequency", main = "custer_occupancy_by_cells")
ggplot(cluster_freq_df, aes(x= 'cluster_id', y= 'cluster_freq')) +
  geom_bar(stat = 'identity')

#Violinplot for the most variable genes on each clusters
VlnPlot(hcc364.seurat_hd5, features = top10)
hcc364_markers <- FindAllMarkers(hcc364.seurat_hd5)
write.csv(hcc364_markers, file = "../hcc364_marker_genes.csv", row.names = FALSE)

#Violinplot for the most variable genes on each clusters
VlnPlot(hcc364.seurat_hdf5, features = top10)
hcc364_markers_hdf5 <- FindAllMarkers(hcc364.seurat_hdf5) # log_FC= 0.25, min.pct = 0.1, only_pos = TRUE, test.use = 'DESeq2', slots = counts)
hcc364_markers_hdf5


str(hcc364.seurat_hdf5)
View(hcc364.seurat_hdf5@meta.data)


# Pairwise
#markers <- FindMarkers(
  #hcc364.seurat_hd5,
  #ident.1 = 0,  # Specify the cluster identity for comparison
  #ident.2 = 1,  # Specify another cluster identity for comparison
  #min.pct = 0.25,  # Minimum percentage of cells expressing the gene in any cluster
  #logfc.threshold = 0.25  # Minimum log-fold change for differential expression
#)

#cluster_markers <- hcc364_markers_hdf5@markers$seurat_cluster   
#marker_genes <- seurat_obj@markers$gene  


DefaultAssay(hcc364.seurat_hdf5)
# #find_conserved_markers_example <- FindConservedMarkers(hcc364_markers_hdf5, 
#                                           ident.3= 3, 
#                                           min.pct=  0.25,
#                                           logfc.threshold =0.25,
#                                           grouping.var= 'cluster')


#colnames(hcc364.seurat_hdf5)

FeaturePlot(hcc364.seurat_hdf5, features = c('AKAP9', "HES1"), min.cutoff = 'q10')
DotPlot(hcc364.seurat_hdf5, features = c("YAP1", "GAPDH", "LAMP1", "AKT1", "FGFR1"))



Idents(hcc364.seurat_hdf5)

#singleR automated cell marker
ref <- celldex::HumanPrimaryCellAtlasData()
ref_2 <-celldex::BlueprintEncodeData()

View(as.data.frame(colData(ref)))

hcc364_counts <- GetAssayData(hcc364.seurat_hdf5, slot = 'counts')
pred <- SingleR(test = hcc364_counts,
                ref = ref,
                labels = ref$label.main)

pred_2 <- SingleR(test = hcc364_counts,
                  ref = ref_2,
                  labels = ref_2$label.main)

pred
pred_2

print(pred$labels, 10)

View(hcc364.seurat_hdf5@meta.data)

hcc364.seurat_hdf5$singleR.labels <- pred$labels[match(rownames(hcc364.seurat_hdf5@meta.data), rownames(pred))]
DimPlot(hcc364.seurat_hdf5, reduction = 'umap', group.by = 'singleR.labels')

hcc364.seurat_hdf5$singleR_2.labels <- pred_2$labels[match(rownames(hcc364.seurat_hdf5@meta.data), rownames(pred_2))]
DimPlot(hcc364.seurat_hdf5, reduction = 'umap', group.by = 'singleR_2.labels')



# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores
pred_2$scores

plotScoreHeatmap(pred)
plotScoreHeatmap(pred_2)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)
plotDeltaDistribution(pred_2)
View(hcc364.seurat_hdf5@meta.data)
hcc364.seurat_hdf5@meta.data$cell_barcode <- rownames_to_column(hcc364.seurat_hdf5@meta.data)

#Pseudobulk analysis

cts_bulk <- AggregateExpression(hcc364.seurat_hdf5, 
                                group.by = 'seurat_clusters',
                                assays = "RNA",
                                slot = 'counts',
                                return.seurat = FALSE)
View(cts_bulk$RNA)
cts_bulk.t <-t(cts_bulk$RNA)
cts_bulk.t <- as.data.frame(cts_bulk.t)
dim(cts_bulk.t)
cts_bulk.t[1:10, 1:10]
cts_bulk.t[is.na(cts_bulk.t)] <-0 
cts_bulk.t[1:7, 1:10]
rownames(cts_bulk.t)[rownames(cts_bulk.t) == '1'] <- "pca_1"
rownames(cts_bulk.t)[rownames(cts_bulk.t) == '0'] <- "pca_0"
for (i in 2:6) {
  current_row_name <- as.character(i)
  new_row_name <- paste0("pca_", current_row_name)
  rownames(cts_bulk.t)[rownames(cts_bulk.t) == current_row_name] <- new_row_name
}
  
View(cts_bulk.t)

# Splitting the dataframe to separate the pca
factors <-rownames(cts_bulk.t)
cts_bulk.t_split <- split.data.frame(cts_bulk.t,
                                     f = factors)
View(cts_bulk.t_split)

cts_bulk_pca_sep <- lapply (cts_bulk.t_split, function(x){
  t(x)
})
cts_bulk_pca_sep$pca_1
View(cts_bulk_pca_sep)

#Followed by DESEQ2 analysis if we had ctrol and stim on the PCAs

#UMAP embedding based LDA

str(hcc364.seurat_hdf5)
hcc364.seurat_hdf5[['umap']]



