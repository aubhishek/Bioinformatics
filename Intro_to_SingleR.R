# script to annotate cell types from 20k Human PBMCs from a healthy female donor
setwd("/Users/az/Desktop/YouTubeTutorials-main/files/Find_all_markers")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Input Data 10X CellRanger .HDF5 format --------------
hdf5_obj <- Read10X_h5(filename = '20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering -----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)
ElbowPlot (pbmc.seurat.filtered)


# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# Things that do not work as a gene list: metadata does not have genes in there!!!
#So we cannot slice it
#features_sc <- subset(pbmc.seurat.filtered@meta.data, subset = seurat_clusters == 3)
#View(features_sc)

#features_sc_2 <- pbmc.seurat.filtered[[pbmc.seurat.filtered$seurat_clusters == 3]]


#gene_list <- unique(pbmc.seurat.filtered$seurat_clusters)
#View(gene_list)
#subset_dotplot <- DotPlot(pbmc.seurat.filtered, features = gene_list)

View (pbmc.seurat.filtered@assays$RNA@counts)


#How to print gene_names
head(pbmc.seurat.filtered[['RNA']])
gene_names <- rownames(pbmc.seurat.filtered@assays$RNA)
print(gene_names)

#How to add gene names on metadata (metadata currently have only barcode)
#1: isolate barcode names from assay (barcode is the key between the two matrix)
assay_barcodes <- colnames(pbmc.seurat.filtered@assays$RNA)
head(assay_barcodes)
length(assay_barcodes)

#isolate gene names
gene_names <- rownames(pbmc.seurat.filtered@assays$RNA)
head(gene_names)
length(gene_names)


#Now merge
updated_metadata <- pbmc.seurat.filtered@metadata[pbmc.seurat.filtered@metadata$barcode %in% assay_barcodes, ]
updated_metadata$gene_names <- gene_names




pbmc.seurat.filtered@meta.data$RNA <- pbmc.seurat.filtered@assays$RNA
View(pbmc.seurat.filtered@meta.data)

DotPlot(pbmc.seurat.filtered, features = c("YAP1", "GAPDH", "LAMP1", "AKT1", "FGFR1"))
DotPlot(pbmc.seurat.filtered, features = features_sc)
DotPlot(pbmc.seurat.filtered, features = pbmc.seurat.filtered$seurat_clusters)

levels(pbm)



# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)




# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))




