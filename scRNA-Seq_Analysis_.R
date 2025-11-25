
BiocManager::install("Seurat")

library(Seurat)
library(dplyr)

## Optimizing Path 
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19")
## Creating the Object 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
## Filtering the data according to the Mitchondrial Content 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

print(pbmc)
# Normalization.
pbmc <- NormalizeData(pbmc)
# Variable Feature: No need of 13,714 genes, Let just proceed with top 2,000 genes.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Scaling data 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# PCA Analysis(Linear Dimensionality Reduction)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Find Neighbours and clusters for Similarity
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Running UMAP(Non_Linear Dimensionality Reduction)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
# Find markers for every cluster 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Showing only top 2 genes for each cluster
pbmc.markers %>%
group_by(cluster) %>%
slice_max(n=2, order_by = avg_log2FC)
new.clusters.ids <- c("Naive CD4 T", "CD14 + Mono", "Memory CD4 T", "Antibodies", "CD8 T", "Monocytes", "NK cells", "Dendritic Cells", "Platelets")
names(new.clusters.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.clusters.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# Renaming one marker, it's not similar
pbmc <- RenameIdents(pbmc, "Antibodies" = "B Cells")

library(ggplot2)

final_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5) +
ggtitle("Single-Cell Analysis of 2,700 PBMCs")

print(final_plot)

