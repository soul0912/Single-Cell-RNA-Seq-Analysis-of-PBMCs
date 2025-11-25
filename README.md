# Single-Cell-RNA-Seq-Analysis-of-PBMCs

This Project demonstrates the Pipeline for Single-Cell Rna Sequencing Anlysis(scRNA-Seq) using **Seurat Package** in R programming:
>>> Processing Raw sequening data of ~2700 Peripheral Blood Mononuclear Cells (PBMCs) provided by 10X Genomics.

## Analyzed By:
**Language** : R
**Packages** : Seurat, dplyr, ggplot2, patchwork.
**Data Source** : 10x Genomics (PBMC 3K Dataset)

## Methodology:

**1. Refining:**
>>> Refining low-quality cells by setting threshold for mitochondrial gene content (<5%) & total feature  count (>200 & <2500)
>>> It reduced the dataset from 2700 to 2638 viable cells.
>>> Log Normalization
>>> Identified the top 2,000 highly variable features to drive downstream clustering.

**2. PCA Analysis:**
>>> Trancriptomic Variance
>>> Generating **UMAP (Uniform Manifold Approximation and Projection)** for non-linear visualization of cell clusters.

**3. Clustering:**
>>>Constructed a K-Nearest Neighbor (KNN) graph.
>>>Applied the Louvain algorithm (Resolution = 0.5) to partition cells into 9 distinct clusters.

**4. Annotation:**
>>>Utilized FindAllMarkers() to identify cluster-specific gene signatures.

## Key Findings:

**T-Cells:** Naive CD4+, Memory CD4+, and CD8+ T-Cells (Markers: CCR7, LEF1, AQP3, CD4OLG, GZMK, GZMH).
**Monocytes:** Markers(CKB, CDKN1C)
**B-Cells:** Identified by VPREB3, LINCO.
**NK Cells:** Identified by SH2D1B, AKR1C3.
**Dendritic Cells & Platelets:** FCER1A, SERPIN, LY6G6F.
   
