# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions

setwd("~/Desktop/Single_cell_project")
set.seed(101)

library(Seurat)
library(tidyverse)

# Load data
ifnb_harmony <- readRDS('ifnb_harmony.rds')
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# Visualize data
clusters <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) + 
        ggtitle("Clusters by Seurat")
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim') + 
        ggtitle("Conditions")

ggsave("Condition_vs_Seurat_clusters.png", plot = condition|clusters, width = 12, height = 5, dpi = 300)


# FindAll markers --------------------------------------------------------------
FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# FindConserved markers --------------------------------------------------------
# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony) <- must be 'RNA', since we use DESeq2
# if not RNA, run this: DefaultAssay(ifnb_harmony) <- 'RNA'
DefaultAssay(ifnb_harmony)

# Identify the cluster 3: compare the cluster 3 with all the clusters
markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                     ident.1 = 3,
                     grouping.var = 'stim')
head(markers_cluster3)
# FCGR3A is expressed in 97.7% of cells in cluster 3 of the control group, and in 20.6% of cells in the other clusters.

# Visualize the top features: FCGR3A
# excludes cells with very low expression
FCGR3A_inCluster3 <- FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')
ggsave("FCGR3A_in_Cluster3_FeaturePlot.png", plot = FCGR3A_inCluster3, width = 6, height = 5, dpi = 300)

# Rename cluster 3 with CD16 Mono (known info)
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')
DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# Cells already have annotations provided in the metadata
View(ifnb_harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual annotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# Setting Idents as Seurat annotations provided (also a sanity check!)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)


# FindMarkers between conditions -----------------------------------------------
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)
View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

plot_cell_cond_label <- DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE) + 
        ggtitle("Seurat Clusters Labeled by Cell Type and Condition")
ggsave("Seurat_Clusters_CellType_Condition.png", plot = plot_cell_cond_label, width = 10, height = 5, dpi = 300)

# Find markers
# To compare the cells of the same cell type between conditions 
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')
head(b.interferon.response)
# IFIT1 is the top DE between control and STIM 


# Plotting conserved features vs DE features between conditions
head(markers_cluster3)
Conserved_vs_DE_Features <- FeaturePlot(ifnb_harmony, features = c('FCGR3A', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')
# 'FCGR3A' the most DE between cluster 3 and the other clusters
# 'IFIT1' the most DE between control and STIM in cluster 3
ggsave("Conserved_vs_DE_Features_in_Cluster3.png", plot = Conserved_vs_DE_Features, width = 6, height = 5, dpi = 300)


