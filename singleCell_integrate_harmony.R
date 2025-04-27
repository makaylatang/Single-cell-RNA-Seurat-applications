# script to integrate across conditions using Harmony
# setwd("~/Desktop/single_cell_project/")

# set seed for reproducibility
set.seed(1234)

#install.packages("remotes")
#remotes::install_github("satijalab/seurat-data")

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# get data -------------------------
AvailableData()
# install dataset
options(timeout = 600)  # timeout set to 600 seconds (10 minutes)
SeuratData::InstallData("ifnb")

# load dataset
ifnb <- LoadData("ifnb")
str(ifnb)


# QC and filtering
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)
# explore QC

# filter
ifnb # 14053 features across 13999 samples within 1 assay 
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)
# 14053 features across 13988 samples within 1 assay 

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')


# run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

before|after
