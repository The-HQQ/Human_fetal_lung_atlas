# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load relevant publically available datasets and Quach & Farrell et al. dataset

He_fetal_lung<-readRDS('/data/He_fetal_lung.rds')
Sountoulidis_fetal_lung<-readRDS('/data/Sountoulidis_fetal_lung.rds')
Cao_fetal_lung<-readRDS('/data/Cao_fetal_lung.rds')

Quach_fetal_lung<-readRDS('/data/full_fetal_lung_dataset.rds')

# Set common identity across all datasets
Idents(He_fetal_lung) <- "He et al. 2022 fetal lung"
He_fetal_lung$dataset <- Idents(He_fetal_lung)

Idents(Sountoulidis_fetal_lung) <- "Sountoulidis et al. 2023 fetal lung"
Sountoulidis_fetal_lung$dataset <- Idents(Sountoulidis_fetal_lung)

Idents(Cao_fetal_lung) <- "Cao et al. 2023 fetal lung"
Cao_fetal_lung$dataset <- Idents(Cao_fetal_lung)

Idents(Quach_fetal_lung) <- "Quach & Farrell et al. 2023 fetal lung"
Quach_fetal_lung$dataset <- Idents(Quach_fetal_lung)

# Prepare all datasets for integration

dataset_list<-c(He_fetal_lung, Sountoulidis_fetal_lung, Cao_fetal_lung, Quach_fetal_lung)

dataset_list <- lapply(X = dataset_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = dataset_list)
dataset_list <- lapply(X = dataset_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integrate all datasets using Seurat rpca
all_fetal.anchors <- FindIntegrationAnchors(object.list = dataset_list, anchor.features = features, reduction = "rpca")

all_fetal.combined <- IntegrateData(anchorset = all_fetal.anchors)
DefaultAssay(all_fetal.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
all_fetal.combined <- ScaleData(all_fetal.combined, verbose = FALSE)
all_fetal.combined <- RunPCA(all_fetal.combined, npcs = 30, verbose = FALSE)
all_fetal.combined <- RunUMAP(all_fetal.combined, reduction = "pca", dims = 1:30)
all_fetal.combined <- FindNeighbors(all_fetal.combined, reduction = "pca", dims = 1:30)
all_fetal.combined <- FindClusters(all_fetal.combined, resolution = 0.5)

saveRDS(all_fetal.combined, '/data/public_quach_integrated_fetal_datsets.rds')

# Run integration with all public datasets without Quach & Farrell et al. dataset

dataset_list<-c(He_fetal_lung, Sountoulidis_fetal_lung, Cao_fetal_lung)

# Prepare all datasets for integration

dataset_list<-c(He_fetal_lung, Sountoulidis_fetal_lung, Cao_fetal_lung, Quach_fetal_lung)

dataset_list <- lapply(X = dataset_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = dataset_list)
dataset_list <- lapply(X = dataset_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integrate all datasets using Seurat rpca
public_fetal.anchors <- FindIntegrationAnchors(object.list = dataset_list, anchor.features = features, reduction = "rpca")

public_fetal.combined <- IntegrateData(anchorset = public_fetal.anchors)
DefaultAssay(public_fetal.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
public_fetal.combined <- ScaleData(public_fetal.combined, verbose = FALSE)
public_fetal.combined <- RunPCA(public_fetal.combined, npcs = 30, verbose = FALSE)
public_fetal.combined <- RunUMAP(public_fetal.combined, reduction = "pca", dims = 1:30)
public_fetal.combined <- FindNeighbors(public_fetal.combined, reduction = "pca", dims = 1:30)
public_fetal.combined <- FindClusters(public_fetal.combined, resolution = 0.5)

saveRDS(public_fetal.combined, '/data/public_fetal_datset.rds')



