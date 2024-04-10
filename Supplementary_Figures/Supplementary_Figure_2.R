# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

options(ggrepel.max.overlaps = Inf)

# Load relevant publically available datasets, Quach & Farrell et al. dataset, and palettes

Sikemma_adult_lung<-readRDS('/data/Sikemma_adult_lung.rds')
Negretti_mouse_lung<-readRDS('/data/Negretti_adult_lung.rds')
Quach_fetal_lung<-readRDS('/data/full_fetal_lung_dataset.rds')
all_fetal_palette<-readRDS('all_fetal_palette.rds')

# Set public datasets to to query and Quach & Farrell et al. dataset to reference

reference<-Quach_fetal_lung

query <- Sikemma_adult_lung

# Find anchors between datasets and map query cells to the references

anchors<- FindTransferAnchors (reference = reference, query = query, dims = 1:30, reference.reduction = "pca")

map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "all_cell_type"), reference.reduction = "pca", reduction.model = "umap")

threshold <- 0.9
map_query$celltype <- ifelse(map_query$prediction.score.max > threshold, map_query$predicted.id, "X")

## Supplementary Note Fig. 1A

DimPlot(map_query, reduction = "ref.umap", group.by = "celltype", label = TRUE, cells = WhichCells(map_query, idents = "X", invert = T), cols = all_fetal_palette, repel = T) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Fig_1A.pdf", width = 16.6, height = 6.6)



# Set Sikemma et al. adult lung to reference and Quach & Farrell et al. dataset to query

query<-Quach_fetal_lung

reference <- Sikemma_adult_lung

# Find anchors between datasets and map query cells to the references

anchors<- FindTransferAnchors (reference = reference, query = query, dims = 1:30, reference.reduction = "pca")

map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "ann_finest_level"), reference.reduction = "pca", reduction.model = "umap")

threshold <- 0.9
map_query$celltype <- ifelse(map_query$prediction.score.max > threshold, map_query$predicted.id, "X")

# Create Unimodal UMAP projection figures

## Supplementary Fig. 1B

DimPlot(map_query, reduction = "ref.umap", group.by = "celltype", label = TRUE, cells = WhichCells(map_query, idents = "X", invert = T), cols = all_fetal_palette, repel = T) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Fig_1B.pdf", width = 11, height = 8.5)

## Supplementary Note Fig. 1C
query <- Negretti_mouse_lung

# Find anchors between datasets and map query cells to the references

anchors<- FindTransferAnchors (reference = reference, query = query, dims = 1:30, reference.reduction = "pca")

map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "all_cell_type"), reference.reduction = "pca", reduction.model = "umap")

# Adjust for Negretti dataset
DimPlot(map_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Fig_1C.pdf", width = 16.6, height = 6.6)

## Frequencies to create the tables

# Rename the output file for each reference dataset
write.csv(table(map_query$predicted.celltype), "Supplementary_Fig_2_Frequency.csv")
