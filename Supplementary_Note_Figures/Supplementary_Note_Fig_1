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

# Because He et al. dataset labels are too long, we encode each annotation with a number

Idents(He_fetal_lung) <- "new_celltype"
annotations<- levels(He_fetal_lung)
names(annotations) <- c(1:144)

He_fetal_lung<-RenameIdents(He_fetal_lung, annotations)

He_fetal_lung$num_ident<- Idents(He_fetal_lung)

# Set public datasets to reference and Quach & Farrell et al. dataset to query

query<-Quach_fetal_lung

# reference <- He_fetal_lung
# reference <- Sountoulidis_fetal_lung
# reference <- Cao_fetal_lung

# Find anchors between datasets and map query cells to the references

anchors<- FindTransferAnchors (reference = reference, query = query, dims = 1:30, reference.reduction = "pca")

# For He dataset
map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "num_ident"), reference.reduction = "pca", reduction.model = "umap")
# For Sountoulidis dataset
map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "indiv_clusters"), reference.reduction = "pca", reduction.model = "umap")
# For Cao dataset
map_query <- MapQuery(anchorset = anchors, reference = reference, query = query,
    refdata = list(celltype = "indiv_clusters"), reference.reduction = "pca", reduction.model = "umap")

# Create Unimodal UMAP projection figures

## Supplementary Note Fig. 1A

DimPlot(map_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Note_Fig_1B.pdf", width = 16.6, height = 6.6)

## Supplementary Note Fig. 1B

DimPlot(map_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Note_Fig_1B.pdf", width = 16.6, height = 6.6)

## Supplementary Note Fig. 1C

DimPlot(map_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
ggsave("Supplementary_Note_Fig_1C.pdf", width = 16.6, height = 6.6)

## Frequencies to create the tables

# Rename the output file for each reference dataset
write.csv(table(map_query$predicted.celltype), "Supplementary_Note_Frequency.csv")
