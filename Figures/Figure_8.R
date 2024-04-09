# Clear workspace
rm(list = ls())

# Load packages
library(foreach)
library(data.table)
library(pheatmap)
library(tidyverse)
library(gplots)
library(Seurat)
library(tidyverse)

# Load datasets

hPSC_fetal_cells<-readRDS('data/hPSC_fetal_lung_cells.rds')
hPSC_fetal_organoids<-readRDS('data/hPSC_fetal_organoids.rds')
adult_lung<-readRDS('data/Adult_lung.rds')
epi_fetal<-readRDS('data/epi_fetal_lung.rds')

# Subset adult lung datatset to adult epithelium

Idents(adult_lung) <- "free_annotation"
adult_lung_epi<-subset(Adult_lung, idents = c("Basal", "Club", "Mucous", "Ciliated", "Alveolar Epithelial Type 1", "Alveolar Epithelial Type 2"))

# Pre-processing hPSC fetal lung and adult lung datasets

DefaultAssay(hPSC_fetal_cells) <- "RNA"
hPSC_fetal_cells <- NormalizeData(hPSC_fetal_cells, normalization.method = "LogNormalize")
hPSC_fetal_cells <- FindVariableFeatures(hPSC_fetal_cells, selection.method = "vst")
hPSC_fetal_cells <- ScaleData(hPSC_fetal_cells)
hPSC_fetal_cells <- RunPCA(hPSC_fetal_cells, features = VariableFeatures(object = hPSC_fetal_cells))
hPSC_fetal_cells <- FindNeighbors(hPSC_fetal_cells, dims = 1:50, reduction = "pca")

# Set up metadata to distinguish datasets

Idents(hPSC_fetal_cells) <- "hPSC fetal lung cells"
hPSC_fetal_cells$development_stage<-Idents(hPSC_fetal_cells)
DefaultAssay(hPSC_fetal_cells)<- "RNA"

Idents(hPSC_fetal_organoids) <- "hPSC fetal lung organoids"
hPSC_fetal_organoids$development_stage<-Idents(hPSC_fetal_organoids)
DefaultAssay(hPSC_fetal_organoids)<- "RNA"

Idents(epi_fetal)<- "sample_week"
epi_fetal$development_stage<-Idents(epi_fetal)
DefaultAssay(epi_fetal) <-"RNA"

Idents(adult_lung_epi) <- "Adult lung"
adult_lung_epi$development_stage <- Idents(adult_lung_epi)

# Integrate tissue-derived lung epithelium

epi_list = c(adult_lung_epi, epi_fetal)
features <- SelectIntegrationFeatures(epi_list, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = epi_list, anchor.features = features)
epi_fetal_adult <- IntegrateData(anchorset = anchors)

DefaultAssay(epi_fetal_adult) <- "integrated"

epi_fetal_adult <- ScaleData(epi_fetal_adult, verbose = FALSE)
epi_fetal_adult <- RunPCA(epi_fetal_adult, npcs = 30, verbose = FALSE)
epi_fetal_adult <- RunUMAP(epi_fetal_adult, reduction = "pca", dims = 1:30)
epi_fetal_adult <- FindNeighbors(epi_fetal_adult, reduction = "pca", dims = 1:30)
epi_fetal_adult <- FindClusters(epi_fetal_adult, resolution = 0.5)

# Find top 50 DEG between tissue-derived lung epithelium

Idents(epi_fetal_adult) <- "development_stage"
DefaultAssay(epi_fetal_adult) <- "RNA"
epi_fetal_adult_DEG<-FindAllMarkers(epi_fetal_adult, only.pos = TRUE, logfc.threshold = 0.5)

epi_fetal_adult_DEG %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) -> top50

features <- top50$gene

# Ensure top genes are present in all datasets

top_genes_fetal<-intersect(top50$gene, rownames(hPSC_fetal_cells))
top_genes_organoids<-intersect(top50$gene, rownames(hPSC_fetal_organoids))

# Calculate average log-normalized gene expression based on development stage

epi_fetal_adult_features_fetal <- AverageExpression(epi_fetal_adult, features = top_genes_fetal, group.by = "development_stage", assays = 'RNA')
epi_fetal_adult_features_fetal_df<-as.data.frame(epi_fetal_adult_features_fetal)

epi_fetal_adult_features_organoids <- AverageExpression(epi_fetal_adult, features = top_genes_organoids, group.by = "development_stage", assays = 'RNA')
epi_fetal_adult_features_organoids_df<-as.data.frame(epi_fetal_adult_features_organoids)

hPSC_fetal_features_fetal<-AverageExpression(hPSC_fetal_cells, features = top_genes_fetal, group.by = "development_stage", assays = 'RNA')
hPSC_fetal_features_fetal_df<-as.data.frame(hPSC_fetal_features_fetal)

hPSC_organoids_features_organoids<-AverageExpression(hPSC_fetal_organoids, features = top_genes_organoids, group.by = "development_stage", assays = 'RNA')
hPSC_organoids_features_organoids_df<-as.data.frame(hPSC_organoids_features_organoids)

# Calculate Spearman coefficient 

fetal_cells_cor.vec<-cor(epi_fetal_adult_features_fetal_df, hPSC_fetal_features_fetal_df, method = 'spearman')
organoids_cor.vec<-cor(epi_fetal_adult_features_organoids_df, hPSC_organoids_features_organoids_df, method = 'spearman')

# Plot heatmap

rownames(fetal_cells_cor.vec)<-rows
colnames(fetal_cells_cor.vec)<-"hPSC fetal lung cells"
fetal_cells_cor.vec<-t(fetal_cells_cor.vec)

rownames(organoids_cor.vec)<-rows
colnames(organoids_cor.vec)<-"hPSC organoids"
organoids_cor.vec<-t(organoids_cor.vec)

y=rbind(fetal_cells_cor.vec, organoids_cor.vec)

y = t(y)

order<- colnames(y)

pheatmap(y, color = brewer.pal(11, 'RdBu'), cluster_cols = FALSE, cluster_rows = FALSE,treeheight_row = 0, treeheight_col = 0, display_numbers = round(y,2), fontsize = 15, number_color = 'white', angle_col = 0)
