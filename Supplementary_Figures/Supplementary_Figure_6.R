# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(scales)
library(viridis)
library(Matrix)
library(rgl)
library(clusterExperiment)
library(gam)

## Supplementary Figure 6B

# Load epithelial dataset & palette
epi_fetal<-readRDS('data/epi_fetal_lung.rds')
epi_palette<-readRDS('palettes/epi_palette.rds')

# Subset based on cell trajectory inferred by trajectory analysis and late gestational ages (GW 18-19)

slingshot_idents<-c("Basal cells (19)", "Club cells (16)", "PNEC (12)", "SCGB3A2+SFTPB+CFTR+ cells (9)", "SCGB3A2+FOXJ1+ cells (13)", "Mature ciliated cells (3)", "Budtip progenitors (1)", "NKX2-1+SOX9+CFTR+ cells (11)", "SOX2lowCFTR+ cells (6)")
epi_slingshot<-subset(epi_fetal, idents = slingshot_idents)
Idents(epi_slingshot)<- 'sample_week'
epi_slingshot_late<-subset(epi_slingshot, idents= c("week_18", "week_19"))

Idents(epi_slingshot_late)<- 'cell_type'

# Process subsetted data
epi_slingshot_late <- FindVariableFeatures(epi_slingshot_late)
epi_slingshot_late <- ScaleData(epi_slingshot_late)
epi_slingshot_late <- RunPCA(epi_slingshot_late)

#Setup object for slingshot
levels(epi_slingshot_late)<-sort(levels(epi_slingshot_late))

dimred <- epi_slingshot_late@reductions$pca@cell.embeddings
clustering <- Idents(epi_slingshot_late)
counts <- as.matrix(epi_slingshot_late@assays$RNA@counts[epi_slingshot_late@assays$RNA@var.features, ])

# Run slingshot analysis

set.seed(1)
sce <- slingshot(Embeddings(epi_slingshot_late, 'pca'), clusterLabels = clustering, start.clus = "Budtip progenitors (1)", reducedDim = 'PCA')
slingshot_df<-data.frame(sce@assays@data@listData[["pseudotime"]])
sce_ds<-SlingshotDataSet(sce)

# Setup palette

epi_palette_slingshot <-epi_palette[names(epi_palette) %in% levels(clustering)]

# Plot 3D UMAP

# Import correct orientation/perspective

pp <-dget('plot3depiView.R')

plot3d(Embeddings(epi_slingshot_late, 'pca'), col = epi_palette_slingshot[clustering], pch = 16, cex = 0.5)
plot3d.SlingshotDataSet(sce_ds, lwd = 5, add = TRUE)
par3d(pp)
