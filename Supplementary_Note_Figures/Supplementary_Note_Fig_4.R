# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load visium objects and RCTD results
vis.gw18<-readRDS('/data/vis.gw18.rds')
vis.gw15<-readRDS('/data/vis.gw15.rds')

# RCTD<-readRDS('/data/RCTD_full_gw15.rds')
# RCTD<-readRDS('/data/RCTD_full_gw18.rds')

# Normalize the cell type proportions; repeat for each GW

barcodes_gw15<-colnames(RCTD@spatialRNA@counts)
weights_gw15 <- RCTD@results$weights
norm_weights_gw15 <-normalize_weights(weights_gw15)
cell_types_gw15<-colnames(norm_weights_gw15)
coords_gw15<-GetTissueCoordinates(vis_gw15)

# Find barcodes for top two visium spots with high SCGB3A2+SFTPB+CFTR+ cell weight

cell_type_interest<- 'SCGB3A2.SFTPB.CFTR..cells'

barcodes_gw15<-c()

barcodes_gw15[[cell_type_interest]]<-row.names(head(norm_weights_gw15[order(norm_weights_gw15[[cell_type_interest]], decreasing = T),], 2))

# Find coordinates of the above barcodes

coords_interest_gw15[[cell_type_interest]][[barcodes_gw15[[cell_type_interest]][1]]]<-coords_gw15[barcodes_gw15[[cell_type_interest]][1],]
coords_interest_gw15[[cell_type_interest]][[barcodes_gw15[[cell_type_interest]][2]]]<-coords_gw15[barcodes_gw15[[cell_type_interest]][2],]

# Define features to plot

features<- c('WNT2', 'WNT3', 'WNT4', 'WNT6', 'WNT2B', 'WNT9A', 'FZD2', 'LRP5', 'FZD3', 'LRP6', 'FZD6', 'FZD7')

# Plot feature plots of zoomed in coords

vis_gw15_subset<-subset(vis.gw15, subset = imagerow > coords_gw15_interest[[cell_type_interest]][[1]][["imagerow"]] - 20 & imagerow < coords_gw15_interest[[cell_type_interest]][[1]][["imagerow"]] + 20 & imagecol > coords_gw15_interest[[cell_type_interest]][[1]][["imagecol"]] - 20 & imagecol < coords_gw15_interest[[cell_type_interest]][[1]][["imagecol"]] + 20)
norm_weights_gw15_subset<-norm_weights_gw15[colnames(vis_gw15_subset),]
p1<-SpatialFeaturePlot(vis_gw15_subset, crop = T, features = features, pt.size.factor = 20)
p1
ggsave(paste('Supplementary_Note_Fig_4B_GW_15.png', sep = ''), height = 15, width = 15, units = 'in', dpi = 300)
}

# Repeat above code for gw 18
