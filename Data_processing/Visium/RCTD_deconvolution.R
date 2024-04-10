# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(spacexr)

# Load visium files and Quach & Farrell et al. dataset

vis.gw18<-readRDS('/data/vis.gw18.rds')
vis.gw15<-readRDS('/data/vis.gw15.rds')
all_fetal <- readRDS('full_fetal_lung_dataset.rds')

# Extract reference information

counts <- all_final[["RNA"]]@counts
cluster <- as.factor(all_final$all_cell_type)
names(cluster) <- colnames(all_final)
nUMI <- all_final$nCount_RNA
names(nUMI) <- colnames(all_final)
reference <- Reference(counts, cluster, nUMI)

# Extract query information from visium data; repeat for each GW

counts <- vis.gw15[["Spatial"]]@counts
coords <- GetTissueCoordinates(vis.gw15)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

# Run RCTD deconvolution function on "full" mode

RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
vis.gw15 <- AddMetaData(vis.gw15, metadata = RCTD@results$results_df)

saveRDS(vis.gw15, '/data/vis.gw15.rds')

# Repeat similar code for gw 18

