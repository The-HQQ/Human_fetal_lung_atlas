# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)

## Insert path for each cellranger count output directory
# path_15 = 'raw_data/visium/gw15/outs/
# path_18 = 'raw_data/visium/gw18/outs/

vis.gw15 <-Load10X_Spatial(path_15)
vis.gw18<-Load10X_Spatial(path_18)

# Quality control

vis.gw18<- subset(vis.gw18, subset = nFeature_Spatial > 300 & percent_mito < 10)
vis.gw15<- subset(vis.gw15, subset = nFeature_Spatial > 300 & percent_mito <10)

# Preprocessing

vis.gw15 <- SCTransform(vis.gw15, assay = "Spatial", verbose = FALSE)
vis.gw18 <- SCTransform(vis.gw18, assay = "Spatial", verbose = FALSE)

preprocessing_process<-function(data.sub){
  data.sub<-RunPCA(data.sub, assay = 'SCT', verbose = F)
  data.sub <- FindNeighbors(data.sub, dims = 1:30, verbose = FALSE, reduction = "pca")
  data.sub <- FindClusters(data.sub, verbose = F)
  data.sub <- RunUMAP(data.sub, reduction = "pca", dims = 1:30, verbose = FALSE,return.model=TRUE)
}

vis.gw15<-preprocessing_process(vis.gw15)
vis.gw18<-preprocessing_process(vis.gw18)

saveRDS(vis.gw18, '/data/vis.gw18.rds')
saveRDS(vis.gw15, '/data/vis.gw15.rds')
