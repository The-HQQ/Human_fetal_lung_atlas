# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
# Set path to raw data

path_gw15<-'/all_data_Xenium/output_gw15/'
path_gw18<-'/all_data_Xenium/output_gw18/'

xenium_gw15<-LoadXenium(path_gw15, fov = "fov")
xenium_gw18<-LoadXenium(path_gw18, fov = "fov")

#Remove cells with 0 counts

xenium_gw15<-subset(xenium_gw15, subset = nCount_Xenium > 0)
xenium_gw18<-subset(xenium_gw18, subset = nCount_Xenium > 0)

xenium_objects<-c(xenium_gw15, xenium_gw18)

xenium_objects<-lapply(X=xenium_objects, FUN = function(x) {
  x<-SCTransform(x, assay = "Xenium")
  x<- RunPCA(x, npcs = 30, features = rownames(x))
  x<-RunUMAP(x, dims = 1:30)
  x<-FindNeighbors(x, reduction = 'pca', dims = 1:30)
  x<-FindClusters(x, resolution = c(0.3)
})

saveRDS(xenium_objects[[1]], '/data/xenium_gw15_processed.rds')
saveRDS(xenium_objects[[2]], '/data/xenium_gw18_processed.rds')
