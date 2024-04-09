# Clear workspace
rm(list = ls())

# Load libraries
library(Seurat, verbose = FALSE)
library(tidyverse, verbose = FALSE)
library(leiden, verbose = FALSE)
library(igraph, verbose = FALSE)

# set_dir <- "/hpf/largeprojects/ccmbio/pkallurkar/scRNA-seq/"

data.big <- readRDS(paste0(set_dir, "data/3000/dim_30/clustered_data.rds"))
DefaultAssay(data.big) <- "integrated"
data.big

# Change the ident to select cluster with resolution 0.01
Idents(data.big) <- data.big$integrated_snn_res.0.01
levels(data.big)

# Select c1 sub-cluster; repeat for each sub-cluster
data.sub <- subset(data.big, subset = integrated_snn_res.0.01 == "1")
data.sub

data.sub <- FindNeighbors(data.sub, dims = 1:30, verbose = FALSE)
data.sub <- FindClusters(data.sub, resolution = c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001), verbose = FALSE, algorithm = 4, method = "igraph", save.SNN = TRUE)
data.sub <- RunUMAP(data.sub, reduction = "pca", dims = 1:30, verbose = FALSE)

clustree(data.big, prefix = "integrated_snn_res.")
ggsave(filename = paste0(set_dir, "results/3000/dim_30/cluster_tree_c1.png"), width = 40, height = 40, units = "cm")

# Save the object
saveRDS(object = data.sub, 
        file = paste0(set_dir, "data/3000/dim_30/c1_cluster.rds"), 
        compress = FALSE)

# Change the assay to RNA for finding the marker genes
DefaultAssay(data.sub) <- "RNA"

# Find the markers
print("Finding markers")
data.sub.markers <- FindAllMarkers(data.sub, only.pos = TRUE, assay = "RNA", verbose = TRUE)

# Order marker genes by fold change
print("Arranging markers by descending values of fold change")
data.sub.markers <- data.sub.markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(data.sub.markers, file = paste0(set_dir, "results/3000/dim_30/c1_markers.csv"), quote = FALSE)
