# Clear workspace
rm(list = ls())

# Load libraries
library(Seurat)
library(tidyverse, verbose = FALSE)
library(leiden, verbose = FALSE)
library(igraph, verbose = FALSE)
library(clustree, verbose = FALSE)

#set_dir <- "/hpf/largeprojects/ccmbio/pkallurkar/scRNA-seq/"

data.big <- readRDS(paste0(set_dir, "data/3000/dim_30/integrated_data.rds"))
DefaultAssay(data.big) <- "integrated"


# Find clusters 
data.big <- FindNeighbors(data.big, dims = 1:30, verbose = FALSE)
data.big <- FindClusters(data.big, resolution = c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001), verbose = FALSE, algorithm = 4, method = "igraph", save.SNN = TRUE)
data.big <- RunUMAP(data.big, reduction = "pca", dims = 1:30, verbose = FALSE)

# Export the clustered data
saveRDS(data.big, file = paste0(set_dir, "data/3000/dim_30/clustered_data.rds"), compress = FALSE)

clustree(data.big, prefix = "integrated_snn_res.")
ggsave(filename = paste0(set_dir, "results/3000/dim_30/cluster_tree_main.png"), width = 40, height = 40, units = "cm")

# Change the assay to RNA for finding the marker genes
DefaultAssay(data.big) <- "RNA"
data.big

# Change the ident to indicate the cluster labels from selected resolution
Idents(data.big) <- data.big$integrated_snn_res.0.01
levels(data.big)

# Find the markers
print("Finding markers")
all_markers <- FindAllMarkers(data.big, logfc.threshold = 0.5, only.pos = TRUE, assay = "RNA", min.pct = 0.25, min.diff.pct = 0.25, verbose = TRUE)


# Order marker genes by fold change
print("Arranging markers by descending values of fold change")
all_markers <- all_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(all_markers, file = paste0(set_dir, "results/3000/dim_30/all_markers_res_0.01.csv"), quote = FALSE)
