library(dplyr)
library(Seurat)
library(patchwork)

rm(list = ls())
#path = "/WongProject/yard/run_cellranger_count/"
## Insert path for each cellranger count output directory
data <- Read10X(data.dir = path)
# Initialize the Seurat object with the raw (non-normalized data).
### GW_10_2 has 5937 single cells, 17271 genes/features
data <- CreateSeuratObject(counts = data, project = "", min.cells = 3, min.features = 200)
data

# creates a new column
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

saveRDS(data, file = 'GW_10_2.rds')
