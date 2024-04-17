# Clear workspace
library(Seurat)

rm(list = ls())
## Insert path for each cellranger count output directory
data <- Read10X(data.dir = path)
# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "", min.cells = 3, min.features = 200)
data

# creates a new column
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

saveRDS(data, file = 'GW_10_2.rds')

# Repeat above steps with each sample
