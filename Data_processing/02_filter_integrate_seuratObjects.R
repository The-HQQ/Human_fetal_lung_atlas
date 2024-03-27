# Clear workspace
rm(list = ls())
# options(future.globals.maxSize = 18000 * 1024^2)

# Import dependencies
library(Seurat, verbose = FALSE)
library(tidyverse, verbose = FALSE)
library(fgsea, verbose = FALSE)
library(progeny, verbose = FALSE)
library(cowplot, verbose = FALSE)
library(ggplot2, verbose = FALSE)

# Sample name
sample_names <- c("GW_10_2", 
                  "GW_11_1", 
                  "GW_12_3", 
                  "GW_13_4", 
                  "GW_13_6",
                  "GW_14_3",
                  "GW_15_5", 
                  "GW_15_5_1", 
                  "GW_16_1", 
                  "GW_16_2", 
                  "GW_18", 
                  "GW_18_1", 
                  "GW_18_2", 
                  "GW_18D_1", "GW_18P_1",
                  "GW_19_0", 
                  "GW_19_2") 

# Predicted gender
gender <- c("Male",
            "Male",
            "Male",
            "Male",
            "Male",
            "Male",
            "Female",
            "Female",
            "Male",
            "Female",
            "Male",
            "Male",
            "Female",
            "Male", "Male",
            "Male",
            "Female")

# Batch number
batch_number <- c(1,
                  2,
                  3,
                  2,
                  4,
                  5,
                  3,
                  6,
                  2,
                  3,
                  9,
                  8,
                  9,
                  5, 5,
                  7,
                  9)

# Groups according to gestation week
sample_group <- c("early", 
                  "early", 
                  "early", 
                  "early", 
                  "early",
                  "mid",
                  "mid", 
                  "mid", 
                  "mid", 
                  "mid", 
                  "late", 
                  "late", 
                  "late", 
                  "late", "late",
                  "late", 
                  "late") 

# Gestational week
sample_week <- c("week_10", 
                  "week_11", 
                  "week_12", 
                  "week_13", 
                  "week_13",
                  "week_14",
                  "week_15", 
                  "week_15", 
                  "week_16", 
                  "week_16", 
                  "week_18", 
                  "week_18", 
                  "week_18", 
                  "week_18", "week_18",
                  "week_19", 
                  "week_19") 

# Directory path
# set_dir <- "/hpf/largeprojects/ccmbio/pkallurkar/scRNA-seq/"

# Set number of features to use for downstream analysis
number_of_features <- 3000
number_of_dimensions <- 30

# List to store the samples
seurat_objects = list()

#------------------ Load and prepare the samples for integration ---------------
# Create a function to calculate mitochondrial percentage
calculate_per_mit = function(seurat_object) {
  seurat_object[["per_mit"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
  return(seurat_object)
}

# Create a function to filter the samples
filter_seurat_object = function(seurat_object) {
  seurat_object <- subset(seurat_object, subset = per_mit < 15 & nFeature_RNA > 200)
  return(seurat_object)
}

# Create a function to load seurat objects, filter them and add information about the samples
load_sample = function(filename, gender, batch_number, sample_group, sample_week) {
  fetal_lung_sample <- readRDS(file = paste0(set_dir, "data/unprocessed data/", filename, ".rds"))
  fetal_lung_sample <- calculate_per_mit(fetal_lung_sample)
  fetal_lung_sample <- filter_seurat_object(fetal_lung_sample)
  fetal_lung_sample <- AddMetaData(object = fetal_lung_sample, metadata = filename, col.name="sample_name")
  fetal_lung_sample <- AddMetaData(object = fetal_lung_sample, metadata = gender, col.name="gender")
  fetal_lung_sample <- AddMetaData(object = fetal_lung_sample, metadata = batch_number, col.name="batch_number")
  fetal_lung_sample <- AddMetaData(object = fetal_lung_sample, metadata = sample_group, col.name="sample_group")
  fetal_lung_sample <- AddMetaData(object = fetal_lung_sample, metadata = sample_week, col.name="sample_week")
  return(fetal_lung_sample)
}

print("Loading samples...")
for (sample_index in 1:length(sample_names)) {
  seurat_objects <- append(seurat_objects, 
                           load_sample(sample_names[sample_index], 
                                       gender[sample_index],
                                       batch_number[sample_index], 
                                       sample_group[sample_index],
                                       sample_week[sample_index]))
}
names(seurat_objects) <- sample_names

# Verify if the filtering worked by plotting violin plot for a sample
sample_index = 1
pdf(paste(set_dir, "results/", number_of_features,"/", 
          "dim_", number_of_dimensions, "/vln_plot_", 
          sample_names[sample_index], ".pdf", sep= ""))
VlnPlot(object = seurat_objects[[sample_index]],
        features = c("nFeature_RNA", "nCount_RNA", "per_mit"),
        pt.size = 0.2,
        ncol = 3)
dev.off()

# Normalize and find top n variable features
print("Normalizing...")
seurat_objects <- lapply(X = seurat_objects,  FUN = function(x) {
  x <-  NormalizeData(x, normalization.method = "LogNormalize", verbose  = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = number_of_features, verbose = FALSE)
})

print("Selecting integration features...")
# Select features that are repeatedly variable across datasets for integration
fetal_lung_features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = number_of_features, verbose = FALSE)

print("Scaling data and then run PCA and UMAP...")
# Scale the data and run PCA on each object
seurat_objects <- lapply(X = seurat_objects, FUN = function(x) {
    x <- ScaleData(x, features = fetal_lung_features, verbose = FALSE)
    x <- RunPCA(x, features = fetal_lung_features, verbose = FALSE)
    x <- RunUMAP(x, reduction = "pca", dims = 1:number_of_dimensions, verbose = FALSE)
})

print("Plotting PCA...")
# Plot PCA for each sample
pca_plot_list = c()
for (object_index in 1:length(seurat_objects)) {
  pca_plot_list[[object_index]] <- DimPlot(seurat_objects[[object_index]], 
                                          reduction = "pca", label = FALSE)
  pca_plot_list[[object_index]] <- pca_plot_list[[object_index]] +
    ggtitle(sample_names[object_index]) +
    theme(
      plot.title = element_text(size = 8),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 8))
}
cowplot::plot_grid(plotlist = pca_plot_list)
ggsave(filename = paste0(set_dir, "results/", number_of_features,"/", "dim_", number_of_dimensions, "/pca_plot.png"), 
       width = 50, height = 40, units = "cm")

print("Plotting UMAP...")
# Plot UMAP for each sample
umap_plot_list = c()
for (object_index in 1:length(seurat_objects)) {
  umap_plot_list[[object_index]] <- DimPlot(seurat_objects[[object_index]], 
                                          reduction = "umap", label = FALSE)
  umap_plot_list[[object_index]] <- umap_plot_list[[object_index]] +
    ggtitle(sample_names[object_index]) +
    theme(
      plot.title = element_text(size = 8),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 8))
}
cowplot::plot_grid(plotlist = umap_plot_list)
ggsave(filename = paste0(set_dir, "results/", number_of_features,"/", "dim_", number_of_dimensions, "/umap_plot.png"), 
       width = 50, height = 40, units = "cm")

print("Saving the data...")
# Export the objects with PCA results
for (seurat_index in 1:length(seurat_objects)) {
  saveRDS(object = seurat_objects[[seurat_index]],
          file = paste0(set_dir, "data/", number_of_features,"/", "dim_", number_of_dimensions, "/", 
                        sample_names[seurat_index], ".rds"),
          compress = FALSE)
}

# ------------------ Integrate the sample ------------------
print("Finding integration anchors...")
# Find integration anchors
anchor_set <- FindIntegrationAnchors(object.list = seurat_objects,
                                     anchor.features = fetal_lung_features,
                                     scale = FALSE,
                                     normalization.method = "LogNormalize",
                                     reduction = "rpca",
                                     dims = 1:30)

print("Integrating the data...")
seurat_integrated <- IntegrateData(anchor_set,
                                   new.assay.name = "integrated", 
                                   normalization.method = "LogNormalize",   
                                   dims = 1:30, 
                                   verbose = FALSE)
rm(seurat_objects)
rm(anchor_set)

print("Scaling the integrated data...")
seurat_integrated <- ScaleData(seurat_integrated, assay = "integrated", verbose = FALSE)

print("Running PCA on the integrated data...")
seurat_integrated <- RunPCA(seurat_integrated, assay = "integrated", verbose = FALSE)

print("Running UMAP on the integrated data...")
seurat_integrated <- RunUMAP(seurat_integrated, assay = "integrated", reduction = "pca", dims = 1:30)

# ------------------ Generate plots --------------------------------------------
print("Generate plots for integrated data...")
# Group by sample
DimPlot(seurat_integrated, reduction = "pca", group.by = "sample_name") + 
  ggtitle(paste("PCA by sample \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/pca_integrated_sample_name.png"), width = 50, height = 40, units = "cm")

# Group by gender
DimPlot(seurat_integrated, reduction = "pca", group.by = "gender") +   
  ggtitle(paste("PCA by gender \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/pca_integrated_gender.png"), width = 50, height = 40, units = "cm")

# Group by batch
DimPlot(seurat_integrated, reduction = "pca", group.by = "batch_number") +
  ggtitle(paste("PCA by batch \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/pca_integrated_batch.png"), width = 50, height = 40, units = "cm")

# Group by sample group
DimPlot(seurat_integrated, reduction = "pca", group.by = "sample_group") + 
  ggtitle(paste("PCA by group \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/pca_integrated_sample_group.png"), width = 50, height = 40, units = "cm")

# Group by sample week
DimPlot(seurat_integrated, reduction = "pca", group.by = "sample_week") + 
  ggtitle(paste("PCA by gestation week \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/pca_integrated_sample_week.png"), width = 50, height = 40, units = "cm")

# Export UMAP plots
# Group by sample
DimPlot(seurat_integrated, reduction = "umap", group.by = "sample_name") + 
  ggtitle(paste("UMAP by sample \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/umao_integrated_sample_name.png"), width = 50, height = 40, units = "cm")

# Group by gender
DimPlot(seurat_integrated, reduction = "umap", group.by = "gender") +   
  ggtitle(paste("UMAP by gender \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/umao_integrated_gender.png"), width = 50, height = 40, units = "cm")

# Group by batch
DimPlot(seurat_integrated, reduction = "umap", group.by = "batch_number") +
  ggtitle(paste("UMAP by batch \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/umao_integrated_batch.png"), width = 50, height = 40, units = "cm")

# Group by sample group
DimPlot(seurat_integrated, reduction = "umap", group.by = "sample_group") + 
  ggtitle(paste("UMAP by group \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/umao_integrated_sample_group.png"), width = 50, height = 40, units = "cm")

# Group by gestation week
DimPlot(seurat_integrated, reduction = "umap", group.by = "sample_week") + 
  ggtitle(paste("UMAP by gestation week \n features = ", number_of_features, sep = "")) + 
  theme(plot.title = element_text(size = 10))
ggsave(filename = paste0(set_dir, "results/", number_of_features, "/", "dim_", number_of_dimensions, "/umao_integrated_sample_week.png"), width = 50, height = 40, units = "cm")


# ------------------ Export data------------------------------------------------
print("Export the integrated data...")
# Export integrated data with umap and pca reduction
saveRDS(object = seurat_integrated,
        file = paste0(set_dir, "data/", number_of_features, "/", "dim_", number_of_dimensions, "/", "integrated_data.rds"),
        compress = FALSE)

