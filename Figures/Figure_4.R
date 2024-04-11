# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load epithelial dataset & palette
epi_fetal<-readRDS('data/epi_fetal_lung.rds')
epi_palette<-readRDS('palettes/epi_palette.rds')

## Fig 4A
p1<-DimPlot(epi_fetal, cols = epi_palette, label = F)
p1$data$num_ident<-epi_fetal@meta.data$num_ident
LabelClusters(p1, id ='num_ident')

ggsave("Fig_4A.pdf", width = 15, height = 9)

## Fig 4B

# Load DEG table for epithelial cluster

DEG_markers<-read.csv('c2_markers.csv')
top3<- DEG_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
features = top3$gene

# DEG features that were chosen
chosen_features = cc('CFTR','SFTPB','NKX2-1', 'SOX2', 'SOX9', 'FOXJ1', 'NRGN', 'PTEN', 'STAT3', 'BRCA1', 'FOXM1','MKI67', 'NUSAP1', 'CDK1', 'TEAD4', 'CPM', 'SFTPA', 'FGFR2', 'MUC1', 'SMAD6', 'SCGB1A1', 'RSPH1', 'DNAH5', 'TUBA1A', 'SFTPA1', 'ABCA3', 'CA2', 'TESC', 'ETV5', 'MFSD2A', 'IGFPB7', 'CHGB', 'CALCA', 'ASCL1')
epi_features <- c(features, chosen_features)

# Change the order of the features to match up with the figure; or can load in the features already ordered
ordered_epi_features <- readRDS('data/epi_DEG_chosen.rds')

# Ensure selected features are scaled

epi_fetal<-ScaleData(epi_fetal, features = ordered_epi_features)

# Generate heatmap
DoHeatmap(epi_fetal, assay = 'RNA', features = ordered_epi_features, size = 4, angle = 0, group.colors = epi_palette, label = F) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30))
ggsave("Fig_4B.pdf", width = 40, height = 18)

## Fig 4C

# Define genes to plot

genes <- c("SOX2", "SOX9", "CFTR")

DotPlot(object = epi_fetal, features=genes) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16)) + scale_size(range = c(2, 5)) +
    guides(size=guide_legend(title = 'Percent Expressed',override.aes=list(shape=21, colour="black", fill="white"))) + scale_color_gradient2(low= "blue", high = "red")

ggsave("Fig_4C.pdf", width = 6, height = 7)

## Fig 4D

# Load in processed xenium objects and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

options(future.globals.maxSize = 3000 * 1024^2)

# Select coordinates to zoom

GW15_coords1 <- Crop(xenium.gw15[["fov"]], x = c(11650, 11950), y = c(5550, 5850), coords = "plot")

xenium.gw15[["zoom"]] <- GW15_coords1
DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"

GW18_coords1 <- Crop(xenium.gw18[["fov"]], x = c(6470, 6770), y = c(1850, 1550), coords = "plot")
GW18_coords2 <- Crop(xenium.gw18[["fov"]], x = c(3700, 4000), y = c(2825, 3125), coords = "plot")
GW18_coords3 <- Crop(xenium.gw18[["fov"]], x = c(3250, 3550), y = c(3750, 4050), coords = "plot")

# Change zoom each time when needed

# xenium.gw18[["zoom"]] <- GW18_coords1
# xenium.gw18[["zoom"]] <- GW18_coords2
# xenium.gw18[["zoom"]] <- GW18_coords3
DefaultBoundary(xenium.gw18[["zoom"]]) <- "segmentation"

# Define features to plot

#features <- c('NKX2-1', 'SOX9', 'CFTR')
#features <- c('SCGB3A2', 'SFTPB', 'CFTR')

# Define gestational week and zoom to plot

# xenium.obj<- xenium.gw15
# xenium.obj<- xenium.gw18
# colours <- xenium.gw15_palette
# colours <- xenium.gw18_palette

# Plot spatial dimplot

ImageDimPlot(xenium.obj, fov = "zoom", group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Fig_4D_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature;
  p <- ImageFeaturePlot(xenium.obj, fov = "zoom", features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Fig_4D_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}

