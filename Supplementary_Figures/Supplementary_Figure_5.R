# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)

options(future.globals.maxSize = 3000 * 1024^2)

# Load in processed xenium objects, epithelial dataset, and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

epi_fetal<-readRDS('data/epi_fetal_lung.rds')
epi_palette<-readRDS('palettes/epi_palette.rds')

# Select coordinates to zoom

## GW 15
## Coordinates 6 - Smaller airway - PNEC
GW15_coords6 <- Crop(xenium.gw15[["fov"]], x = c(12225, 12525), y = c(9125, 9425), coords = "plot")
xenium.gw15[["zoom"]] <- GW15_coords6
DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"

## GW 18
##Coordinates 9 - PNEC
GW18_coords9 <- Crop(xenium.gw18[["fov"]], x = c(3450, 3750), y = c(3150,3450), coords = "plot")
xenium.gw18[["zoom"]] <- GW15_coords9
DefaultBoundary(xenium.gw18[["zoom"]]) <- "segmentation"

features<-c('ASCL1', 'CALCA')

## Supplementary Figure 5A

# Switch when required
# xenium.obj<-xenium.gw15
# xenium.obj<-xenium.gw18

# Plot spatial dimplot

ImageDimPlot(xenium.obj, fov = "zoom", group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Supplementary_Fig_5A_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature;
  p <- ImageFeaturePlot(xenium.obj, fov = "zoom", features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Supplementary_Fig_5A_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}

## Supplementary Figure 5B

features<- c('GRP', 'GHRL')

FeaturePlot(epi_fetal, features = features)

ggsave('Supplementary_Fig_5B.pdf')

## Supplementary Figure 5C

# Define coordinates to zoom

## GW 15
## Coord 8 - PTEN+STAT3+ Cells
GW15_coords8_PTEN <- Crop(xenium.gw15[["fov"]], x = c(11375, 11675), y = c(7750, 8050), coords = "plot")
# xenium.gw15[["zoom"]] <- GW15_coords8_PTEN
# features <- c('PTEN', 'STAT3')

## Coord9 - NRGN+ cells 
GW15_coords9_NRGN <- Crop(xenium.gw15[["fov"]], x = c(10950, 11250), y = c(3000, 3300), coords = "plot")
# xenium.gw15[["zoom"]] <- GW15_coords9_NRGN
# features <- c('NRGN')

DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"
## GW 18
## Coord 15 - PTEN+STAT3+ cells
GW18_coords15_PTEN <- Crop(xenium.gw18[["fov"]], x = c(12850,13150), y = c(4470, 4770), coords = "plot")
# xenium.gw18[["zoom"]] <- GW18_coords8_PTEN
# features <- c('PTEN', 'STAT3')


## Coord 16 - NRGN+ cells
GW18_coords16_NRGN <- Crop(xenium.gw18[["fov"]], x = c(11925, 12225), y = c(6000, 6300), coords = "plot")
# xenium.gw18[["zoom"]] <- GW18_coords16_NRGN
# features <- c('NRGN')

## Switch gestational week and coordinates as needed

# xenium.obj<-xenium.gw15
# xenium.obj<-xenium.gw18

# Plot spatial dimplot

ImageDimPlot(xenium.obj, fov = "zoom", group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Supplementary_Fig_5C_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature;
  p <- ImageFeaturePlot(xenium.obj, fov = "zoom", features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Supplementary_Fig_5C_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}



