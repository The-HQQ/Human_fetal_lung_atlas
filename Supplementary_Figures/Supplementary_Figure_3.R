# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

options(future.globals.maxSize = 3000 * 1024^2)

# Load in processed xenium objects and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

# Select coordinates to zoom and define features to plot for each

## GW 15
## Coord 3 - fibroblasts
GW15_coords3 <- Crop(xenium.gw15[["fov"]], x = c(9000, 9300), y = c(1675, 1975), coords = "plot")
# features <- c('TOP2A', 'CENPF', 'MKI67') 
# features <- c('FOXM1', 'HELLS') 

## Coord 4 - Airway vs vascular SMC
GW15_coords4 <- Crop(xenium.gw15[["fov"]], x = c(9180, 9480), y = c(2100, 2400), coords = "plot")
# features <- c('TAGLN', 'DES', 'MEF2C') 

##Coord 10 - Chondrocytes
GW15_coords10 <- Crop(xenium.gw15[["fov"]], x = c(6550, 6850), y = c(5400, 5700), coords = "plot")
# features <- c('THBS1', 'SOX9')

## GW 18
## Coordinates 10 - fibroblasts
GW18_coords10 <- Crop(xenium.gw18[["fov"]], x = c(13370, 13670), y = c(7000, 7300), coords = "plot")
# features <- c('TOP2A', 'CENPF', 'MKI67') 
# features <- c('FOXM1', 'HELLS') 

## Coord 17 - Airway vs vascular SMC
GW18_coords17 <- Crop(xenium.gw18[["fov"]], x = c(10675, 10975 ), y = c(3380, 3680), coords = "plot")
# features <- c('TAGLN', 'DES', 'MEF2C') 

## Coord 18 - Chondrocytes
GW18_coords18 <- Crop(xenium.gw18[["fov"]], x = c(8100, 8400), y = c(3925, 4225), coords = "plot")
# features <- c('THBS1', 'SOX9')

# Subset xenium based on coordinates; repeat for each field of view and gestational week

xenium.gw15[["zoom"]] <- GW15_coords3
DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"

# Generate spatial plot;repeat for each field of view and gestational week

ImageDimPlot(xenium.gw15, fov = "zoom", group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = xenium.gw15_palette, coord.fixed = T, molecules = features, mols.size = 0.1, mols.cols = c("#00FF00", "#FFFF00", "#FF0088"), nmols = 10000) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
ggsave('Supplementary_Fig_3A_GW_15.pdf')

# Generate spatial plots for each gene; repeat for each field of view and gestational week

for (feature in features) {
  # Generate the plot for the current feature
  p <- ImageFeaturePlot(xenium.gw15, fov = "zoom", features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature
  filename <- paste0("Supplementary_Fig_3A_GW15_", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}

## Repeat above code for Supplementary Fig 3C however changing the coordinates and features

##GW 15

## Coord 5 - Basal, club, ciliated cells
GW15_coords5 <- Crop(xenium.gw15[["fov"]], x = c(8000, 8300), y = c(8700, 9000), coords = "plot")
# features <- c('TP63', 'KRT5', 'KRT14')
# features <- c('SCGB1A1')
# features <- c('FOXJ1', 'TP73')


## Coord 7 - Budtip & tip
GW15_coords7_budtip <- Crop(xenium.gw15[["fov"]], x = c(11700, 11970), y = c(7800, 8100), coords = "plot")
# features <- c('ETV5', 'CA2', 'SFTPC')

## GW 18

## Coordinates 12 - Basal, club, ciliated
```{r}
GW18_coords12 <- Crop(xenium.gw18[["fov"]], x = c(10225, 10525), y = c(4200, 4500), coords = "plot")
# features <- c('TP63', 'KRT5', 'KRT14')
# features <- c('SCGB1A1')
# features <- c('FOXJ1', 'TP73')

## Coord 14 - Budtip prog & tip
GW18_coords14_budtip <- Crop(xenium.gw18[["fov"]], x = c(9600, 9900), y = c(4050, 4350), coords = "plot")
# features <- c('ETV5', 'CA2', 'SFTPC')
