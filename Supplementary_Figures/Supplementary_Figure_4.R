# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)
library(SCopeLoomR)
library(SCENIC)

# Load stromal dataset & palette
epi_fetal<-readRDS('epi_fetal_lung.rds')
epi_palette<-readRDS('epi_palette.rds')

## Supplementary Fig. 4A

pt <- table(epi_fetal$cell_type, epi_fetal$sample_week)
pt <- as.data.frame(pt)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(colour = "black", position = "fill") +
    xlab("Sample") +
    ylab("Proportion") +
    theme(legend.title = element_blank())+
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text = element_text(size=25), axis.title = element_text(size=25)) +
    theme(legend.text=element_text(size=20)) +
    scale_fill_manual(values = epi_palette) +
    NoLegend()

## Supplementary Fig. 4B

# Load in output loom files from pyscenic & list of top transcription factors
epi_loom <- open_loom('c2_scenic_integrated-output.loom')
epi_top_tf <- readRDS('c2_top_tf.rds')

# Read information from output loom files
regulonAUC <- get_regulons_AUC(epi_loom, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(epi_loom)
embeddings <- get_embeddings(epi_loom)
cells<-regulonAUC@colData@rownames

epi_fetal<-subset(epi_fetal, cells = cells)
cellInfo<-data.frame(CellType=Idents(epi_fetal))
epi_top_tf=paste0(epi_top_tf,"_(+)")

# Sort regulons by top transcription factors
regulonAUC <- regulonAUC[rownames(regulonAUC) %in% epi_top_tf,]

# Generate heatmap based on regulon activity

row_names <- levels(epi_fetal)
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity", row_order = c(1:19), column_order = epi_top_tf,row_labels = row_names, row_names_side = 'left')

ggsave("Supplementary_Fig_2B.pdf", width = 40, height = 18)
## Supplementary Fig. 4D

# Load in processed xenium objects and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

options(future.globals.maxSize = 3000 * 1024^2)

# Select coordinates to zoom
## For GW15
## Coord 6 - Smaller airway
GW15_coords6 <- Crop(xenium.gw15[["fov"]], x = c(12225, 12525), y = c(9125, 9425), coords = "plot")
xenium.gw15[["zoom"]] <- GW15_coords6
DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"
                                     
## For GW18
## Coordinates 13 - Smaller airway
GW18_coords13 <- Crop(xenium.gw18[["fov"]], x = c(13750, 14050), y = c(1650, 1950), coords = "plot")
xenium.gw18[["zoom"]] <- GW15_coords13
DefaultBoundary(xenium.gw15[["zoom"]]) <- "segmentation"                                     
# Switch when needed
# xenium.obj<-xenium.gw15 
# xenium.obj<-xenium.gw18                                      
# Plot spatial dimplot

features <- c('TP63', 'SCGB3A2', 'FOXJ1')
ImageDimPlot(xenium.obj, fov = "zoom", group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Supplementary_Fig_4D_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature; change
  p <- ImageFeaturePlot(xenium.obj, fov = "zoom", features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Supplementary_Figure_4D_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}
                                     
