# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)
library(SCopeLoomR)
library(SCENIC)

# Load endothelial dataset and palette
endo_fetal<-readRDS('endo_fetal_lung.rds')
endo_palette<-readRDS('endo_palette.rds')

## Supplementary Note Fig 2A

p1<-DimPlot(endo_fetal, cols = endo_palette)
p1$data$num_ident<-endo_fetal@meta.data$num_ident
LabelClusters(p1, id ='num_ident')

ggsave("Supplementary_Note_Fig_2A.pdf", width = 15, height = 9)

## Supplementary Note Fig 2B

pt <- table(endo_fetal$cell_type, endo_fetal$sample_week)
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
    scale_fill_manual(values = endo_palette) +
    NoLegend()
ggsave("Supplementary_Note_Fig_2B.pdf", width = 15, height = 9)

## Supplementary Note Fig 2C

# Load DEG table for endothelial cluster

DEG_markers<-read.csv('endo_markers.csv')
top3<- DEG_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
features = top3$gene

# DEG features that were chosen
chosen_features = c('NRP2', 'TOP2A', 'MKI67', 'CENPF', 'EGFL6', 'COL6A3', 'REL', 'IRF', 'POU2F2', 'CPE', 'APLNR')
endo_features <- c(features, chosen_features)

# Change the order of the features to match up with the figure; or can load in the features already ordered
ordered_endo_features <- readRDS('ordered_endo_features.rds')

# Generate heatmap
DoHeatmap(endo_fetal, assay = 'RNA', features = ordered_endo_features, size = 4, angle = 90) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30))
ggsave("Supplementary_Note_Fig_2C.pdf", width = 40, height = 18)

## Supplementary Note Fig 2D

# Load in output loom files from pyscenic & list of top transcription factors
endo_loom <- open_loom('c3_scenic_integrated-output.loom')
endo_top_tf <- readRDS('c3_top_tf.rds')

# Read information from output loom files
regulonAUC <- get_regulons_AUC(endo_loom, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(endo_loom)
embeddings <- get_embeddings(endo_loom)
cells<-regulonAUC@colData@rownames

endo_fetal<-subset(endo_fetal, cells = cells)
cellInfo<-data.frame(CellType=Idents(endo_fetal))
endo_top_tf=paste0(endo_top_tf,"_(+)")

# Sort regulons by top transcription factors
regulonAUC <- regulonAUC[rownames(regulonAUC) %in% endo_top_tf,]

# Generate heatmap based on regulon activity

row_names <- levels(endo_fetal)
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity", row_order = c(1:10), column_order = endo_top_tf,row_labels = row_names, row_names_side = 'left')

ggsave("Supplementary_Note_Fig_2D.pdf", width = 40, height = 18)

## Supplementary Note Fig 2H

# Load in processed xenium objects and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

options(future.globals.maxSize = 3000 * 1024^2)

# Select coordinates to zoom

GW15_coords8 <- Crop(xenium.gw15[["fov"]], x = c(12400, 12700), y = c(6580, 6880), coords = "plot")

xenium.gw15[["zoom8"]] <- GW15_coords8
DefaultBoundary(xenium.gw15[["zoom8"]]) <- "segmentation"

GW18_coords8<- Crop(xenium.gw18[["fov"]], x = c(12300, 12600), y = c(3100, 3400), coords = "plot")

xenium.gw18[["zoom8"]] <- GW18_coords8
DefaultBoundary(xenium.gw18[["zoom8"]]) <- "segmentation"

# Define features to plot

#features <- c('CA4', 'KDR', 'STC1', 'THY1')

# Define gestational week and zoom to plot

# xenium.obj<- xenium.gw15
# xenium.obj<- xenium.gw18
# colours <- xenium.gw15_palette
# colours <- xenium.gw18_palette

zoom<-"zoom8"

# Plot spatial dimplot

ImageDimPlot(xenium.obj, fov = zoom, group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Supplementary_Note_Fig_2H_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature; change
  p <- ImageFeaturePlot(xenium.obj, fov = zoom, features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Supplementary_Note_Figure_2H_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}
