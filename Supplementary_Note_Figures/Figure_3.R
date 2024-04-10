# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(enrichR)
library(tidyverse)
library(viridis)
library(SCopeLoomR)
library(SCENIC)

# Load immune dataset and palette
immune_fetal<-readRDS('immune_fetal_lung.rds')
immune_palette<-readRDS('immune_palette.rds')

## Supplementary Note Fig 3A

p1<-DimPlot(immune_fetal, cols = immune_palette)
p1$data$num_ident<-immune_fetal@meta.data$num_ident
LabelClusters(p1, id ='num_ident')

ggsave("Supplementary_Note_Fig_3A.pdf", width = 15, height = 9)

## Supplementary Note Fig 3B

pt <- table(immune_fetal$cell_type, immune_fetal$sample_week)
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
    scale_fill_manual(values = immune_palette) +
    NoLegend()
ggsave("Supplementary_Note_Fig_3B.pdf", width = 15, height = 9)

## Supplementary Note Fig 3C

# Load DEG table for immune cluster

DEG_markers<-read.csv('immune_markers.csv')
top3<- DEG_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
features = top3$gene

# DEG features that were chosen
chosen_features = c('KLF11', 'SMAD6', 'IL1B','IL10','TNF', 'CD68','FLT3', 'CD1C')
immune_features <- c(features, chosen_features)

# Change the order of the features to match up with the figure; or can load in the features already ordered
ordered_immune_features <- readRDS('ordered_immune_features.rds')

# Generate heatmap
DoHeatmap(immune_fetal, assay = 'RNA', features = ordered_immune_features, size = 4, angle = 90) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30))
ggsave("Supplementary_Note_Fig_3C.pdf", width = 40, height = 18)

## Supplementary Note Fig 3D

for (i in 1:length(levels(immune_fetal))){
p1<-DEenrichRPlot(immune_fetal, ident.1 = levels(immune_fetal)[i], ident.2 = levels(immune_fetal)[-i],   enrich.database = "GO_Biological_Process_2023", max.genes = 100, balanced = F)
p1$layers[[2]]<-NULL
p1 + geom_text(aes_string(label = "term", y = 0), 
               size = 8, color = "black", position = position_dodge(1), 
               hjust = 0)
ggsave(paste(levels(immune_fetal)[i], "GO_Bio_pro_2023.png", sep = ""), height = 6, width = 13, units = 'in', dpi = 300)

p2<-DEenrichRPlot(immune_fetal, ident.1 = levels(immune_fetal)[i], ident.2 = levels(immune_fetal)[-i],   enrich.database = "GO_Molecular_Function_2023", max.genes = 100, balanced = F)
p2$layers[[2]]<-NULL
p2 + geom_text(aes_string(label = "term", y = 0), 
               size = 8, color = "black", position = position_dodge(1), 
               hjust = 0)
ggsave(paste(levels(immune_fetal)[i], "GO_Mol_fun_2023.png", sep = ""), height = 6, width = 13, units = 'in', dpi = 300)

p3<-DEenrichRPlot(immune_fetal, ident.1 = levels(immune_fetal)[i], ident.2 = levels(immune_fetal)[-i],   enrich.database = "GO_Cellular_Component_2023", max.genes = 100, balanced = F)
p3$layers[[2]]<-NULL
p3 + geom_text(aes_string(label = "term", y = 0), 
               size = 8, color = "black", position = position_dodge(1), 
               hjust = 0)
ggsave(paste(levels(immune_fetal)[i], "GO_Cell_comp_2023.png", sep = ""), height = 6, width = 13, units = 'in', dpi = 300)
}

## Supplementary Note Fig 3E

# Load in output loom files from pyscenic & list of top transcription factors
immune_loom <- open_loom('immune_scenic_integrated-output.loom')
immune_top_tf <- readRDS('immune_top_tf.rds')

# Read information from output loom files
regulonAUC <- get_regulons_AUC(immune_loom, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(immune_loom)
embeddings <- get_embeddings(immune_loom)
cells<-regulonAUC@colData@rownames

immune_fetal<-subset(immune_fetal, cells = cells)
cellInfo<-data.frame(CellType=Idents(immune_fetal))
immune_top_tf=paste0(immune_top_tf,"_(+)")

# Sort regulons by top transcription factors
regulonAUC <- regulonAUC[rownames(regulonAUC) %in% immune_top_tf,]

# Generate heatmap based on regulon activity

row_names <- levels(immune_fetal)
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
                                     

ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity", row_order = c(1:18), column_order = immune_top_tf,row_labels = row_names, row_names_side = 'left')

ggsave("Supplementary_Note_Fig_3E.pdf", width = 40, height = 18)

## Supplementary Note Fig 3F
                                     
# Load in processed xenium objects and palettes

xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')

xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

# Select coordinates to zoom

GW15_coords9 <- Crop(xenium.gw15[["fov"]], x = c(12100, 12400), y = c(6200, 6500), coords = "plot")

xenium.gw15[["zoom9"]] <- GW15_coords9
DefaultBoundary(xenium.gw15[["zoom9"]]) <- "segmentation"

GW18_coords9<- Crop(xenium.gw18[["fov"]], x = c(11200, 11500), y = c(6700, 7000), coords = "plot")

xenium.gw18[["zoom9"]] <- GW18_coords9
DefaultBoundary(xenium.gw18[["zoom9"]]) <- "segmentation"

# Define features to plot

#features <- c('CD86', 'CD68', 'FCN1', 'CSTA')

# Define gestational week and zoom to plot

# xenium.obj<- xenium.gw15
# xenium.obj<- xenium.gw18
# colours <- xenium.gw15_palette
# colours <- xenium.gw18_palette

zoom<-"zoom9"

# Plot spatial dimplot

ImageDimPlot(xenium.obj, fov = zoom, group.by = "cell_type", axes = TRUE, border.color = "white", border.size = 0.075, cols = colours, coord.fixed = T, molecules = features, mols.size = 0.05, mols.cols = c("#FFFF00", "#00FF00", "#FF0088"), nmols = 10000) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define the filename using the current feature, change based on GW
ggsave(paste0("Supplementary_Note_Fig_3F_GW_15_dimplot.pdf"))

# Plot feature plots

for (feature in features) {
  # Generate the plot for the current feature; change
  p <- ImageFeaturePlot(xenium.obj, fov = zoom, features = feature, max.cutoff = 'q90', size = 0.75, axes = FALSE, coord.fixed = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size = 20))
  
  # Define the filename using the current feature, change based on GW
  filename <- paste0("Supplementary_Note_Figure_3F_GW_15", feature, ".pdf")
  
  # Save the plot to a PDF file
  ggsave(filename, plot = p, device = "pdf", height = 4.5, width = 5, units = "in")
}
                                     
                                  
