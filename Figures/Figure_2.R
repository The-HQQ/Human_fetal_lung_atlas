# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)
library(AUCell)

# Load stromal dataset & palette
str_fetal<-readRDS('data/stromal_fetal_lung.rds')
str_palette<-readRDS('data/stromal_palette.rds')

## Fig 2A

p1<-DimPlot(str_fetal, cols = str_palette)
p1$data$num_ident<-str_fetal@meta.data$num_ident
LabelClusters(p1, id ='num_ident')

ggsave("Fig_2A.pdf", width = 15, height = 9)

## Fig 2B

pt <- table(str_fetal$cell_type, str_fetal$sample_week)
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
    scale_fill_manual(values = str_palette) +
    NoLegend()
ggsave("Fig_2B.pdf", width = 15, height = 9)

## Fig 2C

# Load DEG table for stromal cluster

DEG_markers<-read.csv('c1_markers.csv')
top3<- all_markers_res_0_01 %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
features = top3$gene

# DEG features that were chosen
chosen_features = c('TCF21', 'PLEKHH2', 'FOS', 'EGR1', 'MKI67', 'TIMP3', 'LGALS3', 'MYH11', 'TMEM158')
stromal_features <- c(features, chosen_features)

# Change the order of the features to match up with the figure; or can load in the features already ordered
ordered_stromal_features <- readRDS('ordered_stromal_features.rds')

# Generate heatmap
DoHeatmap(str_fetal, assay = 'RNA', features = ordered_stromal_features, size = 4, angle = 90) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30))
ggsave("Fig_2C.pdf", width = 40, height = 18)

## Fig 2D

# Load in output loom files from pyscenic & list of top transcription factors
str_loom <- open_loom('data/str_scenic_integrated-output.loom')
str_top_tf <- readRDS('data/str_top_tf.rds')

# Read information from output loom files
regulonAUC <- get_regulons_AUC(str_loom, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(str_loom)
embeddings <- get_embeddings(str_loom)
cells<-regulonAUC@colData@rownames

str_fetal<-subset(str_fetal, cells = cells)
cellInfo<-data.frame(CellType=Idents(str_fetal))

# Sort regulons by top transcription factors
regulonAUC <- regulonAUC[rownames(regulonAUC) %in% str_top_tf,]

# Generate heatmap based on regulon activity

row_names <- levels(str_fetal)
regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity", row_order = c(1:19), column_order = str_top_tf,row_labels = row_names, row_names_side = 'left')

ggsave("Fig_2D.pdf", width = 40, height = 18)






