# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load stromal dataset & palette
epi_fetal<-readRDS('epi_fetal_lung.rds')
epi_palette<-readRDS('epi_palette.rds')

## Fig 4A
p1<-DimPlot(epi_fetal, cols = epi_palette, label = F)
p1$data$num_ident<-epi_fetal@meta.data$num_ident
LabelClusters(p1, id ='num_ident')

ggsave("Fig_4A.pdf", width = 15, height = 9)

## Fig 4B

## Fig 2C

# Load DEG table for epithelial cluster

DEG_markers<-read.csv('c2_markers.csv')
top3<- DEG_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
features = top3$gene

# DEG features that were chosen
chosen_features = c('TCF21', 'PLEKHH2', 'FOS', 'EGR1', 'MKI67', 'TIMP3', 'LGALS3', 'MYH11', 'TMEM158')
stromal_features <- c(features, chosen_features)

# Change the order of the features to match up with the figure; or can load in the features already ordered
ordered_stromal_features <- readRDS('ordered_epithelial_features.rds')

# Generate heatmap
DoHeatmap(str_fetal, assay = 'RNA', features = ordered_stromal_features, size = 4, angle = 90) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30))
ggsave("Fig_2C.pdf", width = 40, height = 18)

