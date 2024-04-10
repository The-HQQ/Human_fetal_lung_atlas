# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(clustree)

# Load integrated fetal dataset (including publically available datasets), Quach & Farrell et al. dataset, xenium datasets, and xenium palettes

all_fetal <- readRDS('data/full_fetal_lung_dataset.rds')
public_quach_integrated_fetal_datasets<-readRDS('public_quach_integrated_fetal_datsets.rds')
xenium.gw15 <- readRDS('data/xenium_gw15_processed.rds')
xenium.gw18 <-readRDS('data/xenium_gw18_processed.rds')
xenium.gw15_palette<-readRDS('palettes/xenium_gw15_palette.rds')
xenium.gw18_palette<-readRDS('palettes/xenium_gw18_palette.rds')

## Supplementary Fig 1B

clustree(all_fetal)

ggsave("Supplementary_Fig_1B.pdf")

## Supplementary Fig 1C

DimPlot(public_quach_integrated_fetal_datasets, group.by = 'dataset')

ggsave("Supplementary_Fig_1C.pdf")

## Supplementary Fig 1D,E

# Repeat with gw 18 
# Set other resolutions to null to save space for figure
xenium.gw15@meta.data$SCT_snn_res.0.3 <- NULL
xenium.gw15@meta.data$SCT_snn_res.0.5 <- NULL
xenium.gw15@meta.data$SCT_snn_res.1 <- NULL

clustree(xenium.gw15, prefix = 'SCT_snn_res.')

ggsave("Supplementary_Fig_1D.pdf")
# ggsave("Supplementary_Fig_1E.pdf")

## Supplementary Fig 1F,H

# Repeat with gw 18
DimPlot(xenium.gw15, cols = xenium.gw15_palette)

ggsave("Supplementary_Fig_1F.pdf")
# ggsave("Supplementary_Fig_1H.pdf")

## Supplementary Fig 1G,I

features <- c('KDR', 'THY1', 'EPCAM', 'NKX2-1', 'CD86', 'CD68', 'COL5A2', 'TCF21', 'ACTA2', 'FBN1', 'WNT2', 'FGFR4')

# Repeat with gw 18

DotPlot(xenium.gw15, features = features) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  geom_point(aes(size=pct.exp),shape=21, color="black", stroke=0.5) +
  theme(axis.text.x = element_text(angle = 90))

ggsave("Supplementary_Fig_1G.pdf", width = 7.3, height = 3.5)
# ggsave("Supplementary_Fig_1I.pdf", , width = 7.3, height = 3.5)

## Supplementary Fig 1J,K

# Repeat with gw 18

ImageDimPlot(xenium.gw15, fov = "fov", group.by = "cell_type", axes = TRUE, cols = colours, coord.fixed = T, size = 0.1, dark.background = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none")

ggsave("Supplementary_Fig_1J.pdf", dpi =1000)
# ggsave("Supplementary_Fig_1K.pdf", dpi = 1000)

