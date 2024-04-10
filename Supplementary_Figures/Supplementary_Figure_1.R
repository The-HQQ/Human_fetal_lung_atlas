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
clustree(xenium.gw15)

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

DotPlot(object = xenium.gw15, features=features, cols = "RdBu") +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16)) + scale_size(range = c(2, 5)) +
    guides(size=guide_legend(title = 'Percent Expressed',override.aes=list(shape=21, colour="black", fill="white")))

ggsave("Supplementary_Fig_1G.pdf")
# ggsave("Supplementary_Fig_1I.pdf")

## Supplementary Fig 1J,K

# Repeat with gw 18

ImageDimPlot(xenium.gw15, fov = "fov", cols = xenium.gw15_palette, axes = TRUE)

ggsave("Supplementary_Fig_1J.pdf")
# ggsave("Supplementary_Fig_1K.pdf")

