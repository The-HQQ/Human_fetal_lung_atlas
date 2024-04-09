# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load stromal dataset & palette
epi_fetal<-readRDS('epi_fetal_lung.rds')
epi_palette<-readRDS('epi_palette.rds')

# Supp Fig. 4A

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
