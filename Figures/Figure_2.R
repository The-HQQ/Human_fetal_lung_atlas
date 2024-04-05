# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load stromal dataset & palette
str_fetal<-readRDS('stromal_fetal_lung.rds')
str_palette<-readRDS('stromal_palette.rds')

# Fig 2A

DimPlot(str_fetal, cols = str_palette)

# Fig 2B

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
    NoLegend()
ggsave("~/Fig_1E.pdf", width = 15, height = 9)






