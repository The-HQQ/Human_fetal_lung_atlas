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






