# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load epithelial dataset & palette
epi_fetal<-readRDS('data/epi_fetal_lung.rds')
epi_palette<-readRDS('palettes/epi_palette.rds')

