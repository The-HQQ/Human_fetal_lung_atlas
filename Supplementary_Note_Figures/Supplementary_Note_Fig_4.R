# Clear workspace
rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load visium objects
vis.gw18<-readRDS('/data/vis.gw18.rds')
vis.gw15<-readRDS('/data/vis.gw15.rds')

# Define features to plot

features<-
