# Code to generate figures for Fig. 1

rm(list = ls())

# Load packages
library(Seurat)
library(tidyverse)
library(viridis)

# Load integrated fetal lung dataset object, differential genes, & palette

all_fetal <- readRDS('full_fetal_lung_dataset.rds')
all_markers<- read.csv('all_markers_res_0.01.csv')
all_fetal_palette<-readRDS('all_fetal_palette.rds')

# Fig 1B

DimPlot(all_fetal, group.by = 'cell_type')
ggsave("~/Fig_1B.pdf", width = 11.5, height = 8.5)

# Fig 1C

DimPlot(all_fetal, group.by = 'all_cell_type', cols = all_fetal_palette)
ggsave("~/Fig_1C.pdf", width = 11.5, height = 8.5)

# Fig 1D

top5<- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top_5_features = top5$gene
all_fetal <-ScaleData(all_fetal, features = top_5_features)

DoHeatmap(all_fetal, assay = 'RNA', features = top_5_features, size = 4, angle = 90) +
scale_fill_viridis(option = "D") + guides(color = "none")+ theme(axis.title = element_text(size=30)) +theme(axis.text.y = element_text(size = 30)) +
scale_y_discrete(breaks=top_5_features[seq(1,length(top_5_features),by=2)])

ggsave("~/Fig_1D.pdf", width = 40, height = 18)

# Fig 1E

pt <- table(all_fetal$cell_type, all_fetal$sample_week)
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

# Fig 1F

# Load integrated fetal dataset (including publically available datasets)

public_quach_integrated_fetal_datasets<-readRDS('/data/public_quach_integrated_fetal_datsets.rds')

DimPlot(public_quach_integrated_fetal_datasets, cells.highlight=WhichCells(public_quach_integrated_fetal_datasets, idents = c("Quach et al. 2023")) + scale_color_manual(labels = c("Publicly available dataset", "Quach et al. 2023"), values = c('grey', 'red')) + labs(color = "Fetal lung dataset")

ggsave("~/Fig_1F.pdf", width = 11, height = 8.5)

# Fig 1G

# Load publically available integrated fetal dataset

public_fetal_dataset<- readRDS('public_fetal_dataset.rds')

# Find anchors between query and reference dataset

public_quach.anchors <- FindTransferAnchors(reference = all_fetal, query = public_fetal_dataset, dims = 1:30, reference.reduction = "pca")

predictions <- TransferData(anchorset = public_quach.anchors, refdata = all_fetal$all_cell_type, dims = 1:30)
public_fetal_dataset <- AddMetaData(public_fetal_dataset, metadata = predictions)

saveRDS(public_fetal_dataset, 'public_fetal_dataset.rds')

# Generate frequency table of predicted ids

write.csv(table(public_fetal_dataset$predicted.id), 'public_fetal_dataset_query.csv')


