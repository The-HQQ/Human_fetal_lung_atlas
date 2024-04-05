



# Deconvolution
annotations.df <- RCTD_xenium_8055@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium_8055_processed$predicted.celltype <- annotations
keep.cells <- Cells(xenium_8055_processed)[!is.na(xenium_8055_processed$predicted.celltype)]
xenium_8055_processed <- subset(xenium_8055_processed, cells = keep.cells)



options(future.globals.maxSize = 8000 * 1024^2)
GW18_coords1 <- Crop(xenium_8055_processed[["fov"]], x = c(7330, 7630), y = c(3080, 3380), coords = "plot")

xenium.gw18[["coords1"]] <- GW18_coords1
DefaultBoundary(xenium.gw18[["coords1"]]) <- "segmentation"

xenium_8055_coords1<-subset_opt(xenium_8055_processed,cells = Cells(xenium_8055_processed@images$coords1))

celltype.plot <- ImageDimPlot(xenium_8055_coords1, fov = 'coords1', group.by = 'predicted.celltype', size = 2, cols = all_final_palette, border.color = 'white', border.size = 0.1) 

ggsave('~/scratch/8055_xenium_RCTD_coords3.pdf', width = 11, height = 8.5)
