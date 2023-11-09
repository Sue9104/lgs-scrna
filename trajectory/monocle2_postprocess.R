library(monocle)
library(Seurat)
library(tidyr)
library(dplyr)
outdir <- "/public/home/msu/projects/seq3/20230530/5_trajectory/monocle2"

{
  prefix <- "hvg_seurat"
  cds <- readRDS(glue::glue("{outdir}/cds_{prefix}.rds"))
  p <- plot_ordering_genes(cds) + ggtitle(prefix)
  print(p)
  
  # plot monocle results
  p <- plot_cell_trajectory(cds, color_by = "State") + ggtitle(prefix) 
  print(p)
  p <- plot_cell_trajectory(cds, color_by = "CellType") + ggtitle(prefix)
  print(p)
  p <- plot_cell_trajectory(cds, color_by = "Pseudotime") + ggtitle(prefix)
  print(p)
  p <- plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 1)
  print(p)
  p <- plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~CellType, nrow = 2)
  print(p)
  p <- plot_cell_trajectory(cds, color_by = "CellType") + facet_wrap(~State, nrow = 2)
  print(p)
}

