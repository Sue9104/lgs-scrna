library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(magrittr)
# input
# obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.20230906.rds")
data <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.no_dup.rds")
data_pro <- subset(data, subset = condition == "Proliferative")
data_sec <- subset(data, subset = condition == "Secretory")

indir <- "~/projects/seq3/20230530/8_regulons/trans/hvg_vs_tf_coding/hvt_50000"
auc_mtx  <- read.csv(glue::glue("{indir}/auc_mtx.normal.csv"), row.names = 1)
auc_mtx_z  <- read.csv(glue::glue("{indir}/auc_mtx_Z.csv"), row.names = 1)
bins <- read.csv(glue::glue("{indir}/bin.csv"), row.names = 1)

data[["scenicAUCscore"]] <- Seurat::CreateAssayObject(t(auc_mtx))
data[["scenicAUCZscore"]] <- Seurat::CreateAssayObject(t(auc_mtx_z))
data[["scenicBinary"]] <- Seurat::CreateAssayObject(t(bins))

regulons_Binary_per_celltype <- function(obj, assay = "scenicBinary"){
  t(as.matrix(obj@assays[[assay]]@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(obj)) %>% 
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%    
    summarise(pct = round(sum(score)/n(), 3) ) -> df
  # tfs <- df %>% group_by(source) %>% summarise(std = sd(pct)) %>%
  #   arrange(-abs(std)) %>% head(100) %>% pull(source)
  tfs <- df %>% group_by(cluster) %>% 
    slice_max(pct, n=10) %>% pull(source)
  
  top_acts_mat <- df %>% filter(source %in% tfs) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'pct') %>%
    column_to_rownames('cluster') %>% as.matrix()
  
  pheatmap(top_acts_mat, border_color = NA, color=brewer.pal(n = 9, name = "Blues"), 
           cluster_rows = F)
}

regulons_Binary_per_celltype(data)
regulons_Binary_per_celltype(data_pro)
regulons_Binary_per_celltype(data_sec)



regulons_Z_per_celltype <- function(obj, assay = "scenicAUCZscore"){
  t(as.matrix(obj@assays[[assay]]@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(obj)) %>% 
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%    
    # to avoid lower mean in case of clusters with large cells
    ## top 100 cells
    slice_max(order_by = abs(score), n = 100) %>%
    ## or top 10%
    # slice_max(order_by = abs(score), prop = 0.1) %>%
    summarise(mean = mean(score)) -> df

  # max std
  # tfs <- regulons.df %>% group_by(source) %>% summarise(std = sd(mean)) %>%
  #   arrange(-abs(std)) %>% head(n_top) %>% pull(source)
  # max mean
  tfs <- df %>% group_by(cluster) %>%
    slice_max(abs(mean), n = 10)  %>% distinct(source) %>% pull(source)
  
  top_acts_mat <- df %>% filter(source %in% tfs) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
    column_to_rownames('cluster') %>% as.matrix()
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks,
           cluster_rows = F,clustering_distance_cols = "correlation")
}
regulons_Z_per_celltype(data)
regulons_Z_per_celltype(data_pro)
regulons_Z_per_celltype(data_sec)

# heatmap for all cells
DefaultAssay(data) <- "scenicAUCZscore"
data <- ScaleData(data, features = rownames(data))
markers <- FindAllMarkers(data, logfc.threshold = 0.25, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, abs(pct.1 - pct.2)> 0.2) %>%
  slice_max(abs(avg_log2FC), n = 5) %>%
  ungroup() -> top10
DoHeatmap(data, features = top10$gene, slot="data", size = 6) + NoLegend()

# heatmap for all cells
DefaultAssay(data) <- "scenicBinary"
markers <- FindAllMarkers(data, logfc.threshold = 0.25, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, abs(pct.1 - pct.2)> 0.2) %>%
  slice_max(abs(avg_log2FC), n = 5) %>%
  ungroup() -> top10
DoHeatmap(data, features = top10$gene, slot="data", size = 3) + NoLegend()

