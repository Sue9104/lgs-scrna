library(Seurat)
library(decoupleR)
library(viper)
# library(dorothea)
# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(magrittr)

# input
obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.20230906.rds")
obj.trans <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.no_dup.rds")
com_cells <- intersect(Cells(obj.genes), Cells(obj.trans))
obj.genes_com <- subset(obj.genes, cells = com_cells) 
obj.trans_com <- subset(obj.trans, cells = com_cells) 

# assay = "integrated"
assay = "RNA"
mat <- as.matrix(obj.genes_com@assays[[assay]]@data)
data <- obj.trans_com
outdir <- "~/projects/seq3/20230530/8_regulons/dorothea"

###############################################################
##### PROGENy: curated collection of pathways and target genes
##### Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB
##### p53, POI3K, TGFb, TNFa, Trail, VEGF, WNT
###############################################################
progeny = get_progeny(organism='human')
# progeny = get_progeny(organism='human', top = 10)
progeny %>% head

curated_acts <- run_mlm(mat=mat, net=progeny, .source='source', .target='target',
                       .mor='weight', minsize = 5)
curated_acts
# create new assay
data[['progeny']] <- curated_acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = data) <- "progeny"
data <- ScaleData(data)
data@assays$progeny@data <- data@assays$progeny@scale.data
saveRDS(data, glue::glue("{outdir}/RIF.trans.progeny.{assay}.rds"))


# # top pathways
data_pro <- subset(data, subset = condition == "Proliferative")
data_sec <- subset(data, subset = condition == "Secretory")
plot_pathway_heatmap <- function(obj, n_top = 30, assay = "progeny"){
  # obj <- data_pro
  # n_top <- 20
  # assay <- "progeny"
  pathways.df <- t(as.matrix(obj@assays[[assay]]@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(obj)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    # to avoid lower mean in case of clusters with large cells
    ## top 100 cells
    # slice_max(order_by = abs(score), n = 100) %>%
    ## or top 10%
    slice_max(order_by = abs(score), prop = 0.1) %>%
    summarise(mean = mean(score)) %>% 
    pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
    column_to_rownames('cluster') %>% as.matrix()
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  pheatmap(pathways.df, border_color = NA, color=my_color, breaks = my_breaks,
           cluster_rows = F,clustering_distance_cols = "correlation")
}
plot_pathway_heatmap(data) 
plot_pathway_heatmap(data_pro)
plot_pathway_heatmap(data_sec)

# # plot regulon activity
p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle('Cell types')
p1
p2 <- (FeaturePlot(data_sec, features = pathways, ncol = 5) &
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) 
p2


# plot specific pathway
pathway <- 'MAPK'
pathway_progeny <- progeny %>%
  filter(source == pathway) %>%
  arrange(target) %>%
  mutate(ID = target, color = "black") %>%
  column_to_rownames('target')
inter <- sort(intersect(rownames(mat),rownames(pathway_progeny)))
pathway_progeny <- pathway_progeny[inter, ]
pathway_progeny[['expr']] <- rowMeans(mat[inter,])
pathway_progeny %<>% mutate(color = case_when(expr > 1.5 ~ "blue",
                                              abs(weight) > 5 ~ "red",
                                              .default = "black"))

ggplot(pathway_progeny, aes(x = weight, y = expr, color = color)) + geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  ggrepel::geom_label_repel(aes(label = ID)) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(pathway)

###############################################################
##### dorothEA: TF databases
###############################################################
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A","B","C"))

# Run viper
tf_acts <- decoupleR::run_viper(mat, regulons, 
                     .source='tf', pleiotropy = FALSE)
tf_acts
# create new tf assay
data[['dorothea']] <- tf_acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = data) <- "dorothea"
data <- ScaleData(data)
data@assays$dorothea@data <- data@assays$dorothea@scale.data
saveRDS(data, glue::glue("{outdir}/RIF.trans.dorothea.{assay}.rds"))



data_pro <- subset(data, subset = condition == "Proliferative")
data_sec <- subset(data, subset = condition == "Secretory")
plot_regulons_heatmap <- function(obj, n_top = 50, assay = "dorothea"){
  # obj <- data
  # n_top <- 20
  # assay <- "dorothea"
  regulons.df <- t(as.matrix(obj@assays[[assay]]@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(obj)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>% 
    # to avoid lower mean in case of clusters with large cells
    ## top 100 cells
    slice_max(order_by = abs(score), n = 100) %>%
    ## or top 10%
    # slice_max(order_by = abs(score), prop = 0.1) %>%
    summarise(mean = mean(score))
  # max std
  # tfs <- regulons.df %>% group_by(source) %>% summarise(std = sd(mean)) %>%
  #   arrange(-abs(std)) %>% head(n_top) %>% pull(source)
  # max mean
  tfs <- regulons.df %>% group_by(cluster) %>%
    slice_max(abs(mean), n = 10)  %>% distinct(source) %>% pull(source)
  
  top_acts_mat <- regulons.df %>% filter(source %in% tfs) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
    column_to_rownames('cluster') %>% as.matrix()
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks,
           cluster_rows = F,clustering_distance_cols = "correlation")
}

# 1400 * 400
plot_regulons_heatmap(data)
plot_regulons_heatmap(data_pro)
plot_regulons_heatmap(data_sec)


# # Plot
# pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
# 
# # plot regulon activity
# p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
#   NoLegend() + ggtitle('Cell types')
data_fib <- subset(data, subset = trans_main == "Stromal Fibroblasts")
sign_tfs <- c("FOXP2", "ZEB2", "ZEB1", "STAT4", "FOXL2", "FLI1", "WT1", "SOX9",
              "MAF", "HOXB13", "RFX1", "MEF2A", "MYC", "SP1", "TFAP2C")
p2 <- (FeaturePlot(data_fib, features = sign_tfs, ncol = 5) &
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) 
p2
# DefaultAssay(object = data) <- "RNA"
# p3 <- FeaturePlot(data, features = intersect(rownames(data), tfs), ncol = 5) 
# p1 | p2 | p3

