library(SeuratDisk)
library(Seurat)
library(glue)
library(tidyr)
library(dplyr)
library(magrittr)
library(clustree)
library(pheatmap)
library(purrr)


indir <- "~/projects/seq3/20230530/4_scrna/counts/6_rename_cluster/"

# obj.genes = readRDS(glue("{indir}/tgs-sc_genes/RIF.cluster.celltype_fine.rds"))
# obj.trans = readRDS(glue("{indir}/tgs-sc_isoforms/RIF.cluster.celltype_fine.rds"))
# 
# # #common_cells = intersect(Cells(obj.genes), Cells(obj.trans))
# common_cells = intersect(
#   rownames(obj.genes@meta.data %>% filter(grepl('^eStr|pStr|VSMC', celltype_fine))),
#   rownames(obj.trans@meta.data %>% filter(grepl('^eStr|pro|VSMC', celltype_fine)))
# )
# # common_cells = intersect(
# #   Cells(obj.genes),
# #   rownames(obj.trans@meta.data %>% filter(grepl('^eStr|pro|VSMC', celltype_fine)))
# # )
# obj.genes.common <- subset(obj.genes, cells = common_cells) %>%
#  FindNeighbors(dims = 1:20) %>%
#  FindClusters(resolution = 0.8) %>%
#  RunUMAP(dims = 1:20)
# 
# # cluster
# condition <- sapply(obj.genes.common@meta.data$orig.ident, function(x) ifelse(grepl("Pro",x), "pro", "sec"))
# p1 <- DimPlot(obj.genes.common, group.by = "celltype_fine", label = T) & NoAxes()
# p2 <- DimPlot(obj.genes.common, group.by = "seurat_clusters", label = T) & NoAxes()
# p1 + p2
# 
# markers.genes <- Seurat::FindAllMarkers(obj.genes.common, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# markers.genes.top50 <- markers.genes %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC) %>% distinct(gene) %>% .$gene
# # correlation by avg expression
# avg <- AverageExpression(obj.genes.common, group.by = "seurat_clusters")
# avg <- avg$RNA
# #cg <- names(tail(sort(apply(avg, 1, sd)), 1000))
# p3 <- pheatmap(cor(avg[markers.genes.top50,], method = 'spearman'), main = 'expression')
# # correlation by avg fold change
# data <- map_dfc(levels(obj.genes.common$seurat_clusters),
#                 ~ FoldChange(obj.genes.common, ident.1 = ., features = markers.genes.top50) %>% select(avg_log2FC))
# names(data) <- levels(obj.genes.common$seurat_clusters)
# p4 <- pheatmap(cor(data, method = 'spearman'), main = 'fold change')
# grid.arrange(p1, p2, p3[[4]], p4[[4]], nrow = 2, heights = c(3,2))
# 
# # cells in each clusters
# table(obj.genes.common$seurat_clusters)
# 
# clusters.gene.fine <-  as.character(obj.genes.common$seurat_clusters)
# clusters.gene.fine[clusters.gene.fine %in% c(6)] <- "pStr"
# clusters.gene.fine[clusters.gene.fine %in% c(7)] <- "VSMC"
# clusters.gene.fine[clusters.gene.fine %in% c(0,3)] <- "gStr1"
# clusters.gene.fine[clusters.gene.fine %in% c(1,5)] <- "gStr2"
# clusters.gene.fine[clusters.gene.fine %in% c(2)] <- "gStr3"
# clusters.gene.fine[clusters.gene.fine %in% c(4)] <- "gStr4"
# 
# 
# 
# obj.trans.common <- subset(obj.trans, cells = common_cells) %>%
#  FindNeighbors(dims = 1:20) %>%
#  FindClusters(resolution = 0.8) %>%
#  RunUMAP(dims = 1:20)
# # cluster
# condition <- sapply(obj.trans.common@meta.data$orig.ident, function(x) ifelse(grepl("Pro",x), "pro", "sec"))
# p1 <- DimPlot(obj.trans.common, group.by = "celltype_fine", label = T) & NoAxes()
# p2 <- DimPlot(obj.trans.common, group.by = "seurat_clusters", label = T) & NoAxes()
# p1 + p2
# 
# markers.genes <- Seurat::FindAllMarkers(obj.trans.common, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# markers.genes.top50 <- markers.genes %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC) %>% distinct(gene) %>% .$gene
# # correlation by avg expression
# avg <- AverageExpression(obj.trans.common, group.by = "seurat_clusters")
# avg_raw <- avg$RNA
# #cg <- names(tail(sort(apply(avg, 1, sd)), 1000))
# p3 <- pheatmap(cor(avg_raw[markers.genes.top50,], method = 'spearman'), main = 'raw expression')
# avg_debatch <- avg$integrated
# p4 <- pheatmap(cor(avg_debatch[markers.genes.top50,], method = 'spearman'), main = 'debatch expression')
# # correlation by avg fold change
# data <- map_dfc(levels(obj.trans.common$seurat_clusters),
#                 ~ FoldChange(obj.trans.common, ident.1 = ., features = markers.genes.top50) %>% select(avg_log2FC))
# names(data) <- levels(obj.trans.common$seurat_clusters)
# p5 <- pheatmap(cor(data, method = 'spearman'), main = 'fold change')
# grid.arrange(p1, p2, p3[[4]], p4[[4]], p5[[4]], nrow = 2, heights = c(3,2), 
#              layout_matrix = rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5)))
# 
# table(obj.trans.common$seurat_clusters)
# 
# # show highlight clusters
# DimPlot(obj.trans.common, label= T, label.size = 10,
#         cells.highlight = Seurat::CellsByIdentities(obj.trans.common, idents = c(5)))
# 
# 
# resolution = 0.8
# clusters.isoform.fine <-  as.character(obj.trans.common$seurat_clusters)
# clusters.isoform.fine[clusters.isoform.fine %in% c(7,11)] <- "pro"
# clusters.isoform.fine[clusters.isoform.fine %in% c(8)] <- "VSMC"
# clusters.isoform.fine[clusters.isoform.fine %in% c(0)] <- "tStr1"
# clusters.isoform.fine[clusters.isoform.fine %in% c(4)] <- "tStr2"
# clusters.isoform.fine[clusters.isoform.fine %in% c(3,10)] <- "tStr3"
# clusters.isoform.fine[clusters.isoform.fine %in% c(1)] <- "tStr4"
# clusters.isoform.fine[clusters.isoform.fine %in% c(2)] <- "tStr5"
# clusters.isoform.fine[clusters.isoform.fine %in% c(5)] <- "tStr6"
# clusters.isoform.fine[clusters.isoform.fine %in% c(6,9)] <- "tStr7"
# clusters.isoform.fine[clusters.isoform.fine %in% c(12)] <- "tStr8"
# 
# 
# obj.trans.common %>%
#   AddMetaData(clusters.isoform.fine, col.name = "str_trans_fine") %>%
#   AddMetaData(clusters.gene.fine, col.name = "str_genes_fine") %>%
#   AddMetaData(condition , col.name = "condition") %>%
#   AddMetaData(obj.genes.common@meta.data$celltype_fine, col.name = "celltype_genes") %>%
#   AddMetaData(obj.trans.common@meta.data$celltype_fine, col.name = "celltype_trans") -> object.isoforms
# p6 <- DimPlot(object.isoforms, group.by = "str_trans_fine", label = T) & NoAxes()
# grid.arrange(p1, p2, p6, p3[[4]], p4[[4]], p5[[4]], nrow = 2, heights = c(3,2), 
#              layout_matrix = rbind(c(1,2,3), c(4,5,6)))
# saveRDS(object.isoforms, glue("{indir}/RIF.stromal.trans.rds"))
# 
# obj.genes.common %>%
#   AddMetaData(clusters.isoform.fine, col.name = "str_trans_fine") %>%
#   AddMetaData(clusters.gene.fine, col.name = "str_genes_fine") %>%
#   AddMetaData(condition , col.name = "condition") %>%
#   AddMetaData(obj.genes.common@meta.data$celltype_fine, col.name = "celltype_genes") %>%
#   AddMetaData(obj.trans.common@meta.data$celltype_fine, col.name = "celltype_trans") -> object.genes
# saveRDS(object.genes, glue("{indir}/RIF.stromal.genes.rds"))

# DimPlot(object.isoforms, group.by = "str_trans_fine", label = T)
# DimPlot(object.isoforms, group.by = "str_genes_fine", label = T)
# DimPlot(object.isoforms, group.by = "seurat_clusters", label = T)
# Seurat::DimPlot(object.isoforms, label = T, label.size = 5,
#                 cells.highlight = object.isoforms@meta.data %>% filter(str_trans_fine == 'eStr5') %>% rownames)
# any(Cells(obj.genes.common) != Cells(obj.genes.common))
# any(row.names(obj.genes.common@meta.data) != row.names(obj.trans.common@meta.data))

object.genes <- readRDS(glue("{indir}/RIF.stromal.genes.rds"))
object.genes@assays$integrated@counts = object.genes@assays$RNA@counts
out = glue("{indir}/RIF.stromal.genes.h5Seurat")
SeuratDisk::SaveH5Seurat(object.genes, filename = out)
SeuratDisk::Convert(out, dest = "h5ad", overwrite=TRUE)


object.isoforms <- readRDS(glue("{indir}/RIF.stromal.trans.rds"))
object.isoforms@assays$integrated@counts = object.isoforms@assays$RNA@counts
out = glue("{indir}/RIF.stromal.trans.h5Seurat")
SeuratDisk::SaveH5Seurat(object.isoforms, filename = out)
SeuratDisk::Convert(out, dest = "h5ad", overwrite=TRUE)
