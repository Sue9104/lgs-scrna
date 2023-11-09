library(monocle)
library(Seurat)
library(tidyr)
library(dplyr)
library(patchwork)
outdir <- "/public/home/msu/projects/seq3/20230530/5_trajectory/monocle2"
obj.trans <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.no_dup.rds")
obj.trans <- AddMetaData(obj.trans, Idents(obj.trans), "CellType")
# input1: fib2 + fib3
# {
#   obj.trans <- obj.trans[,obj.trans$CellType %in% c('Fibroblasts-2','Fibroblasts-3')]  
# }

# input2: fib2 + fib3 + VSMC, pseudotime is not same with menstrual stages
# {
#   obj.trans <- obj.trans[,obj.trans$CellType %in% c('Fibroblasts-2','Fibroblasts-3',
#                                                     'Vascular Smooth Muscle Cells')]
# }

# input3: fib2 + fib3 + proStromal, pseudotime is not same with menstrual stages
# {
#   obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.20230906.rds")
#   cells_proFib <- intersect(
#     obj.genes@meta.data %>% filter(genes_fine == "Progenitor Fibroblasts") %>% rownames(),
#     obj.trans@meta.data %>% filter(trans_fine == "Progenitor Cells") %>% rownames()
#   )
#   cells_fib23 <- obj.trans@meta.data %>%
#     filter(trans_fine %in% c('Fibroblasts-2','Fibroblasts-3')) %>% rownames()
#   obj.trans <- obj.trans[,c(cells_proFib, cells_fib23)]
# }
# input1: fib2 + fib3 + fib1
{
  obj.trans <- obj.trans[,obj.trans$CellType %in% c('Fibroblasts-1','Fibroblasts-2','Fibroblasts-3')]  
}

# pbmc3k test data
# obj.trans <- readRDS("~/tests/scrna/pbmc3k_final.rds")
# obj.trans <- obj.trans[,obj.trans$CellType %in% c('CD14+ Mono','FCGR3A+ Mono')]
# obj.trans$condition <- obj.trans$CellType
print(obj.trans)
pdf(glue::glue("{outdir}/pseudotime_monocle2.pdf"), 12, 16)
###############################################################
##### Step0: prepare input
###############################################################
{
  data <- as(as.matrix(obj.trans@assays[["RNA"]]@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = obj.trans@meta.data)
  fData <- data.frame(gene_short_name = row.names(obj.trans@assays[["RNA"]]), 
                      row.names = row.names(obj.trans@assays[["RNA"]]))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle cds
  cds <- newCellDataSet(
    data,
    phenoData = pd,
    featureData = fd,
    # lowerDetectionLimit = 0.5,
    # expressionFamily = negbinomial.size()
  )
  
  #Filtering low-quality cells
  cds <- detectGenes(cds, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
  
  #normalize
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, cores = 20)  
  
  # filter by dispersion
  disp_table <- dispersionTable(cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
  plot_ordering_genes(cds) 
  plot_pc_variance_explained(cds, return_all = F)
  saveRDS(cds, glue::glue("{outdir}/0_0_prepare_input.rds"))
}


# unsupervised clustering
{
  cds <- reduceDimension(cds, max_components = 2, 
                         num_dim = 10, #num_dim based on PCA results
                         reduction_method = 'tSNE', 
                         # residualModelFormulaStr = "~condition + num_genes_expressed",
                         verbose = T)
  cds <- clusterCells(cds, num_clusters = 12) 
  p1 <- plot_cell_clusters(cds, 1, 2 )
  p2 <- plot_cell_clusters(cds, 1, 2, color = "CellType")
  p3 <- plot_cell_clusters(cds, 1, 2, color = "condition")
  print(p1 + p2 + p3)
  plot_cell_clusters(cds, 1, 2, color = "condition") + facet_wrap(~CellType)
  plot_cell_clusters(cds, 1, 2, color = "CellType")  + facet_wrap(~condition)
  table(pData(cds)$Cluster)
  table(pData(cds)$Cluster,pData(cds)$CellType)
  # use seurat celltype instead of cluster
  # pData(cds)$Cluster=pData(cds)$CellType
  saveRDS(cds, glue::glue("{outdir}/0_1_unspervised_cluster.rds"))
}


# semi-supervised clustering
{
  cluster_diff <- function(cds, fc=1.25, n_dim=20) {
    intext <- glue::glue("Clustering: FC = {fc}, dim = {n_dim}")
    print(intext)
    # markers <- FindAllMarkers(obj.trans, assay = "RNA", logfc.threshold = log2(fc), only.pos = TRUE) %>%
    #   group_by(cluster) %>% slice_max(abs(avg_log2FC), n = 50) %>%
    #   distinct(gene) %>% pull(gene)
    markers1 <- FindMarkers(obj.trans, assay = "RNA", logfc.threshold = log2(1.25), 
                            ident.1 = "Fibroblasts-1", only.pos = TRUE) %>% 
      slice_max(abs(avg_log2FC), n = 50) %>% rownames()
    markers2 <- FindMarkers(obj.trans, assay = "RNA", logfc.threshold = log2(1.25), 
                            ident.1 = "Fibroblasts-2", only.pos = TRUE) %>%  
      slice_max(abs(avg_log2FC), n = 50) %>% rownames()
    markers3 <- FindMarkers(obj.trans, assay = "RNA", logfc.threshold = log2(1.25), 
                            ident.1 = "Fibroblasts-3", only.pos = TRUE) %>%  
      slice_max(abs(avg_log2FC), n = 50) %>% rownames()
    markers <- c(markers1, markers2, markers3)
    print(length(markers))
    cds <- setOrderingFilter(cds, markers)
    print(plot_ordering_genes(cds))
    print(plot_pc_variance_explained(cds, return_all = F))
    
    cds <- reduceDimension(cds, max_components = 2, 
                           num_dim = n_dim, #num_dim based on PCA results
                           reduction_method = 'tSNE', 
                           # residualModelFormulaStr = "~condition+num_genes_expressed",
                           verbose = T)
    cds <- clusterCells(cds, num_clusters = 6) 
    # plot_cell_clusters(cds, 1, 2 )
    mytitle <- plot_annotation(intext)
    print(plot_cell_clusters(cds, 1, 2) / plot_cell_clusters(cds, 1, 2, color = "condition") + mytitle)
    p2 <- plot_cell_clusters(cds, 1, 2, color = "CellType")  + facet_wrap(~condition)
    p3 <- plot_cell_clusters(cds, 1, 2, color = "condition") + facet_wrap(~CellType)
    print(p2 / p3  + mytitle)
    
    {
      cds1 <- setOrderingFilter(cds, markers)
      p <- plot_ordering_genes(cds1) + ggtitle(prefix)
      print(p)
      cds1 <- reduceDimension(cds1, max_components = 2, method = 'DDRTree')
      cds1 <- orderCells(cds1)
      p1 <- plot_cell_trajectory(cds1, color_by = "State") + ggtitle(prefix) 
      p2 <- plot_cell_trajectory(cds1, color_by = "CellType") + facet_wrap(~condition, nrow = 1)
      p3 <- plot_cell_trajectory(cds1, color_by = "State") + facet_wrap(~CellType, nrow = 1)
      p4 <- plot_cell_trajectory(cds1, color_by = "Pseudotime") + ggtitle(prefix)
      print((p1+p4)/p2/p3 + mytitle)
    }
  }
  # cluster_diff(cds, fc = 1.25, n_dim = 18)
  
  pdf(glue::glue("{outdir}/0_2_clustering.pdf"), 6,8)
  args <- expand.grid(fc=c(1.25), n_dim=seq(6,21,3))
  mapply(pryr::partial(cluster_diff, cds = cds),
         args[["fc"]], args[["n_dim"]])
  dev.off()

  table(pData(cds)$Cluster,pData(cds)$CellType)
  # saveRDS(cds, glue::glue("{outdir}/0_2_semispervised_cluster.rds"))
}

saveRDS(cds, glue::glue("{outdir}/0_prepare_input.rds"))
cds <- readRDS(glue::glue("{outdir}/0_prepare_input.rds"))


###############################################################
##### Step0: prepare input
###############################################################
monocle_pseudotime <- function(cds, genes, prefix = "hvg_seurat"){
  # plot ordering genes
  print(glue::glue("Input {prefix} Genes: {length(genes)}"))
  cds1 <- setOrderingFilter(cds, genes)
  p <- plot_ordering_genes(cds1) + ggtitle(prefix)
  print(p)
  cg <- as.character(head(genes)) 
  p1 <- plot_genes_jitter(cds1[cg,],
                    grouping = "CellType",
                    color_by = "condition",
                    nrow= 3,
                    ncol = NULL )
  print(p1)
  cg2 <- as.character(tail(genes)) 
  p2 <- plot_genes_jitter(cds1[cg2,],
                    grouping = "CellType",
                    color_by = "condition",
                    nrow= 3,
                    ncol = NULL )
  print(p2)
  
  # step2: Dimension Reduction
  cds1 <- reduceDimension(cds1, max_components = 2, method = 'DDRTree')
  saveRDS(cds1, glue::glue("{outdir}/2_reduce_dimension.{prefix}.rds"))
  # cds1 <- readRDS(glue::glue("{outdir}/2_reduce_dimension.{prefix}.rds"))
  
  # step3: order cells
  cds1 <- orderCells(cds1)
  saveRDS(cds1, glue::glue("{outdir}/3_order_cells.{prefix}.rds"))
  p1 <- plot_cell_trajectory(cds1, color_by = "State") + ggtitle(prefix) 
  p2 <- plot_cell_trajectory(cds1, color_by = "CellType") + facet_wrap(~condition, nrow = 1)
  p3 <- plot_cell_trajectory(cds1, color_by = "State") + facet_wrap(~CellType, nrow = 1)
  p4 <- plot_cell_trajectory(cds1, color_by = "Pseudotime") + ggtitle(prefix)
  print((p1+p4)/p2/p3)

  # step4: pseudotime related genes
  diff_test <- differentialGeneTest(cds1[genes,],
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  write.table(diff_test, glue::glue("{outdir}/diff_pseudotime.{prefix}.tsv"),
              sep = "\t", row.names = T, quote=F)
  saveRDS(diff_test, glue::glue("{outdir}/4_pseudotime_related_genes.{prefix}.rds"))
  sig_genes <- row.names(subset(diff_test, qval < 0.01))
  p <- plot_pseudotime_heatmap(cds1[sig_genes,], num_clusters = 4,
                          cores = 4, show_rownames = T,
                          return_heatmap = F) 
  print(p)
  p <- plot_genes_in_pseudotime(cds1[head(sig_genes),], color_by = "CellType")
  print(p)
  
  # ## branch depended genes
  # BEAM_res <- BEAM(cds1, branch_point = 1, cores = 10)
  # saveRDS(BEAM_res, glue::glue("{outdir}/5_branch_related_genes.{prefix}.rds"))
  # 
  # BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  # BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  # BEAM_genes <- row.names(subset(BEAM_res,qval < 1e-4))
  # plot_genes_branched_heatmap(cds1[BEAM_genes,],
  #                             branch_point = 1, num_clusters = 4,
  #                             cores = 4,
  #                             use_gene_short_name = T, show_rownames = T)
  # # plot_genes_branched_pseudotime(cds_subset,
  # #                                branch_point = 1,
  # #                                color_by = "Hours",
  # #                                ncol = 1)
  return(cds1)
}


GM_state <- function(cds1){
  if (length(unique(pData(cds1)$State)) > 1){
    T0_counts <- table(pData(cds1)$State, pData(cds1)$condition)[,"Proliferative"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

#High Variable Genes
#Seurat
current_assay = obj.trans@active.assay
hvgs <- obj.trans@assays[[current_assay]]@var.features
# hvgs <- intersect(expressed_genes, hvgs)
cds1 <- monocle_pseudotime(cds, hvgs, 
                           prefix = "hvg_seurat")
saveRDS(cds1, glue::glue("{outdir}/cds_hvg_seurat.rds"))

## Monocle
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds2 <- monocle_pseudotime(cds, disp.genes, prefix = "hvg_monocle")
saveRDS(cds2, glue::glue("{outdir}/cds_hvg_monocle.rds"))

# Differential Genes
## cluster-specific
### Seurat
markers <- FindAllMarkers(obj.trans, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE) %>%
  group_by(cluster) %>% slice_max(abs(avg_log2FC), n = 50) %>% 
  # filter(gene %in% expressed_genes) %>% 
  distinct(gene) %>% pull(gene)
cds3 <- monocle_pseudotime(cds, markers, prefix = "markers_seurat")
{
  cdsN <- cds3
  cdsN <- orderCells(cdsN, root_state = GM_state(cdsN))
  plot_cell_trajectory(cdsN, color_by = "Pseudotime")
  plot_cell_trajectory(cdsN, color_by = "condition")
  # plot_cell_trajectory(cdsN, color_by = "State") + facet_wrap(~State, nrow = 1)
  
  ## branch depended genes
  BEAM_res <- BEAM(cdsN, branch_point = 1, cores = 10)
  # saveRDS(BEAM_res, glue::glue("{outdir}/5_branch_related_genes.{prefix}.rds"))
  write.table(BEAM_res, glue::glue("{outdir}/5_branch_related_genes.{prefix}.tsv"),
              sep = "\t", row.names = T, quote=F)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  BEAM_genes <- row.names(subset(BEAM_res,qval < 1e-4))
  
  plot_genes_branched_heatmap(cdsN[sample(BEAM_genes,100),],
                              branch_point = 1, num_clusters = 5,
                              cores = 4,
                              use_gene_short_name = T, show_rownames = T)
  # plot_genes_branched_pseudotime(cdsN[BEAM_genes[1:5],],
  #                                branch_point = 1,
  #                                color_by = "condition",
  #                                ncol = 1)
  
}

# plot_cell_trajectory(HSMM_1, color_by = "Pseudotime")
# plot_cell_trajectory(HSMM_1, color_by = "State") +
#   facet_wrap(~State, nrow = 1)
saveRDS(cds3, glue::glue("{outdir}/cds_celltype_seurat.rds"))

# ### Monocle: Need to perform order cell first
diff_test1 <- differentialGeneTest(cds[expressed_genes,],
                                   fullModelFormulaStr = "~CellType")
write.table(diff_test1, glue::glue("{outdir}/diff_celltype.tsv"),
            sep = "\t", row.names = T, quote=F)
diff_celltype <- diff_test1 %>% filter(qval < 0.01) %>% rownames()
cds4 <- monocle_pseudotime(cds, diff_celltype, prefix = "diff_celltype")
saveRDS(cds4, glue::glue("{outdir}/cds_diff_celltype.rds"))

# ## condition-specific: Need to perform order cell first
diff_test2 <- differentialGeneTest(cds[expressed_genes,],
                                   fullModelFormulaStr = "~condition")
write.table(diff_test2, glue::glue("{outdir}/diff_condition.tsv"),
            sep = "\t", row.names = T, quote=F)
diff_condition <- diff_test2 %>% filter(qval < 0.01) %>% rownames()
cds5 <- monocle_pseudotime(cds, diff_condition, prefix = "diff_condition")
saveRDS(cds5, glue::glue("{outdir}/cds_diff_condition.rds"))

# ## condition + cluster: : Need to perform order cell first 
diff_test3 <- differentialGeneTest(cds[expressed_genes,],
                                  fullModelFormulaStr = "~CellType+condition",
                                  reducedModelFormulaStr = "~condition")
write.table(diff_test3, glue::glue("{outdir}/diff_celltype_condition.tsv"),
            sep = "\t", row.names = T, quote=F)
diff_multiple <- diff_test3 %>% filter(qval < 0.01) %>% rownames()
cds6 <- monocle_pseudotime(cds, diff_multiple, prefix = "diff_multiple")
saveRDS(cds6, glue::glue("{outdir}/cds_diff_multiple.rds"))


## set root


# plot_cell_trajectory(HSMM_1, color_by = "Pseudotime")
# 
# plot_cell_trajectory(HSMM_1, color_by = "State") +
#   facet_wrap(~State, nrow = 1)
# saveRDS(HSMM, glue::glue("{outdir}/monocle2_final.rds"))

dev.off()




