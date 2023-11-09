library(Libra)
library(DEXSeq)
library(tidyr)
library(dplyr)
library(Seurat)
library(magrittr)
library(purrr)
library(glue)
library(patchwork)

outdir <- "~/projects/seq3/20230530/4_scrna/counts/9_couplings/dea"
samples <- c("Pro1", "Pro2", "Pro3", "R2F01", "R2F16")

message("Reading Releation between gene and transcripts ...")
transToGenes <- readr::read_csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv") %>%
  dplyr::select(transcript_id, gene_id, gene_name) %>% 
  mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name))

message("Reading Seurat Object...")
object.isoform <- readRDS("~/projects/seq3/20230530/4_scrna/counts/6_rename_cluster/tgs-sc_isoforms/RIF.cluster.celltype_fine.rds")
DefaultAssay(object.isoform) <- "RNA"

# filter genes: cell > 20 & counts > 50
message("Filtering Genes...")
expr <- object.isoform@assays$RNA@counts
expr <- expr[rowSums(expr >= 1) > 20 & rowSums(expr) > 50, ]
meta <- object.isoform@meta.data %>%
  mutate(label=if_else(grepl("Pro",orig.ident), "pro", "sec"),
         replicate = orig.ident, cell_type = celltype_fine)
matrices = Libra::to_pseudobulk(expr, meta = meta) %T>%
  saveRDS(glue("{outdir}/pseudobulk.input.quick.rds"))
# #######################################################
# # DEA analysis (stage): pro vs sec
# #######################################################
# # DE.isoform.limma_voom <- run_de(
# #   input = expr, meta = meta, n_threads = 20,
# #   de_family = "pseudobulk", de_method = "limma", de_type = "voom")
# 
# message("Differential Expresssion Analysis...")
# results <- bind_rows(
#   run_de(input = expr, meta = meta, n_threads = 20, 
#          de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT"),
#   run_de(input = expr, meta = meta, n_threads = 20, 
#          de_family = "pseudobulk", de_method = "edgeR", de_type = "QLF"),
#   run_de(input = round(expr), meta = meta, n_threads = 20,
#          de_family = "pseudobulk", de_method = "DESeq2", de_type = "LRT"),
#   run_de(input = expr, meta = meta, n_threads = 20,
#          de_family = "pseudobulk", de_method = "limma", de_type = "trend")
# ) %>% rename(transcript_id = gene) %>% left_join(transToGenes, by = "transcript_id") %T>%
#   readr::write_csv(glue("{outdir}/dea_all.info.quick.csv"))
# 
# 
# results %>% filter(p_val_adj < 0.05, abs(avg_logFC) > log(2)) %>% 
#   group_by(cell_type, de_method, de_type) %>% 
#   summarise(nTrans=n(), nGenes = length(unique(.data[["gene_id"]])) )  %T>%
#   readr::write_csv(glue("{outdir}/dea_all.stat.quick.csv"))
# 
# 
# # Idents(object.isoform) <- "celltype_fine"
# # DotPlot(object.isoform, assay = "RNA", 
# #         features = unique(features), cols = "Reds", 
# #         group.by = "orig.ident", idents = "eStr1") + RotatedAxis()
# # VlnPlot(object.isoform, features = features, 
# #         group.by = "orig.ident", idents = "eStr1")
# # FeaturePlot(object.isoform, features = features)
# 
# #######################################################
# # DEA cluster-to-cluster
# #######################################################
# DEAnalysis <- function(input, meta, group1, group2){
#   # group1 <- "eStr1"
#   # group2 <- "eStr2"
#   comp_groups <- paste0(group1, "_", group2)
#   message(comp_groups)
#   mat1 <- matrices[[group1]]
#   mat2 <- matrices[[group2]]
#   meta_1 <- meta %>% filter(cell_type %in% c(group1, group2)) %>%
#     mutate(label = cell_type, cell_type = comp_groups)
#   input <- input[,rownames(meta_1)]
#   results <- bind_rows(
#     run_de(input = input, meta = meta_1, n_threads = 20,
#            de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT"),
#     run_de(input = input, meta = meta_1, n_threads = 20,
#            de_family = "pseudobulk", de_method = "edgeR", de_type = "QLF"),
#     run_de(input = round(input), meta = meta_1, n_threads = 20,
#            de_family = "pseudobulk", de_method = "DESeq2", de_type = "LRT"),
#     run_de(input = input, meta = meta_1, n_threads = 20,
#            de_family = "pseudobulk", de_method = "limma", de_type = "trend")
#   ) %>% rename(transcript_id = gene) %>% left_join(transToGenes, by = "transcript_id") %T>%
#     readr::write_csv(glue("{outdir}/{comp_groups}.dea.info.quick.csv"))
#   
#   results %>% filter(p_val_adj < 0.05, abs(avg_logFC) > log(2)) %>%
#     group_by(cell_type, de_method, de_type) %>%
#     summarise(nTrans=n(), nGenes = length(unique(.data[["gene_id"]])) )  %T>%
#     readr::write_csv(glue("{outdir}/{comp_groups}.dea.stat.quick.csv"))
# }
# 
# # compare between groups
# pmap(
#   list(group1 = c("eStr1", "eStr1", "eStr2", "gEpi"),
#        group2 = c("eStr2", "eStr3", "eStr3", "lEpi")),
#   ~ DEAnalysis(expr, meta, ..1, ..2)
# )

#######################################################
# DEXSeq analysis (stage): pro vs sec
#######################################################
DEXSeqAnalysis <- function(counts, group, feature, prefix, sample_data = NULL){
  if(is.null(sample_data)){
    sample_data <- data.frame(
      condition = c("pro","pro","pro","sec","sec"), 
      disease = c("fib", "normal", "normal", "rif", "normal")
    )
    rownames(sample_data) <-  c("Pro1", "Pro2", "Pro3", "R2F01", "R2F16")
  }
  
  # formulaFullModel    =  ~ sample + exon + disease:exon + condition:exon
  # formulaReducedModel =  ~ sample + exon + disease:exon
  design <- ~ sample + exon + condition:exon
  dxd <- DEXSeqDataSet( counts, sample_data, design,
                        feature, group)
  dxd = estimateSizeFactors( dxd )
  dxd = estimateDispersions( dxd, BPPARAM = SnowParam(workers=10))
  plotDispEsts( dxd, main=prefix )
  dxd = testForDEU( dxd)
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", 
                                 BPPARAM = SnowParam(workers=10))
  dxr = DEXSeqResults( dxd, independentFiltering = T)
  table( dxr$padj < 0.1 )
  pgq <- perGeneQValue(dxr, p = "pvalue")
  results <- data.frame(row.names = names(pgq), q.value = pgq) %>%
    mutate(DIU = if_else(q.value < 0.1, 1, 0)) %>%
    tibble::rownames_to_column("gene_id")
  saveRDS(dxr, glue("{outdir}/{prefix}.rds"))
  readr::write_tsv(as.data.frame(dxr) %>% dplyr::filter(padj < 0.1) %>% dplyr::select(-genomicData), 
                   glue("{outdir}/{prefix}.tsv"))
  readr::write_tsv(results, glue("{outdir}/{prefix}.per_genes.tsv"))
  
  # MA
  padj = 0.1
  ylim = c(-2, 2)
  x <- rowMeans(counts(dxr, normalized = TRUE))
  y <- dxr[, grep("log2fold", colnames(dxr))[1]]
  df <- data.frame(x, y, dxr$padj < padj)
  plotMA(df, ylim, main=prefix)
  
  # PCA
  se <- SummarizedExperiment(log2(counts(dxd, normalized=TRUE) + 1),
                             colData=colData(dxd))
  colnames(se) <- seq(dim(colData(se))[1])
  p1 <- plotPCA( DESeqTransform(se), intgroup = c("exon"))
  #p2 <- plotPCA( DESeqTransform(se), intgroup = c("disease"))
  p3 <- plotPCA( DESeqTransform(se), intgroup = c("condition"))
  p4 <- plotPCA( DESeqTransform(se), intgroup = c("sample"))
  par(mar = c(0, 0, 0, 0), bg = NA)
  (p1 / p3 / p4) + plot_annotation(tag_levels = 'A', title=prefix) 
  
  # # specific gene
  # gene = "ENSG00000065665.21"
  # DEXSeq::plotDEXSeq( dxr, gene, legend = T,
  #                     expression=FALSE, splicing = TRUE,
  #                     cex=1.3, lwd=2)
}

matrices = Libra::to_pseudobulk(expr, meta = meta)
relation <- transToGenes %>% dplyr::select(transcript_id, gene_name) %>% 
  tibble::column_to_rownames(var = "transcript_id")
for (celltype in names(matrices)) {
  # celltype = "cEpi"
  message("DEXSeq analysis: ", celltype)
  counts <- matrices[[celltype]] 
  names(counts) <- stringr::str_split_i(names(counts), ":", 1)
  counts %<>% filter(if_any(all_of(samples), ~ . > 10) | rowSums(.[samples]) > 20)
  trans <- row.names(counts)
  genes <- relation[trans,]
  DEXSeqAnalysis(round(counts), genes, trans, glue("{celltype}.transcript_diu.stage.bambu"))
}

#######################################################
# DEXSeq analysis (stage): cluster vs cluster
#######################################################
cmpDEXSeqAnalysis <- function(matrices, samples, group1, group2){
  # group1 <- "eStr1"
  # group2 <- "eStr2"
  counts <- cbind(matrices[[group1]], matrices[[group2]])
  names(counts) <- paste0(stringr::str_split_i(names(counts), ":", 1), ":", 
                          rep(c(group1,group2), each = length(samples)))
  counts %<>% filter(if_any(all_of(names(counts)), ~ . > 10) | rowSums(.[names(counts)]) > 20)
  trans <- row.names(counts)
  genes <- relation[trans,]
  
  sample_data <- data.frame(
    condition =  rep(c(group1,group2), each = length(samples)), 
    disease = rep(c("fib", "normal", "normal", "rif", "normal"), times = 2)
  )
  rownames(sample_data) <- names(counts)
  prefix = glue("{group1}_{group2}.transcript_diu.cluster.bambu")
  DEXSeqAnalysis(round(counts), genes, trans, prefix, sample_data = sample_data) 
}

pmap(
  list(group1 = c("eStr1", "eStr1", "eStr2", "gEpi"),
       group2 = c("eStr2", "eStr3", "eStr3", "lEpi")),
  ~ cmpDEXSeqAnalysis(matrices, samples, ..1, ..2)
)