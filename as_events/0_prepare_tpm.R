library(dplyr)
library(tidyr)
library(Seurat)


obj.trans <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.no_dup.rds")
trans_all <- rownames(obj.trans@assays$RNA@counts)
counts <- obj.trans@assays$RNA@counts
trans_qc <- rownames(counts)[rowSums(counts) > 20]

info <- read.csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv") 
get_trans_in_same_gene <- function(trans_list){
  geneids <- info %>% filter(tx_id %in% trans_list) %>% .$gene_id
  info %>% filter(gene_id %in% geneids) %>% filter(tx_id %in% trans_all) %>% .$tx_id
}

obj.trans2 <- Seurat::FindVariableFeatures(obj.trans, nfeatures = 2000, assay = "RNA")
hvt <- obj.trans2@assays$RNA@var.features

fine_abbrs <- c("genitor", "pEpi", "pFib", "Fib0", "cEpi", "gEpi","lEpi","Fib1", "Fib2","Fib3", 
                "NK", "T","Mac", "DC", "Endo","VSMC", "?")
names(fine_abbrs) <- c("Progenitor Cells", 
                       "Progenitor Epithelium", "Progenitor Fibroblasts", 
                       "Fibroblasts", "Ciliated Epithelium", 
                       "Glandular Epithelium", "Lumenal Epithelium",
                       "Fibroblasts-1", "Fibroblasts-2","Fibroblasts-3", 
                       "Nature Killer Cells", "T Cells",
                       "Macrophages", "Dendritic Cells", 
                       "Endothelium", "Vascular Smooth Muscle Cells",
                       "Unknown")
fine_abbrs


obj.trans.fib <- subset(obj.trans, 
                        idents = c("Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-1"))
counts.fib <- obj.trans.fib@assays$RNA@counts
str_trans <- rownames(counts.fib)[rowSums(counts.fib) > 10]
str_hvt <- intersect(obj.trans.fib@assays$integrated@var.features, str_trans)

# transform seurat data to TPM
outdir = "~/projects/seq3/20230530/6_as_events/"
## celltype tpm
tpms <- AverageExpression(obj.trans, assays = "RNA", slot = "data")$RNA %>% 
  as.data.frame() %>% mutate(across(everything(), ~exp(.x) - 1)) 
colnames(tpms) <- fine_abbrs[colnames(tpms)]
write.table(tpms, glue::glue("{outdir}/tpm/RIF.celltype_tpm.all.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(str_trans),], glue::glue("{outdir}/tpm/RIF.celltype_tpm.qc.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(hvt),], glue::glue("{outdir}/tpm/RIF.celltype_tpm.hvt.tsv"), quote = F, sep = "\t")

trans.eStr2 <- FindMarkers(obj.trans, assay = "RNA", 
                          ident.1 = "Fibroblasts-2", ident.2 = "Fibroblasts-1", 
                          min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.5), test.use = "wilcox") %>% rownames()
trans_up.eStr2 <- FindMarkers(obj.trans, assay = "RNA", 
                           ident.1 = "Fibroblasts-2", ident.2 = "Fibroblasts-1", 
                           min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.5), test.use = "wilcox") %>% 
  filter(avg_log2FC>0) %>% rownames()
trans_down.eStr2 <- FindMarkers(obj.trans, assay = "RNA", 
                           ident.1 = "Fibroblasts-2", ident.2 = "Fibroblasts-1", 
                           min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.5), test.use = "wilcox") %>% 
  filter(avg_log2FC<0) %>% rownames()
trans.eStr3 <- FindMarkers(obj.trans, assay = "RNA", 
                          ident.1 = "Fibroblasts-3", ident.2 = "Fibroblasts-1", 
                          min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.2), test.use = "wilcox") %>% rownames()
trans_up.eStr3 <- FindMarkers(obj.trans, assay = "RNA", 
                           ident.1 = "Fibroblasts-3", ident.2 = "Fibroblasts-1", 
                           min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.2), test.use = "wilcox") %>% 
  filter(avg_log2FC>0) %>% rownames()
trans_down.eStr3 <- FindMarkers(obj.trans, assay = "RNA", 
                           ident.1 = "Fibroblasts-3", ident.2 = "Fibroblasts-1", 
                           min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = log(1.2), test.use = "wilcox") %>% 
  filter(avg_log2FC<0) %>% rownames()

write.table(tpms[get_trans_in_same_gene(c(trans.eStr2, trans.eStr3)),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(trans.eStr2),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff2.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(trans_up.eStr2),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff2_up.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(trans_down.eStr2),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff2_down.tsv"), quote = F, sep = "\t")

write.table(tpms[get_trans_in_same_gene(trans.eStr3),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff3.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(trans_up.eStr3),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff3_up.tsv"), quote = F, sep = "\t")
write.table(tpms[get_trans_in_same_gene(trans_down.eStr3),], 
            glue::glue("{outdir}/tpm/RIF.celltype_tpm.diff3_down.tsv"), quote = F, sep = "\t")


### split expression by celltype
for (ct in levels(obj.trans)){
  # ct = "Fibroblasts-2"
  ct_abbr = fine_abbrs[ct]
  AverageExpression(subset(obj.trans, idents = ct), assays = "RNA", slot = "data", group.by = "orig.ident")$RNA %>% 
    as.data.frame() %>% mutate(across(everything(), ~exp(.x) - 1)) %>% 
    write.table(glue::glue("{outdir}/tpm/RIF.sample_tpm.{ct_abbr}.tsv"), quote = F, sep = "\t")
}
sapply(colnames(tpms), 
       function(x) write.table(tpms[x], glue::glue("{outdir}/tpm/RIF.celltype_tpm.{x}.tsv"), quote = F, sep = "\t"))

## sample tpm
AverageExpression(obj.trans, assays = "RNA", slot = "data", group.by = "orig.ident")$RNA %>% head
  as.data.frame() %>% 
  mutate(across(everything(), ~exp(.x) - 1)) %>% 
  write.table(glue::glue("{outdir}/tpm/RIF.sample_tpm.tsv"), quote = F, sep = "\t")



# split psi by celltype
files <- Sys.glob(glue::glue("{outdir}/pool_genes/RIF.*.psi"))
for (file in files){
  print(file)
  base <- stringr::str_remove(basename(file), "_all.psi")
  data <- read.table(file, sep = "\t") 
  ### split expression by celltype
  sapply(colnames(data), 
         function(x) write.table(tpms[x], glue::glue("{outdir}/pool_genes/{base}_{x}.psi"), quote = F, sep = "\t"))
}



