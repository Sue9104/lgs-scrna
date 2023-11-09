library(tidyr)
library(dplyr)
obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.20230906.rds")
obj.trans <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.rds")

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


###############################################################
##### convert gene name form gencode 43 to refseq 80
###############################################################

#symbol_gencode = read.csv("/public/home/msu/genomes/hg38/gencode/gencode.v43.genes",
#                          header = F)[["V1"]]
symbol_refseq80 = read.csv("/public/home/msu/genomes/hg38/pyscenic/cisTarget/refseq_r80/mc_v10_clust/gene_based/refseq_genenames",
                          header = F)[["V1"]]
symbol_refseq80
to_refseq80 <- function(texts){
  sapply(texts, function(text) {
    #print(glue::glue("text={text}="))
    genes = unique(stringr::str_split_1(text, "\\s*,\\s*"))
    symbols = intersect(genes, symbol_refseq80)
    ifelse(length(symbols) == 0, "Unknown", symbols)    
  })
}

gene_types = read.csv("~/genomes/hg38/gencode/gencode.v43.gene_type.csv") %>% 
  distinct(gene_name, .keep_all = T) %>% 
  tibble::column_to_rownames("gene_name")
gene_types
symbol_hgcn = readr::read_delim("/public/home/msu/pipelines/lgs-scrna/regulon/hgnc_symbols.tsv",
                  delim = "\t", col_names = T)
symbol_hgcn
gene_symbols <- symbol_hgcn %>%  
  mutate(symbols = paste0(`Approved symbol`, ",",`Previous symbols`,",",`Alias symbols`, ",", `Ensembl ID(supplied by Ensembl)`)) %>% 
  mutate(symbols = stringr::str_replace_all(symbols, "^\\s*,|,\\s*$|\\s*","")) %>% 
  mutate(symbols = stringr::str_replace_all(symbols, ",,",",")) %>% 
  mutate(refseq80 = to_refseq80(symbols)) %>% 
  select(symbols, refseq80) %>% tidyr::separate_longer_delim(symbols, delim = ",") %>% 
  filter(refseq80 != "Unknown", symbols != refseq80) %>% 
  distinct(symbols, .keep_all = T) %>% 
  tibble::column_to_rownames("symbols")
  

outdir <- "~/projects/seq3/20230530/8_regulons"
# drop not useful genes: total - refseq80 - hgncTorefseq80 - protein_coding
genes_targets = rownames(obj.genes@assays$RNA@counts)
genes_drops = tibble(gene = setdiff(genes_targets, symbol_refseq80)) %>% 
  mutate(gene_type = gene_types[gene, "gene_type"]) %>% 
  filter(is.na(gene_type) | gene_type != "protein_coding") %>% .$gene
genes_used = match(setdiff(genes_targets, genes_drops), genes_targets)
gene_expr = obj.genes@assays$RNA@counts[genes_used,]
# convert gene symbol to refseq80
rownames(gene_expr) = tibble(gene = genes_targets[genes_used]) %>% 
  mutate(refseq80 = case_when(gene %in% symbol_refseq80 ~ gene,
                              gene %in% rownames(gene_symbols) ~ gene_symbols[gene, "refseq80"],
                              .default = gene)) %>% .$refseq80
gene_expr %>% as.data.frame() %>% 
  tibble::rownames_to_column("cell_id") %>% 
  readr::write_csv(glue::glue("{outdir}/expr_mat.csv"))

info <- read.csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv") %>% 
  tibble::column_to_rownames("transcript_id")
info %>% head

obj.trans2 <- obj.trans
DefaultAssay(obj.trans2) <- "RNA"
obj.trans2 <- Seurat::FindVariableFeatures(obj.trans2, nfeatures = 20000)
trans_hvg <- Seurat::VariableFeatures(obj.trans2)
#trans_hvg <- obj.trans@assays$integrated@var.features

trans_targets <- rownames(obj.trans@assays$RNA@counts)
trans_keep <- tibble(tran = trans_targets ) %>% 
  mutate(gene = info[tran, "gene_name"], gene_type = info[tran, "gene_type"]) %>% 
  filter(gene %in% genes_targets[genes_used]) %>% 
  filter(tran %in% trans_hvg) %>% 
  .$tran
trans_used = match(trans_keep, trans_targets)
trans_expr = obj.trans@assays$RNA@counts[trans_used,]
rownames(trans_expr) <- trans_targets[trans_used]
# remove duplicate cells ends with .1
trans_expr = trans_expr[, !grepl("\\.1$", colnames(trans_expr))]
trans_expr %>% as.data.frame() %>% 
  tibble::rownames_to_column("cell_id") %>% 
  readr::write_csv(glue::glue("{outdir}/trans/hvg_vs_tf_coding/trans_expr_mat.hvg_20000.csv"))

trans_all = obj.trans@assays$RNA@counts
# remove duplicate cells ends with .1
trans_all = trans_all[, !grepl("\\.1$", colnames(trans_all))]
trans_all %>% as.data.frame() %>% 
  tibble::rownames_to_column("cell_id") %>% 
  readr::write_csv(glue::glue("{outdir}/trans/hvg_vs_tf_coding/trans_expr_mat.csv"))

###############################################################
##### prepare pyscenic database
###############################################################
# TF list in trans
tfs <- read.csv("~/genomes/hg38/pyscenic/TFlist/allTFs_hg38.txt", header = F)[["V1"]]
tfs_trans <- tibble(tran = trans_targets[trans_used]) %>% 
  mutate(gene = info[tran, "gene_name"], 
         tran_type = info[tran, "transcript_type"]) %>% 
  filter(grepl("protein_coding", tran_type)) %>% 
  mutate(refseq80 = case_when(gene %in% symbol_refseq80 ~ gene,
                              gene %in% rownames(gene_symbols) ~ gene_symbols[gene, "refseq80"],
                              .default = gene)) %>% 
  filter(refseq80 %in% tfs) 
tfs_trans %>%   
  write.table("~/projects/seq3/20230530/8_regulons/trans/RIF.allTFs_hg38.trans_coding_info.txt", 
              quote=F, col.names = F, row.names = F)

tfs_trans %>%  select(tran) %>% 
  write.table("~/projects/seq3/20230530/8_regulons/RIF.allTFs_hg38.trans_coding.txt", 
              quote=F, col.names = F, row.names = F)
# ranking database
db_cols <- read.table("~/genomes/hg38/pyscenic/cisTarget/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions.tsv",
                    header = F, col.names = c("seqname", "start", "end"))
db_cols %>% mutate(name = seq(n()), score = 0, strand = ".") %>% 
  write.table("~/genomes/hg38/pyscenic/cisTarget/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")

## only include feather in TSS up+down regions and overlap > 40%
bambuGTF <- GenomicFeatures::makeTxDbFromGFF("~/projects/seq3/20230530/gtf/RIF.bambu_collapse.gtf")
GenomicFeatures::transcripts(bambuGTF) %>% as.data.frame() %>% 
  select(-c(tx_id, width)) -> info_pos
info %>% tibble::rownames_to_column("tx_name") %>% 
  left_join(info_pos, by = "tx_name") %>% 
  write.csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv", row.names = F, quote = F)
### 500bp_up_100bp_down
info_pos %>% mutate(score=0, start0 = start,
                    start = if_else(start0-500<0, 0, start0-500), 
                    end = start0 + 100) %>%  
  select(seqnames, start, end, tx_name, score, strand) %>% 
  arrange(seqnames, start) %>% 
  write.table("~/projects/seq3/20230530/8_regulons/trans/RIF.trans_500bp_up_100bp_down.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")
### 10kbp_up_10kbp_down
info_pos %>% mutate(score=0, start0 = start,
                    start = if_else(start0-5000<0, 0, start0-5000), 
                    end = start0 + 10000) %>%  
  select(seqnames, start, end, tx_name, score, strand) %>% 
  arrange(seqnames, start) %>% 
  write.table("~/projects/seq3/20230530/8_regulons/trans/RIF.trans_5kbp_up_5kbp_down.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")
## read results: bedtools intersect -a ref -b trans -f 0.4 -wao -nonamecheck
read.table("~/projects/seq3/20230530/8_regulons/trans/RIF.hg38_screen.1000bp_up_100bp_down.overlap_40.bed", header = F) %>% 
  filter(V13>0) %>% 
  group_by(V10) %>% filter(V13== max(V13)) %>% distinct(V10, .keep_all = T) %>% 
  write.table("~/projects/seq3/20230530/8_regulons/trans/databases/RIF.hg38_screen.1000bp_up_100bp_down.overlap_40.bed", 
              sep="\t", col.names  = F, row.names = F, quote = F)
## nearest distance among sites
## bedtools closest -a RIF.trans_1000bp_up_100bp_down.sorted.bed -b ~/genomes/hg38/pyscenic/cisTarget/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions.sorted.bed -D ref
data <- read.table("~/projects/seq3/20230530/8_regulons/trans/RIF.trans_screen.1000bp_up_100bp_down.closest.bed", header = F)
data %>% filter(abs(V13)<100, V9 != -1) %>%  
  rowwise() %>% 
  mutate(dist = min(abs(V2-V8), abs(V2-V9), abs(V3-V8), abs(V3-V9)),
         overlap = case_when(V9<V2 | V8 > V3 ~ 0,
                             V8 <= V2 & V9 >= V2 & V9 <= V3 ~ V9 - V2,
                             V8 <= V2 & V9 > V2  & V9 >  V3 ~ V3 - V2,
                             V8 >= V2 & V8 <= V3 & V9 <= V3 ~ V9 - V8,
                             V8 >= V2 & V8 <= V3 & V9 > V3  ~ V3 - V8)) %>% 
  group_by(V4) %>% filter(dist==min(dist)) %>% 
  group_by(V4) %>% filter(overlap==max(overlap)) %>% 
  group_by(V4) %>% arrange(desc(V8)) %>% distinct(V4, .keep_all = T) %>% 
  # select(V2,V3,V8,V9, V13, dist, overlap) %>% head(50) %>% View()
  write.table("~/projects/seq3/20230530/8_regulons/trans/databases/RIF.hg38_screen.1000bp_up_100bp_down.closest.bed", 
            sep="\t", col.names  = F, row.names = F, quote = F)
## to_trans_feather.py

# motif2TF
motif2TF <- read.table("~/genomes/hg38/pyscenic/motif2TF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", 
                       sep = "\t", comment.char = "", header = T)
motif2TF %>% distinct(gene_name) %>% lengths
refseq80_names <- info %>% tibble::rownames_to_column("tx_id") %>% distinct(gene_name) %>% 
  mutate(
    refseq80 = case_when(
      gene_name %in% symbol_refseq80 ~ gene_name,
      gene_name %in% rownames(gene_symbols) ~ gene_symbols[gene_name, "refseq80"],
      .default = gene_name)) 
info %>% tibble::rownames_to_column("tx_id") %>% 
  left_join(refseq80_names, by = "gene_name") %>% 
  select(-gene_name) %>% rename(gene_name = refseq80) %>% 
  select(tx_id, gene_name) %>% 
  right_join(motif2TF, by = "gene_name", relationship = "many-to-many") %>% 
  select(-gene_name) %>% rename(gene_name = tx_id) %>% 
  select(X.motif_id, motif_name,	motif_description,
         source_name,	source_version,	gene_name,
         motif_similarity_qvalue,	similar_motif_id,	similar_motif_description,
         orthologous_identity,	orthologous_gene_name,	orthologous_species,
         description) %>% 
  write.table("~/genomes/hg38/pyscenic/motif2TF/RIF_trans.motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
              col.names = T, row.names = F, quote = F, sep="\t")
## then change colname: sed -i 's|X.motif_id|#motif_id|' 

###############################################################
##### prepare pyscenic database
###############################################################
total_cells <- Cells(obj.trans)
nodup_cells <- setdiff(total_cells, total_cells[grep("\\.1$", total_cells)])
obj.trans2 <- subset(obj.trans, cells = nodup_cells)
saveRDS(obj.trans2, "~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.no_dup.rds")

obj.trans@meta.data %>% mutate(cell_type = fine_abbrs[trans_fine]) %>% 
  select(cell_type) %>%  
  tibble::rownames_to_column("CellID") %>% 
  readr::write_delim("~/projects/seq3/20230530/8_regulons/trans/hvg_vs_tf_coding/celltype.trans_fine.tsv",
                     delim = "\t")
metadata = obj.genes@meta.data %>% select(genes_fine)
tidyr::tibble(cellID = cells) %>% mutate(
  cell_type = ifelse(cellID %in% rownames(metadata), 
                     metadata[cellID, "genes_fine"], "Unknown")
) %>% mutate(cell_type = fine_abbrs[cell_type]) %>% 
  readr::write_delim("~/projects/seq3/20230530/8_regulons/trans/hvg_vs_tf_coding/celltype.genes_fine.tsv",
                     delim = "\t")

###############################################################
##### prepare celltype for input
###############################################################
obj.genes@meta.data %>% select(genes_main) %>% 
  rename(cell_type = genes_main) %>% 
  tibble::rownames_to_column("CellID") %>% 
  readr::write_delim("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/celltype.genes_main.tsv")
obj.genes@meta.data %>% mutate(cell_type = fine_abbrs[genes_fine]) %>% 
  select(cell_type) %>%  
  tibble::rownames_to_column("CellID") %>% 
  readr::write_delim("/public/home/msu/projects/seq3/20230530/8_regulons/genes_fine/celltype.genes_fine.tsv", 
                     delim = "\t")

cells <- obj.genes@meta.data %>% rownames()
metadata = obj.trans@meta.data %>% select(trans_fine)

tidyr::tibble(cellID = cells) %>% mutate(
  cell_type = ifelse(cellID %in% rownames(metadata), 
                     metadata[cellID, "trans_fine"], "Unknown")
) %>% mutate(cell_type = fine_abbrs[cell_type]) %>% 
  readr::write_delim("/public/home/msu/projects/seq3/20230530/8_regulons/trans_fine/celltype.trans_fine.tsv",
                     delim = "\t")


# debug: Less than 80% of the genes in Regulon for ABCF2
ims <- readr::read_delim(glue::glue("{outdir}/expr_mat.adjacencies.tsv"), delim = "\t")
ims %>% filter(TF=="ABCF2") %>% .$target -> targets
targets = rownames(obj.genes@assays$RNA@counts)
leng

un_genes = setdiff(targets, symbol_refseq80)
length(intersect( un_genes, rownames(gene_symbols)))
setdiff( un_genes, rownames(gene_symbols))


###############################################################
##### prepare celltype for input
###############################################################
library(Seurat)
library(SeuratDisk)
SaveH5Seurat(obj, filename = "pbmc3k.h5Seurat")
Convert("pbmc3k.h5Seurat", dest = "h5ad")


