#######################################################
# Couplings with Cell type
#######################################################
library(tidyr)
library(dplyr)
library(glue)
library(purrr)
library(tibble)
celltypes <- c("eStr1","eStr2","eStr3","VSMC","pro", 
               "lEpi","gEpi","cEpi","Endo",
               "Mac", "pDC", "nCD8+T", "NK")
indir <- "~/projects/seq3/20230530/4_scrna/counts/9_couplings/"
object.isoform <- readRDS("~/projects/seq3/20230530/4_scrna/counts/6_rename_cluster/tgs-sc_isoforms/RIF.cluster.celltype_fine.rds")
cells <- obj.isoform@meta.data %>% group_by(celltype_fine) %>% summarise(n=n()) %>% tibble::column_to_rownames(var = "celltype_fine")
cell_ratios <- mutate(cells, ratio = round(n / sum(n), 3))
cell_ratios


count_domP <- function(celltype) {
  genes.domP <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% 
    distinct(gene_id) %>% .$gene_id
  length(genes.domP) 
}
count_domPAS <- function(celltype) {
  genes.domPAS <- readr::read_table(file.path(indir, celltype, "total/pas_dominance.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% 
    distinct(gene_id) %>% .$gene_id
  length(genes.domPAS) 
}
count_LinkASTSS <- function(celltype) {
  genes.LinkASTSS <- readr::read_table(file.path(indir, celltype, "diu_orig/promoter_diu.celltype.per_genes.tsv")) %>% 
    filter(DIU==1) %>% .$gene_id
  length(genes.LinkASTSS) 
}
count_LinkASPAS <- function(celltype) {
  genes.LinkASPAS <- readr::read_table(file.path(indir, celltype, "diu_orig/pas_diu.celltype.per_genes.tsv")) %>% 
    filter(DIU==1) %>% .$gene_id
  length(genes.LinkASPAS) 
}

count_totalGenes <- function(celltype) {
  genes <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
    distinct(gene) %>% .$gene
  length(genes) 
}
count_totalGenes_5 <- function(celltype) {
  genes <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
    group_by(gene) %>% summarise(n = sum(pairs_read_counts)) %>% 
    filter(n > 5) %>% 
    distinct(gene) %>% .$gene
  length(genes) 
}
count_totalGenes_novel <- function(celltype) {
  genes <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    filter(grepl("Bambu", gene_id)) %>% distinct(gene_id) %>% .$gene_id
  length(genes) 
}
count_totalGenes_20 <- function(celltype) {
  genes <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
    group_by(gene) %>% summarise(n = sum(pairs_read_counts)) %>% 
    filter(n > 20) %>% 
    distinct(gene) %>% .$gene
  length(genes) 
}
counts_1 <- data.frame(
  celltypes = celltypes,
  cells = cells[celltypes,], 
  genes = purrr::map_vec(celltypes, count_totalGenes),
  genes_novel = purrr::map_vec(celltypes, count_totalGenes_novel),
  genes_5  = purrr::map_vec(celltypes, count_totalGenes_5),
  genes_20  = purrr::map_vec(celltypes, count_totalGenes_20),
  domP = purrr::map_vec(celltypes, count_domP),
  domPAS = purrr::map_vec(celltypes, count_domPAS),
  AS_TSS_links = purrr::map_vec(celltypes, count_LinkASTSS),
  AS_PAS_links = purrr::map_vec(celltypes, count_LinkASPAS)
)
counts_1
counts_1 %>% mutate(genes_novel = genes_novel*10, domP=domP*10, domPAS=domPAS*10) %>% 
  pivot_longer(-celltypes, names_to = "type", values_to = "counts" ) %>% 
  ggplot(aes(x=celltypes, y = counts, group = type, color = type)) + geom_line() + theme_bw() + ggsci::scale_color_aaas()
  

count_gene_type <- function(celltype) {
  db <- readr::read_table("~/pipelines/lgs-scrna/coupling/RIF.tss_tes.txt")
  info <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    mutate(transcript_id = if_else(pairType=="novelPair", pairs_id, transcript_id)) %>% 
    separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
    mutate(gene_id = gene, promoter_id = paste(gene_id, promoter, sep = ":"),
           tes_id=paste(gene_id, pas, sep = ":")) %>% 
    left_join(db, by="gene_id") %>% 
    mutate(promoter_type=case_when(promoter == distal_tss ~ "distal",
                                   promoter == proximal_tss ~ "proximal",
                                   TRUE ~ 'intermediate'),
           utr_type=case_when(pas==proximal_tes ~ "proximal",
                              pas==distal_tes ~ "distal",
                              TRUE ~ "intermediate")) 
  genes <- length(info %>% distinct(gene) %>% .$gene)
  status_tss_tes <- info %>% group_by(gene_id) %>% 
    summarise(tss=length(unique(promoter)), tes=length(unique(pas))) %>% 
    mutate(category=case_when(tss==1 & tes==1 ~ "TSS-PAS",
                              tss>1 & tes==1 ~ "ATSS-PAS",
                              tss==1 & tes>1 ~ "TSS-APA",
                              tss>1 & tes>1 ~ "ATSS-APA"))  
  counts_tss_tes <- status_tss_tes %>% group_by(category) %>% summarise(n=n()) %>% 
    pivot_wider(names_from = category, values_from = n) %>% 
    mutate(celltype = celltype, total_genes = genes, .before = 1) %>% 
    mutate(across(-c(celltype, total_genes), ~ ./total_genes * 100))
  counts_tss_tes
}

purrr::map_df(celltypes, count_gene_type)

count_trans_type <- function(celltype) {
  db <- readr::read_table("~/pipelines/lgs-scrna/coupling/RIF.tss_tes.txt")
  info <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.LATER.info.txt")) %>% 
    mutate(transcript_id = if_else(pairType=="novelPair", pairs_id, transcript_id)) %>% 
    separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
    mutate(gene_id = gene, promoter_id = paste(gene_id, promoter, sep = ":"),
           tes_id=paste(gene_id, pas, sep = ":")) %>% 
    left_join(db, by="gene_id") %>% 
    mutate(promoter_type=case_when(promoter == distal_tss ~ "distal",
                                   promoter == proximal_tss ~ "proximal",
                                   TRUE ~ 'intermediate'),
           utr_type=case_when(pas==proximal_tes ~ "proximal",
                              pas==distal_tes ~ "distal",
                              TRUE ~ "intermediate")) 
  transcripts <- length(info %>% distinct(transcript_id) %>% .$transcript_id)
  status_tss_tes <- info %>% group_by(gene_id) %>% 
    summarise(tss=length(unique(promoter)), tes=length(unique(pas))) %>% 
    mutate(category=case_when(tss==1 & tes==1 ~ "TSS-PAS",
                              tss>1 & tes==1 ~ "ATSS-PAS",
                              tss==1 & tes>1 ~ "TSS-APA",
                              tss>1 & tes>1 ~ "ATSS-APA")) 
  status_singleTSS <- status_tss_tes %>% filter(category == "TSS-APA") %>% 
    left_join(info, by = "gene_id") %>% mutate(n1 = n()) %>% 
    group_by(utr_type) %>% summarise(singleTSS = unique(n1), n=n()/unique(n1)*100) %>% 
    pivot_wider(id_cols = singleTSS, names_from = utr_type, values_from = n, 
                names_prefix = "PAS_")
  
  status_singlePAS <- status_tss_tes %>% filter(category == "ATSS-PAS") %>% 
    left_join(info, by = "gene_id") %>% mutate(n1 = n()) %>% 
    group_by(promoter_type) %>% summarise(singlePAS = unique(n1), n=n()/unique(n1)*100) %>% 
    pivot_wider(id_cols = singlePAS, names_from = promoter_type, values_from = n, 
                names_prefix = "TSS_")
  
  counts_single_type <- bind_cols(status_singleTSS, status_singlePAS) %>% 
    mutate(celltype = celltype, total_transcripts = transcripts, .before = 1) 
  counts_single_type
}
purrr::map_df(celltypes, count_trans_type)

indir <- "~/projects/seq3/20230530/4_scrna/counts/9_couplings/"
celltype  <- "eStr2"
count_diu <- function(celltype) {
  file <- file.path(indir, celltype, "diu_orig/exon_diu.stage.per_genes.tsv")
  input_genes <- length(readLines(file)) - 1
  diu_genes <- length(readLines( pipe(glue::glue("awk '$3==1' {file} "))))
  positive <- length(readLines(file.path(indir, celltype, "diu_orig/exon_diu.stage.tsv"))) - 1 
  exons_diu.tib <- tidyr::tibble(exonsInputGenes=input_genes, exonsDIUGenes=diu_genes, exonsDIU = positive)
  
  file <- file.path(indir, celltype, "diu_orig/transcript_diu.stage.per_genes.tsv")
  input_genes <- length(readLines(file)) - 1
  diu_genes <- length(readLines( pipe(glue::glue("awk '$3==1' {file} "))))
  positive <- length(readLines(file.path(indir, celltype, "diu_orig/transcript_diu.stage.tsv"))) - 1
  trans_diu.tib <- tidyr::tibble(transInputGenes=input_genes, transDIUGenes=diu_genes, transDIU = positive)
  diu <- bind_cols(exons_diu.tib, trans_diu.tib) %>% mutate(celltype = celltype, .before = 1) 
  diu
}
# node5
a <- purrr::map_df(setdiff(celltypes,c("cEpi", "Mac", "pDC", "nCD8+T", "NK")), count_diu)
# node3 
b <- purrr::map_df(setdiff(celltypes,c("eStr2", "cEpi", "Mac", "pDC", "nCD8+T", "NK")), count_diu)
# node6
purrr::map_df(setdiff(celltypes,c("pDC", "cEpi")), count_diu)

count_diu_bambu <- function(celltype) {
  file <- file.path(indir, glue("dea/{celltype}.transcript_diu.stage.bambu.per_genes.tsv"))
  input_genes <- length(readLines(file)) - 1
  # diu_genes <- length(readLines( pipe(glue::glue("awk '$3==1' {file} "))))
  diu_genes <- length(readr::read_delim(file) %>% filter(DIU == 1) %>% .$gene_id)
  positive <- length(readLines(file.path(indir, glue("dea/{celltype}.transcript_diu.stage.bambu.tsv")))) - 1
  ratio <- cell_ratios[celltype, "ratio"]
  trans_diu.tib <- tidyr::tibble(transInputGenes=input_genes, transDIUGenes=diu_genes, 
                                 scaledGenes_1 = round(diu_genes/ratio),
                                 scaledGenes_2 = round(diu_genes*sum(transInputGenes)/transInputGenes),
                                 transDIU = positive,
                                 scaledTrans_1 = round(positive/ratio),
                                 scaledTrans_2 = round(positive*sum(transInputGenes)/transInputGenes))
  diu <- trans_diu.tib %>% mutate(celltype = celltype, cells = cells[celltype,], .before = 1) 
  diu
}
aa <- purrr::map_df(celltypes, count_diu_bambu)
aa
aa %>% dplyr::select(celltype, cells, transInputGenes, transDIUGenes, scaledGenes_1 ) %>% 
  pivot_longer(-celltype, names_to = "type", values_to = "counts" ) %>% 
  ggplot(aes(x=celltype, y = counts, group = type, color = type)) + geom_line() + theme_bw() + ggsci::scale_color_aaas()


relations <- readr::read_csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv")
locations <- readr::read_delim("~/projects/seq3/20230530/gtf/RIF.tss_tes.txt")
info <- relations %>% left_join(locations, by = "gene_id")
count_diu_intersection <- function(celltype) {
  message(celltype)
  transInfo <- readr::read_delim(glue("{indir}/{celltype}/total/promoter_dominance.bambu.LATER.info.txt")) %>% 
    filter(! is.na(transcript_id)) %>% 
    mutate(promoter_id = stringr::str_split_i(promoter_id,":",2), tes_id = stringr::str_split_i(tes_id,":",2)) %>% 
    dplyr::select(transcript_id, promoter_type, utr_type) 
  transDomP <- readr::read_table(file.path(indir, celltype, "total/pas_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(transcript_id) %>% .$transcript_id
  transDomPAS <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(transcript_id) %>% .$transcript_id
  
  transDIU <- readr::read_delim(file.path(indir, glue("dea/{celltype}.transcript_diu.stage.bambu.tsv"))) %>% 
    rename(transcript_id = featureID) %>% dplyr::select(transcript_id) %>% 
    left_join(transInfo, by = "transcript_id") %>% 
    mutate(domP = if_else(transcript_id %in% transDomP, 1, 0), domPAS = if_else(transcript_id %in% transDomPAS, 1, 0))
  transDIU.counts <- length(transDIU$transcript_id)
  transDIU.PAS <- transDIU %>% filter(! is.na(utr_type)) %>% mutate(utr_type = paste0("PAS_",utr_type)) %>% 
    group_by(utr_type) %>% summarise(n=round(n()/transDIU.counts,3)) %>% 
    pivot_wider(names_from = "utr_type", values_from = n)
  transDIU.TSS <- transDIU %>% filter(! is.na(promoter_type)) %>% mutate(promoter_type = paste0("TSS_",promoter_type)) %>% 
    group_by(promoter_type) %>% summarise(n=round(n()/transDIU.counts,3)) %>% 
    pivot_wider(names_from = "promoter_type", values_from = n)
  cell_ratio = cell_ratios[celltype, "ratio"]
  transDIU.domP <- sum(transDIU$domP)
  transDIU.domPAS <- sum(transDIU$domPAS)
  bind_cols(transDIU.TSS, transDIU.PAS) %>% 
    mutate(celltype = celltype, cell_ratio = cell_ratio, transDIU = transDIU.counts, 
           domP = transDIU.domP, per_domP = round(transDIU.domP / transDIU.counts, 3),
           domPAS = transDIU.domPAS, per_domPAS = round(transDIU.domPAS / transDIU.counts, 3), .before = 1)  
}
a <- purrr::map_df(celltypes, count_diu_intersection)
p1 <- a %>% dplyr::select(celltype, cell_ratio, per_domP, per_domPAS ) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()                  
p2 <- a %>% dplyr::select(celltype, cell_ratio, PAS_distal, PAS_other, PAS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()
p3 <- a %>% dplyr::select(celltype, cell_ratio, TSS_distal, TSS_intermediate, TSS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas() 
p1 / (p2 + p3)

count_diu_clusters <- function(g1, g2){
  message(g1, g2)
  # g1 <- "eStr1"
  # g2 <- "eStr2"
  file <- file.path(indir, glue("dea/{g1}_{g2}.transcript_diu.cluster.bambu.per_genes.tsv"))
  input_genes <- length(readLines(file)) - 1
  diu_genes <- length(readLines( pipe(glue::glue("awk '$3==1' {file} "))))
  diu_genes.detail <- paste0(readLines( pipe(glue::glue("awk '$3==1' {file} | cut -f1"))), collapse = ",")
  positive <- length(readLines(file.path(indir, glue("dea/{g1}_{g2}.transcript_diu.cluster.bambu.tsv")))) - 1 
  exons_diu.tib <- tidyr::tibble(Cluster1 = g1, Cluster2 = g2, InputGenes=input_genes, 
                                 DIUGenes=diu_genes, DIUTrans = positive, detail = diu_genes.detail)
  exons_diu.tib  
}

pmap_df(list(group1 = c("eStr1", "eStr1", "eStr2", "gEpi"),
          group2 = c("eStr2", "eStr3", "eStr3", "lEpi")) ,
     ~ count_diu_clusters(..1, ..2))

#######################################################
# DEA
#######################################################
# cluster vs cluster
count_dea_clusters <- function(g1, g2){
  message(g1, g2)
  # g1 <- "eStr1"
  # g2 <- "eStr2"
  readr::read_csv(file.path(indir, glue("dea/{g1}_{g2}.dea.info.csv"))) %>% 
    filter(p_val_adj < 0.05, abs(avg_logFC) > log(1.5)) %>% 
    mutate(method=paste0(de_method,'_',de_type)) %>% 
    group_by(cell_type, method) %>% summarise(DEATrans = n(), DEAGenes = length(unique(gene_name))) %>% 
    pivot_wider(names_from = method, values_from = c(DEATrans, DEAGenes))
}
purrr::pmap_df(list(group1 = c("eStr1", "eStr1", "eStr2", "gEpi"),
             group2 = c("eStr2", "eStr3", "eStr3", "lEpi")) ,
        ~ count_dea_clusters(..1, ..2)) %>% 
  column_to_rownames("cell_type") -> cmpClusters
as_tibble(cbind(method = names(cmpClusters), t(cmpClusters)))
## genes
as_tibble(cbind(method = names(cmpClusters), t(cmpClusters))) %>% 
  filter(grepl('Genes', method)) %>% 
  mutate(method = stringr::str_replace(method, "DEAGenes_", ""))
## trans
as_tibble(cbind(method = names(cmpClusters), t(cmpClusters))) %>% 
  filter(grepl('Trans', method)) %>% 
  mutate(method = stringr::str_replace(method, "DEATrans_", ""))

# stage: pro vs sec
a <- readr::read_csv(file.path(indir, glue("dea/dea_all.info.csv"))) %>% 
  filter(p_val_adj < 0.05, abs(avg_logFC) > log(1.5)) %>% 
  dplyr::mutate(method=paste0(de_method,'_',de_type)) %>% 
  group_by(cell_type, method) %>% summarise(DEATrans = n(), DEAGenes = length(unique(gene_name))) 
## genes
a %>% dplyr::select(cell_type, method, DEAGenes) %>% 
  pivot_wider(names_from = method, values_from = DEAGenes)
## transcripts
a %>% dplyr::select(cell_type, method, DEATrans) %>% 
  pivot_wider(names_from = method, values_from = DEATrans)

b <- a %>% mutate(cell_ratio= cell_ratios[cell_type, "ratio"]) %>% 
  rowwise() %>% dplyr::mutate(across(where(is.integer), ~ round(./cell_ratio))) %>% 
  dplyr::select(-cell_ratio)
# aa: dexseq results
p1 <- aa %>% select(celltype, scaledGenes_1) %>% mutate(method = "DEXSeq") %>% 
  rename(cell_type = "celltype", DEAGenes = "scaledGenes_1") %>% 
  bind_rows(select(b, -c("DEATrans"))) %>% 
  rename(counts = "DEAGenes") %>% 
  replace(is.na(.), 0) %>% 
  ggplot(aes(x=cell_type, y = counts, group = method, color = method)) + geom_line() + theme_bw() + 
    ggsci::scale_color_aaas() + ggtitle("Scaled Differential Genes")
p2 <- aa %>% select(celltype, scaledTrans_1) %>% mutate(method = "DEXSeq") %>% 
  rename(cell_type = "celltype", DEATrans = "scaledTrans_1") %>% 
  bind_rows(select(b, -c("DEAGenes"))) %>% 
  rename(counts = "DEATrans") %>% 
  replace(is.na(.), 0) %>% 
  ggplot(aes(x=cell_type, y = counts, group = method, color = method)) + geom_line() + theme_bw() + 
  ggsci::scale_color_aaas() + ggtitle("Scaled Differential Transcripts")
p1 + p2

gene = "ENSG00000077232.19" #DNAJC10
p1 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr1/diu_orig/transcript_quantify.bambu.rds"))), 
                 type = "annotation", gene_id = gene)
p1
p2 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr2/diu_orig/transcript_quantify.bambu.rds"))), 
                 type = "annotation", gene_id = gene)
p2

p3 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr3/diu_orig/transcript_quantify.bambu.rds"))), 
                 type = "annotation", gene_id = gene)
p3

gene = "ENSG00000166598.16" #HSP90B1
p1 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr1/diu_orig/transcript_quantify.bambu.rds"))), 
                       type = "annotation", gene_id = gene)
p1
p2 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr2/diu_orig/transcript_quantify.bambu.rds"))), 
                       type = "annotation", gene_id = gene)
p2

p3 <- bambu::plotBambu(readRDS(file.path(indir, glue("eStr3/diu_orig/transcript_quantify.bambu.rds"))), 
                       type = "annotation", gene_id = gene)
p3
p1 + p2 + p3
#######################################################
# Intersection
#######################################################
relations <- readr::read_csv("~/projects/seq3/20230530/gtf/RIF.genes_transcripts.csv")
locations <- readr::read_delim("~/projects/seq3/20230530/gtf/RIF.tss_tes.txt")
info <- relations %>% left_join(locations, by = "gene_id")
count_dea_intersection <- function(celltype, method, type) {
  message(celltype)
  transInfo <- readr::read_delim(glue("{indir}/{celltype}/total/promoter_dominance.bambu.LATER.info.txt")) %>% 
    filter(! is.na(transcript_id)) %>% 
    mutate(promoter_id = stringr::str_split_i(promoter_id,":",2), tes_id = stringr::str_split_i(tes_id,":",2)) %>% 
    dplyr::select(transcript_id, promoter_type, utr_type) 
  transDomP <- readr::read_table(file.path(indir, celltype, "total/pas_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(transcript_id) %>% .$transcript_id
  transDomPAS <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(transcript_id) %>% .$transcript_id
  
  transDIU <- readr::read_delim(file.path(indir, glue("dea/dea_all.info.csv"))) %>% 
    dplyr::filter(cell_type == celltype, de_method == method, de_type == type, 
                  p_val_adj < 0.05, abs(avg_logFC) > log(1.5)) %>% 
    dplyr::select(transcript_id) %>% 
    left_join(transInfo, by = "transcript_id") %>% 
    mutate(dom = if_else(transcript_id %in% union(transDomP,transDomPAS), 1, 0), 
           domP = if_else(transcript_id %in% transDomP, 1, 0), 
           domPAS = if_else(transcript_id %in% transDomPAS, 1, 0))
  transDIU.counts <- length(transDIU$transcript_id)
  transDIU.PAS <- transDIU %>% filter(! is.na(utr_type)) %>% mutate(utr_type = paste0("PAS_",utr_type)) %>% 
    group_by(utr_type) %>% summarise(n=round(n()/transDIU.counts,3)) %>% 
    pivot_wider(names_from = "utr_type", values_from = n)
  transDIU.TSS <- transDIU %>% filter(! is.na(promoter_type)) %>% mutate(promoter_type = paste0("TSS_",promoter_type)) %>% 
    group_by(promoter_type) %>% summarise(n=round(n()/transDIU.counts,3)) %>% 
    pivot_wider(names_from = "promoter_type", values_from = n)
  cell_ratio = cell_ratios[celltype, "ratio"]
  transDIU.dom <- sum(transDIU$dom)
  transDIU.domP <- sum(transDIU$domP)
  transDIU.domPAS <- sum(transDIU$domPAS)
  bind_cols(transDIU.TSS, transDIU.PAS) %>% 
    mutate(celltype = celltype, cell_ratio = cell_ratio, transDIU = transDIU.counts, 
           dom = transDIU.dom, per_dom = round(transDIU.dom / transDIU.counts, 3),
           domP = transDIU.domP, per_domP = round(transDIU.domP / transDIU.counts, 3),
           domPAS = transDIU.domPAS, per_domPAS = round(transDIU.domPAS / transDIU.counts, 3), .before = 1)  
}
a <- purrr::map_df(celltypes, count_dea_intersection, "edgeR", "LRT")
p1 <- a %>% dplyr::select(celltype, cell_ratio, per_dom, per_domP, per_domPAS ) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()                  
p2 <- a %>% dplyr::select(celltype, cell_ratio, PAS_distal, PAS_other, PAS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()
p3 <- a %>% dplyr::select(celltype, cell_ratio, TSS_distal, TSS_intermediate, TSS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas() 
p1 / (p2 + p3)


b <- purrr::map_df(celltypes, count_dea_intersection, "edgeR", "QLF")
p1 <- b %>% dplyr::select(celltype, cell_ratio, per_dom, per_domP, per_domPAS ) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()                  
p2 <- b %>% dplyr::select(celltype, cell_ratio, PAS_distal, PAS_other, PAS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()
p3 <- b %>% dplyr::select(celltype, cell_ratio, TSS_distal, TSS_intermediate, TSS_proximal) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas() 
p1 / (p2 + p3)

count_dea_intersection_genes <- function(celltype, method, type) {
  message(celltype)
  genesDomP <- readr::read_table(file.path(indir, celltype, "total/pas_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(gene_id) %>% .$gene_id
  genesDomPAS <- readr::read_table(file.path(indir, celltype, "total/promoter_dominance.bambu.LATER.result.txt")) %>% 
    filter(tss.status == "ATSS", apa.status == "APA") %>% distinct(gene_id) %>% .$gene_id
  genesTSSAS <- readr::read_delim(glue("{indir}/{celltype}/total/tss_exon_couplings.per_junction.LASER.result.txt")) %>% 
    distinct(gene_id) %>% .$gene_id
  genesTESAS <- readr::read_delim(glue("{indir}/{celltype}/total/tes_exon_couplings.per_junction.LASER.result.txt")) %>% 
    distinct(gene_id) %>% .$gene_id
  genesDIU <- readr::read_delim(file.path(indir, glue("dea/dea_all.info.csv"))) %>% 
    dplyr::filter(cell_type == celltype, de_method == method, de_type == type, 
                  p_val_adj < 0.05, abs(avg_logFC) > log(1.5)) %>% 
    dplyr::distinct(gene_id) %>% 
    mutate(counts = 1,
           dom = if_else(gene_id %in% union(genesDomP,genesDomPAS), 1, 0), 
           domP = if_else(gene_id %in% genesDomP, 1, 0), 
           domPAS = if_else(gene_id %in% genesDomPAS, 1, 0),
           link = if_else(gene_id %in% union(genesTSSAS,genesTESAS), 1, 0), 
           linkTSS = if_else(gene_id %in% genesTSSAS, 1, 0), 
           linkTES = if_else(gene_id %in% union(genesTSSAS,genesTESAS), 1, 0))
  transDIU.counts <- length(transDIU$gene_id)

  cell_ratio = cell_ratios[celltype, "ratio"]
  plyr::colwise(sum)(select(genesDIU, -gene_id), na.rm = TRUE) %>% 
    mutate(celltype = celltype, cell_ratio = cell_ratio, .before = 1) %>% 
    mutate(per_dom = round(dom/counts, 3), per_domP = round(domP/counts, 3), per_domPAS = round(domPAS/counts, 3),
           per_link = round(link/counts, 3), per_linkTSS = round(linkTSS/counts, 3), per_linkTES = round(linkTES/counts, 3))
}
a <- purrr::map_df(celltypes, count_dea_intersection_genes, "edgeR", "LRT")
a %>% replace(is.na(.), 0) %>% 
  dplyr::select(celltype, cell_ratio, per_dom, per_domP, per_domPAS, per_linkTSS, per_linkTES) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()                  



b <- purrr::map_df(celltypes, count_dea_intersection_genes, "edgeR", "QLF")
b %>% replace(is.na(.), 0) %>% 
  dplyr::select(celltype, cell_ratio, per_dom, per_domP, per_domPAS, per_linkTSS, per_linkTES) %>% 
  tidyr::pivot_longer(-celltype, names_to = "type", values_to = "ratios" ) %>% 
  ggplot(aes(x=celltype, y = ratios, group = type, color = type)) + geom_line() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + ggsci::scale_color_aaas()                  
