source("~/softwares/LATER/R/LATER-class.R")

library(glue)
library(dplyr)
library(tidyr)
celltype = "cEpi"
indir <- "~/projects/seq3/20230530/4_scrna/counts/9_couplings/"
isoformData <- readRDS("~/pipelines/lgs-scrna/coupling/LATER.isoformdata.rds")
relation <- isoformData@isoform_database %>% select(transcript_id, gene_id, pairs_id)

domAnalysis <- function(celltype, pvalue = 0.05){
  outdir <- file.path(indir, celltype, "total")
  message(outdir)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  mat <- readRDS(file.path(indir, "dea/pseudobulk.input.rds"))[[celltype]]
  names(mat) <- stringr::str_split_i(names(mat), ':', 1)
  samples <- names(mat)
  mat %>% tibble::rownames_to_column("transcript_id") %>%
    left_join(relation, by = "transcript_id") %>% filter(!is.na(pairs_id)) %>%
    group_by(pairs_id, gene_id) %>%
    summarise(across(all_of(samples), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(read_counts = rowSums(across(all_of(samples)), na.rm=TRUE)) %>%
    dplyr::select(pairs_id, read_counts, gene_id) -> counts

  dd <- edgeR::DGEList(counts = counts$read_counts)
  dge <- edgeR::calcNormFactors(dd)
  counts[["cpm"]] <- as.numeric(edgeR::cpm(dge))

  # countData <- LATER(isoformCounts =  countData@isoformCounts)
  countData <- LATER(isoformCounts = counts)
  dominance_result <- LATER::estimatePromoterDominance(countData, isoformData, method = "chisq")
  saveRDS(dominance_result, file.path(outdir, "promoter_dominance.bambu.LATER.rds"))
  readr::write_tsv(dominance_result@dominance,
                   file.path(outdir, "promoter_dominance.bambu.LATER.info.txt"))
  biasGenes <- dominance_result@result %>% filter(p.adj.chisq < pvalue) %>% .$gene_id
  dominancePromoters <- dominance_result@dominance %>%
    filter(pairs_read_counts > 20) %>%
    filter(gene_id %in% biasGenes) %>%
    filter(promoterDominance > 0.2, endDominance > 0.6)
  readr::write_tsv(dominancePromoters,
                   file.path(outdir, "promoter_dominance.bambu.LATER.result.txt"))
  dominancePAS <- dominance_result@dominance %>%
    filter(pairs_read_counts > 20) %>%
    filter(gene_id %in% biasGenes) %>%
    filter(endDominance > 0.2, promoterDominance > 0.6)
  readr::write_tsv(dominancePAS,
                   file.path(outdir, "pas_dominance.bambu.LATER.result.txt"))

}


celltypes <- c("eStr1","eStr2","eStr3","VSMC","pro",
               "lEpi","gEpi","cEpi","Endo",
               "Mac", "pDC", "nCD8+T", "NK")
purrr::map(celltypes, ~ domAnalysis(.))

