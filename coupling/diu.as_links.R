library(argparser)
p <- arg_parser("DIU for ATSS, APA and exon")
p <- add_argument(p, "indir", help="input directory", type="character")
p <- add_argument(p, "outdir", help="output directory", type="character")
p <- add_argument(p, "--bams", help="require full length bam from later, default is {indir}/{sample}/{sample}.full_length.bam")
p <- add_argument(p, "--gtf", help="annotation gtf",
                  default="~/projects/seq3/20230530/3_collapse/3_2_bambu_by_CB/default/extended_annotations.exclude_nostrand.gtf")
p <- add_argument(p, "--exonGTF", help="annotated exons gtf",
                  default="~/projects/seq3/20230530/gtf/RIF.exons.gtf")
p <- add_argument(p, "--fa", help="genome fasta file",
                  default="~/genomes/hg38/ensembl/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa")

argv <- parse_args(p)

library(parallel)
library(DEXSeq)
library(GenomicFeatures)
library(GenomicAlignments)
library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(bambu)
library(patchwork)

gtf <- argv$gtf
exonGTF <- argv$exonGTF
fa <- argv$fa
indir <- argv$indir
outdir <- argv$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
samples <- c("Pro1", "Pro2", "Pro3", "R2F01", "R2F16")
pdf(glue::glue("{outdir}/DEXSeq.diu.pdf"))
#######################################################
# DEXSeq analysis
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
  dxd = estimateDispersions( dxd, BPPARAM = SnowParam(workers=4))
  plotDispEsts( dxd, main=prefix )
  dxd = testForDEU( dxd)
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", 
                                 BPPARAM = SnowParam(workers=4))
  dxr = DEXSeqResults( dxd, independentFiltering = T)
  table( dxr$padj < 0.1 )
  pgq <- perGeneQValue(dxr, p = "pvalue")
  results <- data.frame(row.names = names(pgq), q.value = pgq) %>%
    mutate(DIU = if_else(q.value < 0.1, 1, 0)) %>%
    tibble::rownames_to_column("gene_id")
  saveRDS(dxr, glue::glue("{outdir}/{prefix}.rds"))
  readr::write_tsv(as.data.frame(dxr) %>% dplyr::filter(padj < 0.1) %>% dplyr::select(-genomicData), 
                   glue::glue("{outdir}/{prefix}.tsv"))
  readr::write_tsv(results, glue::glue("{outdir}/{prefix}.per_genes.tsv"))
  
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
#######################################################
# DEXSeq for differential analysis
# ATSS, APA, ATSS-APA
# input from LATER
#######################################################
files <- purrr::map_chr(samples, ~glue::glue("{indir}/{.x}/promoter_dominance.LATER.info.txt"))
data <- readr::read_delim(files, id = "sample") %>%
  mutate(sample=stringr::str_split_i(dirname(sample), "/", -1),
         transcript_id = if_else(pairType=="novelPair", pairs_id, transcript_id)) %>%
  separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>%
  mutate(gene_id = gene) %>%  
  dplyr::select(sample, gene_id, transcript_id, promoter, pas, pairs_read_counts)
stat <- data %>% group_by(sample, gene_id, promoter, pas) %>%
  summarise(counts = sum(pairs_read_counts)) %>% ungroup()
data %>% pivot_wider(
  id_cols = c(gene_id, transcript_id, promoter, pas),
  names_from = sample, values_from = pairs_read_counts) %>%
  drop_na() -> data
# cond: promoter or pas
for (cond in list("promoter", "pas", c("promoter", "pas"))){
  message("DIU analysis: ", cond, "...")
  cond_data <- stat %>% group_by_at(c("sample","gene_id", all_of(cond))) %>%
    summarise(counts=sum(counts)) %>%
    pivot_wider(id_cols = c(gene_id, cond), names_from = sample, values_from = counts) %>%
    right_join(dplyr::select(data, -all_of(samples)), by = c("gene_id", cond)) %>% ungroup()
  cond <- paste(cond, collapse = "_")
  counts <- as.matrix(cbind(dplyr::select(data, all_of(samples)),
                            dplyr::select(cond_data, all_of(samples))))
  # only work for current RIF analysis
  sample_data <- data.frame(
    condition=rep(c("isoform", cond), each = length(samples)),
    stage=rep(c("pro","pro","pro","sec","sec"), times = 2)``
  )
  rownames(sample_data) <-  c(samples, paste(cond,samples, sep ="."))
  group <- data$gene_id
  if (cond != "promoter_pas") {
    cols <- unique(c(cond, "promoter", "pas"))
  } else {
    cols <- c("promoter", "pas")
  }
  feature <- data[, cols]
  feature <- apply( feature, 1, paste, collapse = "" )

  prefix <- glue::glue("{cond}_diu.celltype")
  DEXSeqAnalysis(counts, group, feature, prefix, sample_data = sample_data)
}
#######################################################
# DEXSeq for differential isoform analysis
# Exon usage
#######################################################
message("DIU analysis: stage ...")
# prepare annotation
# txdb = makeTxDbFromGFF(gtf)
# flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
# names(flattenedAnnotation) =
#   sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
# seqlevelsStyle(flattenedAnnotation) = "UCSC"

# get exon counts
f <- function(x) system(glue::glue("~/miniconda3/envs/seq3/bin/featureCounts {indir}/{x}/{x}.full_length.bam -a {exonGTF} -o {outdir}/{x}.full_length.exon_counts.tsv --extraAttributes exon_number -O -f --readExtension5 50 --readExtension3 150 -L --minOverlap 20; "))
mclapply(samples, f, mc.cores = length(samples))

files <- purrr::map_vec(samples, ~glue::glue("{outdir}/{.}.full_length.exon_counts.tsv"))
exon_counts <- readr::read_delim(files, id = "sample", delim = "\t", comment = "#",
                  col_names = c("gene_id", "seqnames", "start", "end", "strand", 
                                "length", "exon_number", "counts")) %>%
  filter(gene_id != "Geneid") %>%
  mutate(sample = stringr::str_split_i(basename(sample), "\\.", 1),
         exon = paste0("E", sprintf("%02d", as.integer(exon_number)))) %>%
  dplyr::select(sample, gene_id, exon, counts) %>%
  pivot_wider(id_cols = c(gene_id, exon), names_from = sample, values_from = counts) %>%
  mutate(across(all_of(samples), as.integer)) %>%
  filter(if_any(samples, ~ . > 10) | rowSums(.[samples]) > 20)

# build dexseq
counts <- as.matrix(exon_counts %>% dplyr::select(all_of(samples)))
DEXSeqAnalysis(counts, exon_counts$gene_id, exon_counts$exon, "exon_diu.stage")


#######################################################
# DEXSeq for differential isoform analysis
# isoform usage
#######################################################
message("Differential isoform analysis: stage ...")
# quantification using bambu
#bambuAnnotations <- prepareAnnotations(gtf)
bambuAnnotations <- readRDS("~/pipelines/lgs-scrna/coupling/RIF.bambu.annotation.rds")
bams <- purrr::map_vec(samples, ~glue::glue("{indir}/{.}/{.}.full_length.bam"))
se.quantOnly <- bambu(reads = bams, annotations = bambuAnnotations, genome = fa, 
                      discovery = FALSE, ncore = 1)
saveRDS(se.quantOnly, glue::glue("{outdir}/transcript_quantify.bambu.rds"))
# gene = "ENSG00000008952.17"
# plotBambu(se.quantOnly, type = "annotation", gene_id = gene)

# prepare input
relation <- as.data.frame(mcols(se.quantOnly)) %>% dplyr::select(TXNAME,GENEID) 
transcript_counts <- as.data.frame(assays(se.quantOnly)$counts) %>% 
  tibble::rownames_to_column("TXNAME") %>% 
  pivot_longer(-TXNAME, names_to = "sample", values_to = "counts") %>% 
  mutate(sample = stringr::str_split_i(basename(sample), "\\.", 1)) %>% 
  pivot_wider(id_cols = TXNAME, names_from = sample, values_from = counts) %>%
  left_join(relation, by = "TXNAME") %>% 
  mutate(across(all_of(samples), as.integer)) %>% 
  filter(if_any(samples, ~ . > 10) | rowSums(.[samples]) > 20)

counts <- as.matrix(transcript_counts %>% dplyr::select(all_of(samples)))
DEXSeqAnalysis(counts, transcript_counts$GENEID, transcript_counts$TXNAME, 
               "transcript_diu.stage")

dev.off()




