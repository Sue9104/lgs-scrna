suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))

gtf <- "~/projects/seq3/20230530/3_collapse/3_2_bambu_by_CB/default/extended_annotations.exclude_nostrand.gtf"
ref_custom <- rtracklayer::import.gff(gtf)
ref_custom_exons <- ref_custom[ref_custom$type == "exon"]

strandSort <- function(x) {
  c(
    GenomicRanges::sort(x[x@strand == "+"], decreasing = FALSE),
    GenomicRanges::sort(x[x@strand == "-"], decreasing = TRUE)
  )
}
# ref_custom_exons_2 <- ref_custom_exons %>% group_by(transcript_id) %>% 
#   mutate(max_exon = max(as.integer(exon_number))) %>%
#   filter(exon_number != 1, exon_number != max_exon) %>% ungroup()


strandSort(
  GenomicRanges::makeGRangesFromDataFrame(
    reshape::melt(GenomicRanges::reduce(
      GenomicRanges::split(ref_custom_exons %>% ungroup(), ~gene_id),
      min.gapwidth = 5)),
    keep.extra.columns = TRUE
  )) %>% group_by(value.group_name) %>%
  mutate(gene_id = value.group_name, exon_number =  sequence(n()), type = "exon") %>% 
  ungroup() %>% plyranges::select(-value.group, -value.group_name) -> exons  

rtracklayer::export.gff2(exons, "~/projects/seq3/20230530/gtf/RIF.exons.gtf")

gene = "ENSG00000000419.14"
ref_1 <- ref_custom_exons %>% filter(gene_id == gene) 
strandSort(
  GenomicRanges::makeGRangesFromDataFrame(
    reshape::melt(GenomicRanges::reduce(
      GenomicRanges::split(ref_1 %>% ungroup(), ~gene_id),
      min.gapwidth = 5)),
    keep.extra.columns = TRUE
  )) %>% group_by(value.group_name) %>%
  mutate(gene_id = value.group_name, exon_number =  sequence(n()), type = "exon") %>% 
  ungroup() %>% plyranges::select(-value.group, -value.group_name) -> a

read_bam("~/projects/seq3/20230530/4_scrna/counts/9_couplings/test.bam", 
         index = "~/projects/seq3/20230530/4_scrna/counts/9_couplings/test.bam.bai")


