library(argparser)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
p <- arg_parser("Relation between TSS ASS TES")
p <- add_argument(p, "bam", help="input bam file", type="character")
p <- add_argument(p, "outdir", help="output directory", type="character")
p <- add_argument(p, "--gtf", help="transcript gtf",
                  default="~/projects/seq3/20230530/3_collapse/3_2_bambu_by_CB/default/extended_annotations.exclude_nostrand.gtf")
p <- add_argument(p, "--junction", help="short reads junctions", type="character",
                  default="~/projects/seq3/20230530/RIF.SJ.out.tab")
p <- add_argument(p, "--pvalue", help="sigificant value",
                  default=0.05)
argv <- parse_args(p)

outdir <- argv$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
outdir <- normalizePath(outdir)
junction_path <- argv$junction
bam_path <- argv$bam
sample <- basename(outdir)
## read reference gtf
#reference_annotation <- rtracklayer::import.gff(annot_path)
## only using protein coding genes
#protein_coding_exons <- reference_annotation[reference_annotation$type == "exon" & reference_annotation$gene_type == "protein_coding"]
#protein_coding_genes <- reference_annotation[reference_annotation$type == "gene" & reference_annotation$gene_type == "protein_coding"]


# LATER: TSS vs TES
####ref_custom <- rtracklayer::import.gff(argv$gtf)
####ref_custom <-  ref_custom[ref_custom$transcript_id != "BambuTx2068"] # remove unkonwn strand
####ref_custom_exons <- ref_custom[ref_custom$type == "exon"]
####isoformData <- LATER::prepareIsoformDatabase(ref_custom_exons, tss.window = 50, tes.window = 150)
####message("Reference Slots: ", slotNames(isoformData))
isoformData <- readRDS("~/pipelines/lgs-scrna/coupling/LATER.isoformdata.rds")

# Counting 5'-3' isoforms
countData <- LATER::countLinks(argv$bam, isoformData)
readr::write_tsv(LATER::readAssignments(countData),
                 file.path(outdir, "gene_model.read_assignments.txt"))
# using qname to extract from bam
#Sys.setenv(PATH = paste("~/miniconda3/envs/seq3/bin/", Sys.getenv("PATH"), sep=":"))
cmd = glue::glue("~/miniconda3/envs/seq3/bin/samtools view -N {outdir}/gene_model.read_assignments.txt -o {outdir}/{sample}.full_length.bam {bam_path}; ")
message(cmd)
system(cmd)

# Statistical testing of 5’-3’ couplings.
dominance_result <- LATER::estimatePromoterDominance(countData, isoformData, method = "chisq")
saveRDS(dominance_result, file.path(outdir, "promoter_dominance.LATER.rds"))
readr::write_tsv(dominance_result@dominance,
                 file.path(outdir, "promoter_dominance.LATER.info.txt"))
biasGenes <- dominance_result@result %>% filter(p.adj.chisq < argv$pvalue) %>% .$gene_id
dominancePromoters <- dominance_result@dominance %>%
  filter(pairs_read_counts > 20) %>%
  filter(gene_id %in% biasGenes) %>%
  filter(promoterDominance > 0.2, endDominance > 0.6)
readr::write_tsv(dominancePromoters,
                 file.path(outdir, "promoter_dominance.LATER.result.txt"))
dominancePAS <- dominance_result@dominance %>%
  filter(pairs_read_counts > 20) %>%
  filter(gene_id %in% biasGenes) %>%
  filter(endDominance > 0.2, promoterDominance > 0.6)
readr::write_tsv(dominancePAS,
                 file.path(outdir, "pas_dominance.LATER.result.txt"))
####referenece_junctions <- LASER::create_reference_junctions(
####  junction_path, min.jcounts = 2, ref_custom, type="short")
referenece_junctions <- readRDS("~/pipelines/lgs-scrna/coupling/LASER.refjunctions.RIF.rds")

## Full-length read filtering
exonlinks.counts<- LASER::read_to_junctions(bam_path, referenece_junctions, argv$gtf)
dat <- exonlinks.counts %>% dplyr::select(read_id , new_junID, tes_id,  promoter_id, pairs_id)
readr::write_tsv(dat,
                 file.path(outdir, "gene_junction.read_assignments.txt"))

# Calculate exon-junction couplings to 5'/3'
couplings <- LASER::calculate_exon_couplings(exonlinks.counts,  referenece_junctions )
saveRDS(couplings, file.path(outdir, "exon_couplings.LASER.rds"))
# results: counts > 20 & residuals > 0.7
couplings.tss <- couplings$TSS.couplingsPerJunction %>%
  filter(observed > 20, abs(residuals) > 0.7) %>%
  mutate(gene_id = stringr::str_split_i(pairs_id,":", 3))
couplings.tes <- couplings$TES.couplingsPerJunction %>%
  filter(observed > 20, abs(residuals) > 0.7) %>%
  mutate(gene_id = stringr::str_split_i(pairs_id,":", 3))
readr::write_tsv(couplings$TSS.couplingsPerGene,
                 file.path(outdir, "tss_exon_couplings.per_gene.LASER.txt"))
readr::write_tsv(couplings$TSS.couplingsPerJunction,
                 file.path(outdir, "tss_exon_couplings.per_junction.LASER.info.txt"))
readr::write_tsv(couplings.tss,
                 file.path(outdir, "tss_exon_couplings.per_junction.LASER.result.txt"))

readr::write_tsv(couplings$TES.couplingsPerGene,
                 file.path(outdir, "tes_exon_couplings.per_gene.LASER.txt"))
readr::write_tsv(couplings$TES.couplingsPerJunction,
                 file.path(outdir, "tes_exon_couplings.per_junction.LASER.txt"))
readr::write_tsv(couplings.tes,
                 file.path(outdir, "tes_exon_couplings.per_junction.LASER.result.txt"))

# plot specific gene
#couplings$TSS.couplingsPerJunction %>%
#  filter(grepl("FBgn0266521", pairs_id)) %>%
#  mutate(promoter_id=stringr::str_split_fixed(.data$pairs_id,":" ,n = 3)[,3],
#         junction=stringr::str_split_fixed(.data$pairs_id,":" ,n = 3)[,2]) %>%
#  ggplot(., aes(x=junction, y=residuals, color=promoter_id)) +
#  geom_point(size=3) +
#  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  facet_grid(cols=vars(promoter_id)) + scale_color_manual(values=c("#FF355E","#004F98", "#679267"))+
#  geom_hline(yintercept=0)
