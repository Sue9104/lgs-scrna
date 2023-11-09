library(bambu)
library(tidyverse)
library(glue)
samples <- c("Pro1", "Pro2", "Pro3", "R2F01", "R2F16")
bamdir <- "/public/home/msu/projects/seq3/20230530/3_collapse/3_1_split_bam_by_CB/"
bams <- Sys.glob(file.path(bamdir, "*.bam"))
print(bams)

outdir <- "/public/home/msu/projects/seq3/20230530/3_collapse/3_2_bambu_by_CB/default"
fa.file <- "/public/home/msu/genomes/hg38/ensembl/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa"
gtf.file <- "/public/home/msu/genomes/hg38/gencode/gencode.v43.annotation.gtf"

bambuAnnotations <- prepareAnnotations(gtf.file)
se.p <- bambu(reads = bams, annotations = bambuAnnotations, genome = fa.file, discovery = FALSE, quant = FALSE, ncore = 16)
saveRDS(se.p, glue("{outdir}/RIF.bamsboo.rcFiles.rds"))

se <- bambu(reads = se.p, annotations = bambuAnnotations, genome = fa.file, trackReads = T, ncore = 16,
            opt.discovery = list(remove.subsetTx = F, min.sampleNumber = 5))
saveRDS(se, glue("{outdir}/RIF.bamsboo.default.rds"))
writeBambuOutput(se, path = outdir)
bamTrackReads <- metadata(se)$readToTranscriptMaps[[1]] %>% rowwise() %>% mutate(equalMatches=paste(equalMatches, collapse=','), compatibleMatches=paste(compatibleMatches, collapse=',')) %>% ungroup()
write.table(bamTrackReads, glue("{outdir}/reads_to_trans.tsv"), sep = "\t", quote = F, row.names = F)
transInfo <- as.data.frame(rowData(se)) %>% rowwise() %>% mutate(eqClassById=paste(eqClassById, collapse=','))
write.table(transInfo, glue("{outdir}/transcript_info.tsv"), sep = "\t", quote = F, row.names = F)
transExons <- as.data.frame(rowRanges(se))
write.table(transExons, glue("{outdir}/transcript_exons.tsv"), sep = "\t", quote = F, row.names = F)

#ndrs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
#for (ndr in ndrs){
#  se <- bambu(reads = se.p, annotations = bambuAnnotations, genome = fa.file, trackReads = T, ncore = 5, NDR = ndr,
#              opt.discovery = list(remove.subsetTx = F))
#  hndr <- ndr * 100
#  saveRDS(se, glue("RIF.bamsboo.NDR_{hndr}.rds"))
#  outdir.hndr <- glue("./NDR_{hndr}")
#  writeBambuOutput(se, path = outdir.hndr)
#  bamTrackReads <- metadata(se)$readToTranscriptMaps[[1]] %>% rowwise() %>% mutate(equalMatches=paste(equalMatches, collapse=','), compatibleMatches=paste(compatibleMatches, collapse=',')) %>% ungroup()
#  write.table(bamTrackReads, glue("{outdir.hndr}/reads_to_trans.tsv"), sep = "\t", quote = F, row.names = F)
#  transInfo <- as.data.frame(rowData(se)) %>% rowwise() %>% mutate(eqClassById=paste(eqClassById, collapse=','))
#  write.table(transInfo, glue("{outdir.hndr}/transcript_info.tsv"), sep = "\t", quote = F, row.names = F)
#  transExons <- as.data.frame(rowRanges(se))
#  write.table(transExons, glue("{outdir.hndr}/transcript_exons.tsv"), sep = "\t", quote = F, row.names = F)
#}

