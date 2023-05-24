library(argparser)
library(tidyverse)
library(UpSetR)
p <- arg_parser("Generate Matrix for transcript, gene and protein")
p <- add_argument(p, "indir", help="analysis directory", type="character")
p <- add_argument(p, "outdir", help="output directory")
argv <- parse_args(p)
if (!dir.exists(argv$outdir)) {dir.create(argv$outdir, recursive = TRUE) }
indir <- argv$indir
outdir <- argv$outdir

# read relation
info <- read.table(file.path(indir, "Data/transcript_total_info.tsv"), sep='\t', header=T) %>% group_by(gene) %>% mutate(nTrans=n()) %>% ungroup
head(info)
# DEA
deg <- read.table(file.path(indir, "DEA/result_gene.tsv"), sep='\t', header=T)  %>%
    mutate(DEGpvalue=1-prob) %>% mutate(DEG=if_else((prob>0.99) & (abs(log2FC)>log2(2)), 1, 0)) %>%
    rename(DEGlog2FC=log2FC) %>% select(gene, DEG, DEGlog2FC, DEGpvalue)
det <- read.table(file.path(indir, "DEA/result_trans.tsv"), sep='\t', header=T)  %>%
    mutate(DETpvalue=1-prob)  %>% mutate(DET=if_else((prob>0.95) & (abs(log2FC)>log2(1.5)), 1, 0)) %>%
    rename(DETlog2FC=log2FC) %>% select(transcript, DET, DETlog2FC, DETpvalue)
dep <- read.table(file.path(indir, "DEA/result_protein.tsv"), sep='\t', header=T)  %>%
    mutate(DEPpvalue=1-prob)  %>% mutate(DEP=if_else((prob>0.95) & (abs(log2FC)>log2(1.5)), 1, 0)) %>%
    rename(DEPlog2FC=log2FC) %>% select(protein, DEP, DEPlog2FC, DEPpvalue)
dea <- info %>% left_join(deg, by='gene') %>% left_join(det, by='transcript') %>% left_join(dep, by='protein')
nDET <- dea %>% group_by(gene, DET) %>% summarise(nDET=n()) %>% filter(DET==1) %>% select(gene, nDET)
dea <- dea %>% left_join(nDET, by="gene")
head(dea)

# DIU
diu <- read.table(file.path(indir, "DIU/result_trans.tsv"), sep='\t', header=T)  %>%
    mutate(DIUswitching=if_else(podiumChange,1,0)) %>%
    mutate(DIU=if_else(qValue<0.1,1,0)) %>%
    mutate(DIUqvalue=round(qValue,3)) %>%
    rename(DIUchange=totalChange)  %>%
    select(gene, DIU, DIUswitching, DIUchange, DIUqvalue)
head(diu)

# DPA
dpa <- read.table(file.path(indir, "DPA/result.tsv"), sep='\t', header=T)  %>%
    mutate(DPAswitching=if_else(podiumChange,1,0)) %>%
    mutate(DPA=if_else(qValue<0.1,1,0)) %>%
    mutate(DPAqvalue=round(qValue,3)) %>%
    rename(DPAchange=totalChange)  %>%
    select(gene, DPA, DPAswitching, DPAchange, DPAqvalue)
head(dpa)

# Summary
total <- dea %>% left_join(diu, by='gene')  %>% left_join(dpa, by='gene')
write.table(total, file.path(outdir, "result.orig.tsv"), sep='\t', quote=F, row.names=F)
sink(file.path(outdir, "result.summary.txt"))
cat("# Transcripts\n")
nTrans <- dim(total)[1]
nTrans.fsm <- (total %>% filter(startsWith(category, 'full')) %>% lengths)[["transcript"]]
nTrans.ism <- (total %>% filter(startsWith(category, 'incomplete')) %>% lengths)[["transcript"]]
nTrans.nic <- (total %>% filter(startsWith(category, 'novel_in_catalog')) %>% lengths)[["transcript"]]
nTrans.nnc <- (total %>% filter(startsWith(category, 'novel_not_in_catalog')) %>% lengths)[["transcript"]]
nTrans.fusion <- (total %>% filter(startsWith(category, 'fusion')) %>% lengths)[["transcript"]]
nTrans.genic <- (total %>% filter(startsWith(category, 'genic')) %>% lengths)[["transcript"]]
cat(c("Transcripts","FSM","ISM","NIC","NNC","fusion","genic\n"), sep="\t")
cat(c(nTrans,nTrans.fsm,nTrans.ism,nTrans.nic,nTrans.nnc,nTrans.fusion,nTrans.genic,"\n"), sep="\t")
cat(c(nTrans,
      sprintf("%d(%.2f%%)", nTrans.fsm, nTrans.fsm/nTrans*100),
      sprintf("%d(%.2f%%)", nTrans.ism, nTrans.ism/nTrans*100),
      sprintf("%d(%.2f%%)", nTrans.nic, nTrans.nic/nTrans*100),
      sprintf("%d(%.2f%%)", nTrans.nnc, nTrans.nnc/nTrans*100),
      sprintf("%d(%.2f%%)", nTrans.fusion, nTrans.fusion/nTrans*100),
      sprintf("%d(%.2f%%)", nTrans.genic, nTrans.genic/nTrans*100),"\n"), sep="\t")

cat("\n# Genes\n")
nGenes <- (total %>% distinct(gene) %>% lengths)[["gene"]]
nGenes.multi <-(total %>% group_by(gene, .drop = T) %>% summarise(n=n()) %>% filter(n>1) %>% distinct(gene) %>% lengths)[["gene"]]
nGenes.novel <-(total %>% filter(startsWith(gene, "novelGene")) %>% distinct(gene) %>% lengths)[["gene"]]
cat(c("Genes","Genes:MultiTrans","Genes:Novel\n"), sep="\t")
cat(c(nGenes,nGenes.multi,nGenes.novel,"\n"), sep="\t")
cat(c(nGenes,sprintf("%d(%.2f%%)", nGenes.multi, nGenes.multi/nGenes*100),
      sprintf("%d(%.2f%%)", nGenes.novel, nGenes.novel/nGenes*100),"\n"), sep="\t")

cat("\n# Results\n")
nDEG <- (total %>% filter(DEG==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDET <- (total %>% filter(DET==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDEP <- (total %>% filter(DEP==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDIU <- (total %>% filter(DIU==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDIUswitching <- (total %>% filter(DIUswitching==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDPA <- (total %>% filter(DPA==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDPAswitching <- (total %>% filter(DPAswitching==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDEGDIU <- (total %>% filter(DEG==1, DIU==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDIUboth <- (total %>% filter(DIU==1, DIUswitching==1) %>% distinct(gene) %>% lengths)[["gene"]]
nDPAboth <- (total %>% filter(DPA==1, DPAswitching==1) %>% distinct(gene) %>% lengths)[["gene"]]
cat(c("DEG","DET","DEP","DIU","DIUswitching","DPA","DPAswitching","DEG+DIU","DIU+DIUswitching","DPA+DPAswitching\n"), sep="\t")
cat(c(nDEG,nDET,nDEP,nDIU,nDIUswitching,nDPA,nDPAswitching,nDEGDIU,nDIUboth,nDPAboth,"\n"), sep="\t")
cat(c(sprintf("%d(%.2f%%)", nDEG, nDEG/nGenes*100),
      sprintf("%d(%.2f%%)", nDET, nDET/nGenes*100),
      sprintf("%d(%.2f%%)", nDEP, nDEP/nGenes*100),
      sprintf("%d(%.2f%%)", nDIU, nDIU/nGenes*100),
      sprintf("%d(%.2f%%)", nDIUswitching, nDIUswitching/nGenes*100),
      sprintf("%d(%.2f%%)", nDPA, nDPA/nGenes*100),
      sprintf("%d(%.2f%%)", nDPAswitching, nDPAswitching/nGenes*100),
      sprintf("%d(%.2f%%)", nDEGDIU, nDEGDIU/nGenes*100),
      sprintf("%d(%.2f%%)", nDIUboth, nDIUboth/nGenes*100),
      sprintf("%d(%.2f%%)", nDPAboth, nDPAboth/nGenes*100), "\n"), sep="\t")
sink()


## UpSetR
gene_stat <- total %>% select(gene, DEGlog2FC, DEGpvalue,DIUchange, DIUqvalue,DPAchange,
                              DEG, DIU, DIUswitching, DPA, DPAswitching, DPAqvalue) %>%
    distinct(gene, .keep_all = TRUE) %>% replace(is.na(.), 0) %>%
    mutate(DPAqvalue = if_else(DPAqvalue == 1, 0.999, DPAqvalue))
head(gene_stat)
pdf(file.path(outdir, "gene.stat.pdf"))
if (nDIU > 1){
    upset(
        as.data.frame(gene_stat), sets = c("DEG", "DIU", "DIUswitching", "DPA", "DPAswitching"),
        queries = list(
            list(query = intersects, params = list("DIU","DIUswitching"), active = T),
            list(query = intersects, params = list("DIU"), active = T),
            list(query = intersects, params = list("DEG","DIU","DIUswitching"), active = T)
        ),
        attribute.plots=list(
            gridrows=100,
            plots=list(
                #list(plot=histogram, x="DEGlog2FC",queries=T),
                list(plot=scatter_plot, x="DEGlog2FC", y="DEGpvalue",queries=T),
                list(plot=scatter_plot, x="DIUchange", y="DIUqvalue",queries=T),
                list(plot=scatter_plot, x="DPAchange", y="DPAqvalue",queries=T)),
            ncols = 1
        ),
        order.by = "freq",
        #empty.intersections = "on"
    )
} else {
    upset(
        as.data.frame(gene_stat), sets = c("DEG", "DIU", "DIUswitching", "DPA", "DPAswitching"),
        order.by = "freq",
        #empty.intersections = "on"
    )
}
dev.off()

trans_stat <- total %>% select(transcript,
                              DEGlog2FC, DEGpvalue,DETlog2FC, DETpvalue,
                              DEPlog2FC, DEPpvalue,DIUchange, DIUqvalue,
                              DEG, DET, DEP, DIU, DIUswitching,DPA, DPAswitching,
                              DPAchange, DPAqvalue) %>%
    replace(is.na(.), 0) %>%
    mutate(DPAqvalue = if_else(DPAqvalue == 1, 0.999, DPAqvalue))
head(trans_stat)
pdf(file.path(outdir, "trans.stat.pdf"))
if (nDIU>1){
    upset(
        as.data.frame(trans_stat), sets = c("DEG", "DET", "DEP", "DIU", "DIUswitching", "DPA", "DPAswitching"),
        queries = list(
            list(query = intersects, params = list("DIU","DIUswitching"), active = T),
            list(query = intersects, params = list("DIU"), active = T)
            #list(query = intersects, params = list("DEG","DIU","DIUswitching"), active = T)
        ),
        attribute.plots=list(
            gridrows=100,
            plots=list(
                #list(plot=histogram, x="DIUchange",queries=T),
                #list(plot=histogram, x="DEGlog2FC",queries=T),
                list(plot=scatter_plot, x="DEGlog2FC", y="DEGpvalue",queries=T),
                list(plot=scatter_plot, x="DETlog2FC", y="DETpvalue",queries=T),
                list(plot=scatter_plot, x="DIUchange", y="DIUqvalue",queries=T),
                list(plot=scatter_plot, x="DPAchange", y="DPAqvalue",queries=T)
            ),
            ncols = 2
        ),
        order.by = "freq",
        #empty.intersections = "on"
    )
} else {
    upset(
        as.data.frame(trans_stat), sets = c("DEG", "DET", "DEP", "DIU", "DIUswitching", "DPA", "DPAswitching"),
        order.by = "freq"
    )
}

dev.off()



