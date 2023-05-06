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
info <- read.table(file.path(indir, "Data/transcript_total_info.tsv"), sep='\t', header=T)
head(info)
# DEA
deg <- read.table(file.path(indir, "DEA/result_gene.tsv"), sep='\t', header=T)  %>%
    mutate(DEGpvalue=1-prob) %>% mutate(DEG=if_else(prob>0.95, 1, 0)) %>%
    rename(DEGlog2FC=log2FC) %>% select(gene, DEG, DEGlog2FC, DEGpvalue)
det <- read.table(file.path(indir, "DEA/result_trans.tsv"), sep='\t', header=T)  %>%
    mutate(DETpvalue=1-prob)  %>% mutate(DET=if_else(prob>0.95, 1, 0)) %>%
    rename(DETlog2FC=log2FC) %>% select(transcript, DET, DETlog2FC, DETpvalue)
dep <- read.table(file.path(indir, "DEA/result_protein.tsv"), sep='\t', header=T)  %>%
    mutate(DEPpvalue=1-prob)  %>% mutate(DEP=if_else(prob>0.95, 1, 0)) %>%
    rename(DEPlog2FC=log2FC) %>% select(protein, DEP, DEPlog2FC, DEPpvalue)
dea <- info %>% left_join(deg, by='gene') %>% left_join(det, by='transcript') %>% left_join(dep, by='protein')
head(dea)

# DIU
diu <- read.table(file.path(indir, "DIU/result_trans.tsv"), sep='\t', header=T)  %>%
    mutate(DIUswitching=if_else(podiumChange,1,0)) %>%
    mutate(DIU=if_else(qValue<0.05,1,0)) %>%
    mutate(DIUqvalue=round(qValue,3)) %>%
    rename(DIUchange=totalChange)  %>%
    select(gene, DIU, DIUswitching, DIUchange, DIUqvalue)
head(diu)

# DPA
dpa <- read.table(file.path(indir, "DPA/result.tsv"), sep='\t', header=T)  %>%
    mutate(DPAswitching=if_else(podiumChange,1,0)) %>%
    mutate(DPA=if_else(qValue<0.05,1,0)) %>%
    mutate(DPAqvalue=round(qValue,3)) %>%
    rename(DPAchange=totalChange)  %>%
    select(gene, DPA, DPAswitching, DPAchange, DPAqvalue)
head(dpa)

# Summary
total <- dea %>% left_join(diu, by='gene')  %>% left_join(dpa, by='gene')
write.table(total, file.path(outdir, "total_result.tsv"), sep='\t', quote=F, row.names=F)
## UpSetR
gene_stat <- total %>% select(gene, DEGlog2FC, DEGpvalue,DIUchange, DIUqvalue,DPAchange, DPAqvalue,
                              DEG, DIU, DIUswitching, DPA, DPAswitching) %>%
    distinct(gene, .keep_all = TRUE) %>% replace(is.na(.), 0)
head(gene_stat)
png(file.path(outdir, "gene.stat.png"))
upset(
    as.data.frame(gene_stat), sets = c("DEG", "DIU", "DIUswitching", "DPA", "DPAswitching"),
    queries = list(
        list(query = intersects, params = list("DIU","DIUswitching"), active = T),
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
dev.off()
trans_stat <- total %>% select(transcript,
                              DEGlog2FC, DEGpvalue,DETlog2FC, DETpvalue,
                              DEPlog2FC, DEPpvalue,DIUchange, DIUqvalue,
                              DEG, DET, DEP, DIU, DIUswitching,DPA, DPAswitching,
                              DPAchange, DPAqvalue) %>%
    replace(is.na(.), 0)
head(trans_stat)
png(file.path(outdir, "trans.stat.png"))
upset(
    trans_stat, sets = c("DEG", "DET", "DEP", "DIU", "DIUswitching", "DPA", "DPAswitching"),
    queries = list(
        list(query = intersects, params = list("DIU","DIUswitching"), active = T),
        list(query = intersects, params = list("DEG","DIU","DIUswitching"), active = T)
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
dev.off()
gene_num <- dim(select(total,gene)  %>% distinct(gene) )[1]
stat <- total %>% select(gene, DEG, DET, DEP, DIU, DIUswitching, DPA, DPAswitching) %>%
    select(DEG, DET, DEP, DIU, DIUswitching, DPA, DPAswitching) %>%
    filter(DEG != "") %>% group_by(DEG, DET, DEP, DIU, DIUswitching,DPA, DPAswitching) %>%
    summarise(count=n())  %>% mutate(percentage=round(count/gene_num, 3) )
write.table(stat, file.path(outdir, "total_result.stat.tsv"), sep='\t', quote=F, row.names=F)
