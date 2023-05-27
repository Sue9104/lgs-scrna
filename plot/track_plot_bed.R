library(argparser)
library(glue)
library(scriptName)
currDir <- normalizePath(dirname(current_filename())) 
p <- arg_parser("IGV Track Plot like IGV")
p <- add_argument(p, "count", help="input transcript counts file", type="character")
p <- add_argument(p, "gff", help="input gff file", type="character")
p <- add_argument(p, "region", help="format as chr1:1-1000", type="character")
p <- add_argument(p, "outdir", help="outdir", type="character")
p <- add_argument(p, "prefix", help="output prefix", type="character")
p <- add_argument(p, "--bin", help="dirname of executable gffread",
                  default=normalizePath("~/miniconda3/envs/trackplot/bin/"))
p <- add_argument(p, "--strand", help="mapping strand",
                  default="+")
argv <- parse_args(p)

prefix <- argv$prefix
outdir <- argv$outdir
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
outdir <- normalizePath(outdir)
# find executable path
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH=glue("{old_path}:{argv$bin}"))
print(Sys.getenv("PATH"))
# read region
posInfo <- unlist(strsplit(argv$region, ':|-'))
chrN <- posInfo[1]
chrStart <- as.numeric(posInfo[2])
chrEnd <- as.numeric(posInfo[3])
chrStrand <- ifelse(argv$strand=="+",1,-1)

library(Sushi)
library(GenomicFeatures)
#library(rtracklayer)
library(tidyverse)

# prepare gff
#outGff <- "/public/home/msu/softwares/flair/test/test_output/flair.gff3"
outGff <- glue("{outdir}/{prefix}.gff3")
cmd3 <- glue("gffread --keep-genes {argv$gff} -r {argv$region} -o {outGff}")
print(cmd3)
system(cmd3)

#countsFile <- "~/softwares/flair/test/test_output/test.collapse.isoform.read.map.counts.tsv"
transInfo <- makeTxDbFromGFF(outGff) 
transCounts <- read.table(argv$count, sep="\t", header = F, row.names = 1)
trans <- transcripts(transInfo, columns=c("tx_name"), use.names=FALSE)@elementMetadata$tx_name

input1.df = c()
input2.df = c()
i = 1
for (tran in rev(trans)){
  print(tran)
  #tran <- "HISEQ:1287:HKCG7BCX3:1:1106:11025:52034"
  nTrans = transCounts[tran,]
  info <- as.data.frame(exons(transInfo, filter=list(tx_name = tran))@ranges) %>% 
    dplyr::rename(stop=end) %>% mutate(gene = tran, type = "exon", chrom = chrN, strand = chrStrand, score=nTrans) %>% 
    dplyr::select(chrom, start, stop, gene, score, strand, type)
  nExons <- dim(info)[1]
  #info[1,"type"] = "utr"
  #info[nExons,"type"] = "utr"
  input1.df <- rbind(input1.df, info)
  infoTrans <- info %>% dplyr::slice(rep(1:n(), nTrans))
  infoTrans["gene"] = rep(paste0(tran, "_", seq(nTrans)), each = nExons)
  input2.df <- rbind(input2.df, infoTrans)
}
# plot transcript structure
pdf(glue("{outdir}/{prefix}.pdf"), height = 12, width = 9)
layout(matrix(c(1,2,3),3,1,byrow = TRUE), heights = c(2,10,1))
par(mar=c(0,1,1,1))
pg = plotGenes(input1.df,chrN,chrStart,chrEnd,
               types = input1.df$type, 
               colorby = input1.df$score, bentline = FALSE,
               labeltext=TRUE,
               plotgenetype="box")
par(mar=c(0,1,0,1))
pg = plotGenes(input2.df,chrN,chrStart,chrEnd,
               types = input2.df$type, col = "grey",
               labeltext=FALSE, plotgenetype="box")
par(mar=c(1,1,0.5,1))
labelgenome(chrN,chrStart,chrEnd,n=4,scale="Mb")
dev.off()
