library(argparser)
library(glue)
library(scriptName)
currDir <- normalizePath(dirname(current_filename())) 
p <- arg_parser("IGV Track Plot like IGV")
p <- add_argument(p, "bam", help="input bam file", type="character")
p <- add_argument(p, "gff", help="input gff file", type="character")
p <- add_argument(p, "region", help="format as chr1:1-1000", type="character")
p <- add_argument(p, "outdir", help="outdir", type="character")
p <- add_argument(p, "prefix", help="output prefix", type="character")
p <- add_argument(p, "--bin", help="dirname of executable gffread and samtools",
                  default=normalizePath("~/miniconda3/envs/trackplot/bin/"))
p <- add_argument(p, "--mane", help="MANE gff file",
                  default=glue("{currDir}/data/MANE.GRCh38.v1.0.ensembl_genomic.gff"))
p <- add_argument(p, "--genome", help="genome version",
                  default="hg38")
p <- add_argument(p, "--strand", help="mapping strand",
                  default="+")
argv <- parse_args(p)
print(argv$mane)

library(Gviz)
library(tidyverse)


outdir <- argv$outdir
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
outdir <- normalizePath(outdir)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH=glue("{old_path}:{argv$bin}"))
print(Sys.getenv("PATH"))

if (system2("command", args = c("-v", "samtools")) == 1) stop("samtools NOT FOUND!!!")
if (system2("command", args = c("-v", "gffread")) == 1) stop("gffread NOT FOUND!!!")
# Preparation bam and gff3
## subset bam in specified region
prefix <- argv$prefix
outBam <- glue("{outdir}/{prefix}.bam")
cmd = glue("samtools view -h {argv$bam} {argv$region} | rg -v 'SN:KI|SN:GL'| samtools view -bS - > {outBam} && samtools index {outBam}")
print(cmd)
system(cmd)
## convert gff to gff3
refGff <- glue("{outdir}/{prefix}.MANEv1.gff3")
cmd2 <- glue("gffread --keep-genes {argv$mane} -r {argv$region} -o {refGff}")
print(cmd2)
system(cmd2)
outGff <- glue("{outdir}/{prefix}.gff3")
cmd3 <- glue("gffread --keep-genes {argv$gff} -r {argv$region} -o {outGff}")
print(cmd3)
system(cmd3)
if (system2("command", args = c("-v", "rmats2sashimiplot")) == 1) stop("rmats2sashimiplot NOT FOUND!!!")
region2 <- str_replace(str_replace(argv$region,":",paste0(":",argv$strand,":")),"-",":")
cmd4 <- glue("rmats2sashimiplot --b1 {outBam} -c {region2}:{outGff} --l1 {prefix} --exon_s 1 --intron_s 5 -o {outdir}/sashimi")
print(cmd4)
system(cmd4)

#idxTrack <- IdeogramTrack(genome="hg38", chromosome=chrN)
#axisTrack <- GenomeAxisTrack()

posInfo <- unlist(strsplit(argv$region, ':|-'))
chrN <- posInfo[1]
start <- as.numeric(posInfo[2])
end <- as.numeric(posInfo[3])
# reference
refTrack <- GeneRegionTrack(refGff, name = "MANE", fill = "#8282d2", stackHeight = 0.1,
                             transcriptAnnotation = "transcript",
)
relations <- read.table(refGff, comment.char = "#", header = F, sep = "\t") %>%
  dplyr::filter(V3=="gene") %>%
  separate_wider_delim(V9, ";", names = c("transcript", "gene_name")) %>%
  mutate(transcript = str_replace(transcript, "ID=","")) %>%
  mutate(gene_name = str_replace(gene_name, "gene_name=",""))
transcript(refTrack) <- tibble(transcript = transcript(refTrack)) %>% left_join(relations, by="transcript") %>%
  mutate(trans=if_else(startsWith(transcript,"ENSG"), gene_name, transcript)) %>% .$trans

# transcripts
transTrack <- GeneRegionTrack(outGff, name = "Transcripts", color = "#960000", fill = "#960000",
                           geneSymbols = TRUE, transcriptAnnotation = "transcript", stackHeight = 0.2
)

# reads
readTrack <- AlignmentsTrack(outBam, isPaired = FALSE, name = "Reads", stacking = "full",
                             type = c("coverage", "sashimi", "pileup"),
                             coverageHeight = 0.08, sashimiHeight = 0.08,
                             min.height = 0, minCoverageHeight = 0, minSashimiHeight = 0)

#plotTracks(c(idxTrack, axisTrack, refTrack, transTrack, readTrack), type = c("coverage", "sashimi", "pileup"),
plotTracks(c(refTrack, transTrack, readTrack),
           chromosome = chrN, from = start, to = end,
           sizes = c(1,4,10), title.width = 0.5, margin = 1,
           background.panel = "#FFFEDB",
           background.title = "darkblue"
)
