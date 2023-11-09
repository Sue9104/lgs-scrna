library(DropletUtils)
library(Matrix)
library(glue)
library(tidyverse)
library(tictoc)

library(argparser)
p <- arg_parser("Convert Bambu Results to 10X MTX format")
p <- add_argument(p, "indir", type="character",
                  help="input bambu result directory, including counts_gene|transcript.txt")
p <- add_argument(p, "outdir", help="output directory")
p <- add_argument(p, "--gtf", help="genome gtf used in bambu",
                  default=normalizePath("~/genomes/hg38/gencode/gencode.v43.annotation.gtf"))
argv <- parse_args(p)

# read GTF
tic("Reading annotation")
cmd = glue("grep -P '\\tgene\\t' {argv$gtf} | cut -f9| cut -d';' -f1,3|sed 's|^gene_id ||g;s|; gene_name |,|g' ")
gencodeSymbol <- read.csv(pipe(cmd), header = F, col.names = c("gene_id", "gene_name"), row.names = 1)
message("gencode data size (Mb): ", pryr::object_size(gencodeSymbol) / 1024 / 1024 )
toc()

indir <- normalizePath(argv$indir)
outdir <- argv$outdir
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

levels <- c("counts", "fullLengthCounts", "uniqueCounts")
for (l in levels) {
    message("level: ", l)
    tic("Reading isoform expression")
    dataTrans <- read.table(glue("{indir}/{l}_transcript.txt"), header = T, row.names = 1)
    toc()
    samples <- unique(gsub(pattern = "\\.CB_\\w+", replacement = "", x = colnames(dataTrans)[-1]))
    message("samples: ", paste(samples, collapse = " "))
    tran_symbol <- gencodeSymbol[dataTrans[,"GENEID"], "gene_name"]
    tran_symbol[is.na(tran_symbol)] <- row.names(dataTrans)[which(is.na(tran_symbol))]
    message("transcript data size (Mb): ", pryr::object_size(dataTrans) / 1024 / 1024 )
    message("Memory used (Mb): ", pryr::mem_used()/1024/1024)
    tic("Write 10X in isoform level")
    for (sample in samples) {
      message(sample)
      sampledir <- glue("{outdir}/{l}/{sample}/")
      if (!dir.exists(sampledir)) { dir.create(sampledir, recursive = TRUE) }
      sample_data <- as.matrix(dataTrans[,grepl(sample, colnames(dataTrans))])
      colnames(sample_data) <- sub("\\w+.CB_", "", colnames(sample_data))
      sampleoutdir <- glue("{sampledir}/isoforms")
      unlink(sampleoutdir, recursive = TRUE)
      message("sample data dim: ", paste(dim(sample_data), collapse=" "))
      write10xCounts(path = sampleoutdir,
                     as(sample_data, "sparseMatrix"),
                     barcodes = colnames(sample_data),
                     gene.id = rownames(sample_data),
                     gene.symbol = tran_symbol)
    }
    toc()


    #dataGenes <- read.table(file.path(indir,"counts_gene.txt"), header = T, row.names = 1)
    dataGenes <- dataTrans %>% group_by(GENEID) %>% summarise(across(everything(), list(sum))) %>% column_to_rownames("GENEID")
    rm(dataTrans)
    message("gene data size (Mb): ", pryr::object_size(dataGenes) / 1024 / 1024 )
    gene_symbol <- gencodeSymbol[row.names(dataGenes), "gene_name"]
    gene_symbol[is.na(gene_symbol)] <- row.names(dataGenes)[which(is.na(gene_symbol))]
    message("Memory used (Mb): ", pryr::mem_used()/1024/1024)
    tic("Write 10X in gene level")
    for (sample in samples) {
      sampledir <- glue("{outdir}/{l}/{sample}/")
      if (!dir.exists(sampledir)) { dir.create(sampledir, recursive = TRUE) }
      sample_data <- as.matrix(dataGenes[,grepl(sample, colnames(dataGenes))])
      colnames(sample_data) <- sub("\\w+.CB_", "", colnames(sample_data))
      message("sample data dim: ", paste(dim(sample_data), collapse=" "))
      sampleoutdir <- glue("{sampledir}/genes")
      unlink(sampleoutdir, recursive = TRUE)
      write10xCounts(path = sampleoutdir,
                     as(sample_data, "sparseMatrix"),
                     barcodes = colnames(sample_data),
                     gene.id = rownames(sample_data),
                     gene.symbol = gene_symbol)
    }
    toc()
    message("Done: ", l, "\n")
}

