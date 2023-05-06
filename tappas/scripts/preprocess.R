library(argparser)
library(tidyverse)
library(sqldf)
p <- arg_parser("Generate Matrix for transcript, gene and protein")
p <- add_argument(p, "gff3", help="transcript annotation gff3", type="character")
p <- add_argument(p, "design", help="group design tsv", type="character")
p <- add_argument(p, "raw", help="input raw transcript matrix tsv", type="character")
p <- add_argument(p, "infile", help="input filtered transcript matrix tsv", type="character")
p <- add_argument(p, "outdir", help="output directory")
argv <- parse_args(p)
if (!dir.exists(argv$outdir)) {dir.create(argv$outdir, recursive = TRUE) }
design <- read.table(argv$design, sep='\t', header=T)
samples <- as.vector(design$sample)
groups <- as.vector(unique(design$group))


# extract info from annotation
annotation <- read.csv.sql(
  argv$gff3, sep="\t", header = T,
  sql="select * from file where type in ('gene','transcript', 'protein', 'CDS',
      'polyA_Site', 'PAS', '3UTR_Length', '5UTR_Length')")
print(head(annotation))
gene.tib <- annotation %>% filter(type == "gene") %>%
  mutate(gene=str_replace(str_extract(attributes,"ID=\\w+"), "ID=","")) %>%
  select(seqid, gene)
head(gene.tib)
protein.tib <- annotation %>% filter(type == "protein") %>% mutate(proteinLength=end) %>%
  mutate(protein=str_replace(str_replace(str_extract(attributes,"ID=.*?;"), "ID=",""),';','')) %>%
  select(seqid, protein, proteinLength)
head(protein.tib)
trans.tib <- annotation %>% filter(type == "transcript") %>% mutate(transLength=end)  %>%
  mutate(category=str_replace(str_extract(attributes,"primary_class=\\w+"), "primary_class=",""))  %>%
  select(seqid, category, transLength, strand)
head(trans.tib)
polya.tib <- annotation %>% filter(type == "polyA_Site") %>%
  mutate(polyAPos=end) %>%
  select(seqid, polyAPos)
head(polya.tib)
pas.tib <- annotation %>% filter(type == "PAS") %>%
  distinct(seqid, .keep_all = TRUE) %>%
  mutate(PASstart = start) %>% mutate(PASend = end) %>%
  mutate(LengthPAS = as.numeric(end) - as.numeric(start)) %>%
  select(seqid, PASstart, PASend, LengthPAS)
head(pas.tib)
utr5.tib <- annotation %>% filter(type == "5UTR_Length") %>%
  mutate(Length5 = as.numeric(end) - as.numeric(start)) %>%
  select(seqid, Length5)
utr3.tib <- annotation %>% filter(type == "3UTR_Length") %>%
  mutate(Length3=as.numeric(end) - as.numeric(start)+1) %>%
  select(seqid, Length3)
head(utr3.tib)
cds.tib <- annotation %>% filter(type == "CDS", source == "tappAS")  %>%
    mutate(LengthCDS= if_else(start == '.', 0, as.numeric(end) - as.numeric(start) )) %>%
    select(seqid, LengthCDS)
head(cds.tib)
# headers: transcript	category	transLength	gene protein	proteinLength	polyAPos Length3	Length5	LengthCDS
anno <- trans.tib %>% left_join(gene.tib, by='seqid') %>%
  left_join(protein.tib, by='seqid') %>%
  left_join(polya.tib, by='seqid') %>%
  left_join(pas.tib, by='seqid') %>%
  left_join(utr3.tib, by='seqid') %>%
  left_join(utr5.tib, by='seqid') %>%
  left_join(cds.tib, by='seqid') %>% rename(transcript=seqid)
print(head(anno))
write.table(anno, file.path(argv$outdir, "trans_info.tsv"),
            sep="\t", quote = F, row.names = F)
# dpa_info.tsv
write.table(anno %>% filter(polyAPos != "") %>% select(gene, transcript, polyAPos, strand),
            file.path(argv$outdir, "dpa_info.tsv"),
            sep="\t", quote = F, row.names = F,
            col.names = c("Gene",	"Trans", "GenPos", "Strand"))
# structural_info.tsv
write.table(anno %>% select(transcript, Length3, Length5, LengthCDS, polyAPos, transLength),
            file.path(argv$outdir, "structural_info.tsv"),
            sep="\t", quote = F, row.names = F,
            col.names = c("#SeqName",	"Length3", "Length5", "LengthCDS","PosPAS","TotalLength"))
# result_gene_trans.tsv
write.table(anno %>% select(gene, transcript),
            file.path(argv$outdir, "result_gene_trans.tsv"),
            sep="\t", quote = F, row.names = F, col.names = F)
# gene_transcripts.tsv
write.table(anno %>% select(gene, transcript, transLength),
            file.path(argv$outdir, "gene_transcripts.tsv"),
            sep="\t", quote = F, row.names = F,
            col.names = c("geneName", "transcript", "length"))
# gene_proteins.tsv
write.table(anno %>% filter(protein != "") %>% mutate(proteinP=paste0(gene,'_',protein)) %>% select(gene, proteinP, proteinLength),
            file.path(argv$outdir, "gene_proteins.tsv"),
            sep="\t", quote = F, row.names = F,
            col.names = c("geneName", "protein", "length"))

GenerateMatrixFile <- function(infile, prefix = "") {
  expMatrix <- read.table(infile, row.names=1, sep="\t", header=TRUE)
  expMatrix <- expMatrix[,samples]
  write.table(expMatrix, file.path(argv$outdir, paste0("transcript_matrix",prefix,".tsv")),
              sep="\t", quote = F, row.names = T)
  mean.mat <- sapply(groups, function(x) {
                      s = design$sample[design$group == x];
                      rowMeans(expMatrix[,s])})
  colnames(mean.mat) <- groups
  head(mean.mat)
  write.table(mean.mat, file.path(argv$outdir, paste0("transcript_mean_matrix",prefix,".tsv")),
              sep="\t", quote = F, row.names = T)
  mean.tib <- as.data.frame(mean.mat) %>% rownames_to_column('transcript') %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write.table(anno %>% right_join(mean.tib, by = 'transcript') ,
              file.path(argv$outdir, paste0("transcript_total_info",prefix,".tsv")),
              sep="\t", quote = F, row.names = F)
  # utrl_info.tsv
  utrl_tmp.tib <- select(anno, gene,transcript,Length3,Length5,strand) %>% filter(Length3 != "")
  utrl_info.tib <- utrl_tmp.tib %>% right_join(mean.tib, by = 'transcript')
  head(utrl_info.tib)
  write.table(utrl_info.tib,
              file.path(argv$outdir, paste0("utrl_info",prefix,".tsv")),
              sep="\t", quote = F, row.names = F,
              col.names = c(c("Gene", "Trans", "UTR3", "UTR5", "Strand"), as.vector(groups)))

  gene_trans_protein <- anno %>% column_to_rownames('transcript') %>% select(gene, protein)
  data <- merge(expMatrix, gene_trans_protein, by = 0, all.x = TRUE) %>%
    column_to_rownames('Row.names')
  head(data)
  gene.df <- data %>% select(-protein)  %>% group_by(gene) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames('gene')
  write.table(gene.df, file.path(argv$outdir, paste0("gene_matrix",prefix,".tsv")),
              sep="\t", quote = F, row.names = T)

  protein.df <- data %>% filter(protein != "") %>%
    select(-gene) %>% group_by(protein) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames('protein')
  write.table(protein.df, file.path(argv$outdir, paste0("protein_matrix",prefix,".tsv")),
              sep="\t", quote = F, row.names = T)
  protein2.df <- data %>% filter(protein != "")  %>% group_by(gene,protein) %>%
    summarise(across(everything(), sum)) %>%
    mutate(new = paste0(gene,'_', protein)) %>% column_to_rownames('new') %>%
    select(-gene, -protein)
  write.table(protein2.df, file.path(argv$outdir, paste0("gene_protein_matrix",prefix,".tsv")),
              sep="\t", quote = F, row.names = T)
}
GenerateMatrixFile(argv$raw, prefix="_raw")
GenerateMatrixFile(argv$infile, prefix="")

