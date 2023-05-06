# Differential Expression Analysis R script
#
# Differential expression analysis will be performed by using two different methods:
#   edgeR package and
#   NOISeq package

#### NOISeq differential expression - 'data' has already been normalized

de_noiseq <- function(data, factors, replicates=c("technical", "biological", "no")) {
  library("NOISeq")
  data = as.matrix(data)
  factors[,1]=as.factor(factors[,1]) #Must be always a factor.
  data_object =  readData(data = data, factors = factors)
  factor=colnames(factors)[1] ##Because we have only one factor defined. If more than one possible as input, we need to add a new argument to the function: factor
  QCreport(data_object, samples = NULL, factor = factor,
           file = file.path(outdir, paste0("QC.",dataType,".pdf")))
  if(replicates == "biological") {
       result = noiseqbio(data_object, norm="n", k=NULL, lc = 0 , factor = factor, filter=0)
  }
  else if(replicates == "technical") {
      result = noiseq(data_object, norm="n", k=NULL, lc = 0 , factor = factor, replicates=replicates)
  }
  else if(replicates == "no") {
      result = noiseq(data_object, norm="n", k=NULL, lc = 0 , factor = factor, replicates=replicates)
  }
  de <- degenes(result, q = 0.95)
  #print(de)
  DE.plot(result, q = 0.95, graphic = "expr", log.scale = TRUE)
  DE.plot(result, q = 0.95, graphic = "MD")

  # NEED chromosome position for readData
  #DE.plot(result, chromosomes = NULL, q = 0.8, graphic = "distr")
  #DE.plot(result, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE, q = 0.8, graphic = "chrom")
  # return results but remove theta column, not used
  if("theta" %in% colnames(result@results[[1]])){
    return (subset(result@results[[1]], select=-c(theta)))
  }else{
    return (subset(result@results[[1]], select=-c(M,D)))
  }
}

#### edgeR Differential Expression - 'data' contains raw counts

de_edgeR <- function(data, factors) {
  library(edgeR)
  tmm_factors = calcNormFactors(round(data))
  myedgeR = DGEList(counts = round(data),
                                  group = as.factor(factors[,1]),
                                  norm.factors = tmm_factors,
                                  remove.zeros = FALSE)
  myedgeR = estimateDisp(myedgeR)
  res = exactTest(myedgeR, pair=sort(levels(as.factor(myfactors[,1])), decreasing = T))
  res.sort = topTags(res, n = nrow(data))[[1]]
  return(res.sort)
}

#### validate command line arguments
library(argparser)
p <- arg_parser("Different Expression Analysis")
p <- add_argument(p, "indir", help="input directory", type="character")
p <- add_argument(p, "outdir", help="outdir directory", type="character")
p <- add_argument(p, "--dataType", help="gene, trans or protein", default="gene")
p <- add_argument(p, "--method", help="NOISEQ or EDGER", default='NOISEQ')
p <- add_argument(p, "--reps", help="biological or technical or no", default='biological')
argv <- parse_args(p)

dataType = argv$dataType
method = argv$method
reps = argv$reps
indir = argv$indir
outdir = argv$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if(nchar(outdir) == 0 || nchar(indir) == 0 || nchar(dataType) == 0 || nchar(reps) == 0 || nchar(method) == 0)
  stop("Missing command line argument.")

#### determine what type of data to run DEA for
deTrans = FALSE
deProteins = FALSE
deGenes = FALSE
if(dataType == "trans") {
    deTrans = TRUE
}
if(dataType == "protein") {
    deProteins = TRUE
}
if(dataType == "gene") {
    deGenes = TRUE
}

#### Read expression factors

cat("\nReading factors file data...")
myfactors=read.table(file.path(indir, "exp_factors.txt"), row.names=1, sep="\t", header=TRUE)
pdf(paste0(outdir, '/DEA.', dataType, '.pdf'))
#### check if DE for transcripts requested
if(deTrans) {
  #run DEA function
  cat("\nCalculating transcripts differential expression...\n")
  result_trans <- NULL
  if(method == "EDGER") {
    cat("\nReading raw counts transcript matrix file data...")
    expMatrix = read.table(file.path(indir, "transcript_matrix_raw.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " raw counts transcripts expression data rows")
    result_trans = de_edgeR(data=expMatrix, factors=myfactors)
  } else {
    cat("\nReading normalized transcript matrix file data...")
    expMatrix = read.table(file.path(indir, "transcript_matrix.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " normalized transcripts expression data rows")
    result_trans = de_noiseq(data=expMatrix, factors=myfactors, replicates=reps)
  }
  # add transcript (already in rownames but added for write.matrix)
  result_trans[,"transcript"] = rownames(result_trans)
  results = result_trans[complete.cases(result_trans[,1]),]

  cat("\nSaving transcript DEA results to file...")
  write.table(results, file.path(outdir, "result_trans.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
}

#### check if DE for proteins requested
if(deProteins) {
  #run DEA function
  cat("\nCalculating CDS differential expression...\n")
  result_proteins <- NULL
  if(method == "EDGER") {
    cat("\nReading raw counts CDS matrix file data...")
    expMatrix = read.table(file.path(indir, "protein_matrix_raw.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " raw counts CDS expression data rows")
    result_proteins = de_edgeR(data=expMatrix, factors=myfactors)
  } else {
    cat("\nReading normalized CDS matrix file data...")
    expMatrix = read.table(file.path(indir, "protein_matrix.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " normalized CDS expression data rows")
    result_proteins = de_noiseq(data=expMatrix, factors=myfactors, replicates=reps)
  }
  # add protein (already in rownames but added for write.matrix)
  result_proteins[,"protein"] = rownames(result_proteins)
  results = result_proteins[complete.cases(result_proteins[,1]),]

  cat("\nSaving proteins DEA results to file...")
  write.table(results, file.path(outdir, "result_protein.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
}

#### check if DE for genes requested
if(deGenes) {
  #run DEA function
  cat("\nCalculating genes differential expression...\n")
  result_genes <- NULL
  if(method == "EDGER") {
    cat("\nReading raw counts genes matrix file data...")
    expMatrix = read.table(file.path(indir, "gene_matrix_raw.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " raw counts genes expression data rows")
    result_genes = de_edgeR(data=expMatrix, factors=myfactors)
  } else {
    cat("\nReading normalized genes matrix file data...")
    expMatrix = read.table(file.path(indir, "gene_matrix.tsv"), row.names=1, sep="\t", header=TRUE)
    cat("\nRead ", nrow(expMatrix), " normalized genes expression data rows")
    result_genes = de_noiseq(data=expMatrix, factors=myfactors, replicates=reps)
  }
  # add gene (already in rownames but added for write.matrix)
  result_genes[,"gene"] = rownames(result_genes)
  results = result_genes[complete.cases(result_genes[,1]),]

  cat("\nSaving genes DEA results to file...")
  write.table(results, file.path(outdir, "result_gene.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
}
dev.off()
#### write completion file
cat("\nWriting DEAnalysis completed file...")
filedone <- file(file.path(outdir, paste("done_", dataType, ".txt", sep="")))
writeLines("end", filedone)
close(filedone)
cat("\nAll done.\n")
