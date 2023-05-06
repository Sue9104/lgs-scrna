# Generate scatter plot for DIU Analysis
#
# WARNING: This script is intended to be used directly by the tappAS application ONLY
#          Some of the file names are hard coded to avoid having to pass all names as arguments
#
# Script arguments:
#   input matrix: full path to matrix file
#   data type: TRANS, PROTEIN, GENE
#   output plot: full path to output plot image file
#
# Note:   All arguments are required. Expression matrix must be in raw counts.
#
# Written by Pedro Salguero - Lorena de la Fuente

library("ggplot2")
library("ggrepel")

# handle command line arguments - don't want to use external package so users don't have to install it
library(argparser)
p <- arg_parser("Differential Isoform Usage Density Plot")
p <- add_argument(p, "infile", help="input diu matrix, ", type="character")
p <- add_argument(p, "outfile", help="output pdf file", type="character")
p <- add_argument(p, "--dataType", help="gene, trans or protein", default="gene")
argv <- parse_args(p)

infile = argv$infile
outfile = argv$outfile
datatype = argv$dataType
if(nchar(outfile) == 0 || nchar(infile) == 0 || nchar(datatype) == 0)
  stop("Missing command line argument.")

# Read matrix
cat("\nReading data matrix...")
diu_matrix=read.table(file.path(infile), sep="\t", header=TRUE)
diu_matrix=diu_matrix[order(-diu_matrix$totalChange),]

myPalette = c("#f2635c","#7cafd6","#a2ca72","#f7c967","#f9f784","#a16cc1","#a796ff","#ff96eb")

# Create graph
ggp = ggplot(diu_matrix) +
  ggtitle("Gene Summary DIU Results") +
  geom_point(aes(x=log2(MeanExpression2/MeanExpression1), y=totalChange, color=podiumChange, shape=DS, size=DS)) +
  labs(x = "log2(FC)") +
  scale_x_continuous(limits = c(-7.5,7.5)) +
  scale_y_continuous(limits = c(0,max(diu_matrix$totalChange)*1.1), expand = c(0,0)) +

  theme_classic() +
  scale_shape_manual(values = c(19,1), name = "DIU") + scale_size_manual(values = c(2,2), name = "DIU") +
  scale_color_manual(values = myPalette[c(3,2)], name = "Mayor Switching") +
  ylab(label = "% Usage Change") +

  geom_vline(xintercept = 0.5, color=myPalette[1], size=1, linetype="dashed") + geom_vline(xintercept= -0.5, size=1, color=myPalette[1], linetype="dashed") +

  geom_text_repel(
    data = subset(diu_matrix[1:5,]), aes(x = log2(MeanExpression2/MeanExpression1), y = totalChange, label=Gene),
    #data = subset(podium.results2, gene %in% c("Klc1","Cadm1")), aes(-log2(q.value),log2FC,label=gene),
    size = 5,
    max.iter = 10000,
    label.padding = max(apply(diu_matrix,2,nchar)[,1])*2,
    box.padding = 1,
    point.padding = 0.3
  )+

  theme(plot.title = element_text(size=18,  hjust = 0.5, margin=margin(5,5,10,0)),
        axis.title.x = element_text(size=18, margin=margin(5,0,0,0)),
        axis.text.x  = element_text(margin=margin(7,0,0,0), size=16),
        axis.title.y = element_text(size=18,  margin=margin(0,15,0,0)),
        axis.text.y  = element_text(vjust=0.5, size=16)) +

  theme(legend.text = element_text(size = 15), legend.title = element_text(size=15), legend.key.size = unit(0.5, "cm"))

ggsave(outfile, plot=ggp, device="png", bg="transparent", width=10, height=7.5, dpi=120)
cat("\nAll done.\n")
