# Generate heatmap plot for DPA Analysis
#
# WARNING: This script is intended to be used directly by the tappAS application ONLY
#          Some of the file names are hard coded to avoid having to pass all names as arguments
#
# Script arguments:
#   input DPA folder: DPA result matrix
#   input matrix: DPA significance level
#   output plot: full path to output plot image file
#
# Note:   All arguments are required. Expression matrix must be in raw counts.
#
# Written by Pedro Salguero - Ángeles Arzalluz

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
# handle command line arguments - don't want to use external package so users don't have to install it
library(argparser)
p <- arg_parser("Differential PolyAdenylation Plot")
p <- add_argument(p, "indir", help="input directory", type="character")
p <- add_argument(p, "outfile", help="outdir directory", type="character")
p <- add_argument(p, "--sig", help="p value thresold", default=0.05)
argv <- parse_args(p)

indir = argv$indir
outfile = argv$outfile
sig = argv$sig
if(nchar(outfile) == 0 || nchar(indir) == 0 || nchar(sig) == 0)
  stop("Missing command line argument.")

# user-defined significance level
sig <- as.numeric(sig)

# internal tappAS file with DPA results
dpa_genes <- read_tsv(file.path(indir, "result.tsv")) %>% filter(qValue <= sig) %>% arrange(desc(totalChange))

# tappAS-retrieved DPAU levels for each sample
dpau <- read_tsv(file.path(indir, "dpa_DPAU_ByCondition.tsv"))
cn <- colnames(dpau)[2:(length(colnames(dpau))-1)]
dpau <- dpau %>%
    filter(Gene %in% dpa_genes$gene) %>%
    gather(key = CellType, value = DPAU, cn, -Gene, -Switching) %>%
    mutate(DPAU = DPAU / 100) %>%
    mutate(Switching_factor = ifelse(Switching, "Switching\n", "No Switching\n")) %>%
    arrange(CellType, desc(DPAU)) %>%
    mutate(Gene = ordered(Gene) %>% fct_reorder2(desc(CellType), DPAU))

# HEATMAP
# calculate dendrogram
dpau_mat <- dpau %>%
    select(Gene, CellType, DPAU) %>%
    spread(CellType, DPAU) %>%
    column_to_rownames("Gene") %>%
    as.matrix

dpau_dendro <- hclust(dist(dpau_mat)) %>% as.dendrogram
heatmap_order <- order.dendrogram(dpau_dendro)

dpau <- dpau %>%
mutate(Gene = ordered(Gene) %>%
fct_relevel(rownames(dpau_mat)[heatmap_order]))

# PLOT HEATMAP
ggp <- ggplot(dpau) + geom_tile(aes(x = Gene, y = CellType, fill = DPAU), color = "black") +
    ylab("Cell type") + xlab("DPA genes") +
    scale_fill_gradient2("DPAU", low = "#1e8780", mid = "#fdffba", high = "#a23427", midpoint = 0.5) +
    theme(plot.title = element_text(face = "plain"),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(), strip.background = element_rect(fill = "white"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.1,0.5),
        legend.direction = "horizontal",
        legend.text = element_text(angle = 45, size=8)) +
    facet_grid(cols = vars(Switching_factor), scales = "free_x", space = "free_x")

legend <- get_legend(ggp)

ggp_nolegend <- ggplot(dpau) + geom_tile(aes(x = Gene, y = CellType, fill = DPAU)) +
  ylab("Cell type") + xlab("DPA genes") +
  scale_fill_gradient2("DPAU", low = "#1e8780", mid = "#fdffba", high = "#a23427", midpoint = 0.5) +
  theme(plot.title = element_text(face = "plain"),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6, hjust = 1),
        axis.ticks.x = element_blank(), strip.text = element_text(size = 18),
        axis.line = element_blank(), strip.background = element_rect(fill = "white"),
        plot.subtitle = element_text(hjust = 0.5), legend.position = "none") +
  facet_grid(cols = vars(Switching_factor), scales = "free_x", space = "free_x")

ggp <- plot_grid(ggp_nolegend, legend, nrow = 2, rel_heights = c(0.9, 0.1))

ggsave(outfile, plot=ggp, device="png", bg="white", width=20, height=5, dpi=150)
cat("\nAll done.\n")
