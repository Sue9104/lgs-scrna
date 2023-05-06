library(argparser)
library(tidyverse)
p <- arg_parser("Generate Matrix for transcript, gene and protein")
p <- add_argument(p, "design", help="group design tsv", type="character")
p <- add_argument(p, "outdir", help="output directory")
argv <- parse_args(p)
if (!dir.exists(argv$outdir)) {dir.create(argv$outdir, recursive = TRUE) }
design <- read.table(argv$design, sep='\t', header=T)
design[["group"]] <- as.numeric(factor(design$group))
# exp_factors.txt
write.table(design %>% column_to_rownames('sample'),
            file.path(argv$outdir, "exp_factors.txt"),
            sep="\t", quote = F, row.names = T,
            col.names = c("Replicate"))


