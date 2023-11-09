library(tidyr)
library(dplyr)
library(glue)
library(ggplot2)

obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.20230906.rds")
obj.trans <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.trans.20230906.rds")
genes_fine <- read.table(glue::glue("{indir}/genes_fine/celltype.genes_fine.tsv"), 
                         sep="\t", row.names = 1, header = T)
trans_fine <- read.table(glue::glue("{indir}/trans_fine/celltype.trans_fine.tsv"), 
                         sep="\t", row.names = 1, header = T)
trans_fine %>% head
##############################################
## pyscenic 
##############################################


indir <- "/public/home/msu/projects/seq3/20230530/8_regulons"
auc <- read.csv(glue::glue("{indir}/genes_fine/auc_mtx_Z.csv"), row.names = 1)



rss_genes <- readr::read_csv(glue::glue("{indir}/genes_fine/cellType-regulon_specificity_scores.csv")) %>% 
  rename(cell_type = `...1`) %>% 
  tidyr::gather(-cell_type, key="Regulon", value = "RSS") %>% 
  group_by(cell_type) %>% arrange(desc(RSS)) %>% 
  mutate(Regulon2 = factor(interaction(cell_type, Regulon, drop=TRUE))) 
rss_genes %>% 
  #filter(cell_type == "Macrophages") %>%
  ggplot(aes(x = reorder(Regulon, -RSS), y = RSS)) + geom_point() + 
  ggrepel::geom_label_repel(data = rss_genes %>% group_by(cell_type) %>% 
               arrange(desc(RSS)) %>% 
               group_modify( ~head(.x, 10L)),
             aes(label = Regulon), size = 3, direction = "both") + 
  facet_wrap(cell_type~., scales = "free_x", drop = TRUE, ncol = 4) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))


rss_trans <- readr::read_csv(glue::glue("{indir}/trans_fine/cellType-regulon_specificity_scores.csv")) %>% 
  rename(cell_type = `...1`) %>% 
  tidyr::gather(-cell_type, key="Regulon", value = "RSS") %>% 
  group_by(cell_type) %>% arrange(desc(RSS)) 
rss_trans %>% 
  #filter(cell_type == "Macrophages") %>%
  ggplot(aes(x = reorder(Regulon, -RSS), y = RSS)) + geom_point() + 
  ggrepel::geom_label_repel(data = rss_trans %>% group_by(cell_type) %>% 
                              arrange(desc(RSS)) %>% 
                              group_modify( ~head(.x, 10L)),
                            aes(label = Regulon), size = 3, direction = "both") + 
  facet_wrap(cell_type~., scales = "free_x", drop = TRUE, ncol = 4) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))


top5_regulons <- rss_genes %>% group_by(cell_type) %>% 
  arrange(desc(RSS)) %>% 
  group_modify( ~head(.x, 5L)) %>% distinct(Regulon) %>% .$Regulon
top5_regulons

# Z-score
auc_mtx_Z <- read.csv(glue::glue("{indir}/genes_fine/auc_mtx_Z.csv"), row.names = 1)
celltypes = data.frame(
  genes = genes_fine[rownames(auc_mtx_Z), "cell_type"],
  trans = trans_fine[rownames(auc_mtx_Z), "cell_type"],
  row.names = rownames(auc_mtx_Z)
)
library(RColorBrewer)
cellnames <- sort(unique(c(celltypes$genes, celltypes$trans)))
celltype_colors <- colorspace::rainbow_hcl(length(cellnames), l = 80, c = 80)
names(celltype_colors) <- cellnames
celltype_colors 

### all not significant
pheatmap::pheatmap(auc_mtx_Z, 
                   annotation_row = celltypes, 
                   show_rownames = F, show_colnames = F,
                   cluster_rows = F, cluster_cols = F,
                   annotation_colors = list(genes = celltype_colors, trans = celltype_colors))  

auc_mtx <- read.csv(glue::glue("{indir}/genes_fine/auc_mtx.normal.csv"), row.names = 1)
total_mean <- auc_mtx %>% summarise(across(where(is.numeric), ~ mean(.x)))
total_sd <- auc_mtx %>% summarise(across(where(is.numeric), ~ sd(.x)))

auc_mtx_z_celltype <- auc_mtx %>% tibble::rownames_to_column("cell") %>% 
  mutate(genes_fine = genes_fine[cell, "cell_type"]) %>% 
  group_by(genes_fine) %>% 
  summarise(across(where(is.numeric), ~ mean(.x)))
total_mean <- total_mean[rep(1, nrow(auc_mtx_z_celltype)),]
rownames(total_mean) <- seq(nrow(auc_mtx_z_celltype))
total_sd <- total_sd[rep(1, nrow(auc_mtx_z_celltype)),]
rownames(total_sd) <- seq(nrow(auc_mtx_z_celltype))

auc_mtx_z_celltype[, -1] <- (auc_mtx_z_celltype[,-1] - total_mean) / total_sd
auc_mtx_z_celltype %>% gather(-genes_fine, "Regulon", "Z") %>% spread(genes_fine, Z)
  filter(if_any(where(is.numeric), ~ .>3))


