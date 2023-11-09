library(argparser)
p <- arg_parser("DIU for ATSS, APA and exon")
p <- add_argument(p, "indir", help="input directory", type="character")
p <- add_argument(p, "outdir", help="output directory", type="character")
argv <- parse_args(p)
indir <- argv$indir
outdir <- argv$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

library(tidyr)
library(dplyr)
library(ggplot2)
library(UpSetR)
#######################################################
# Intergration
#######################################################
celltype <- basename(indir)
pdf(glue::glue("{outdir}/{celltype}_stat.pdf"))
db <- readr::read_table("~/pipelines/lgs-scrna/coupling/RIF.tss_tes.txt")
info <- readr::read_table(file.path(indir, "total/promoter_dominance.LATER.info.txt")) %>% 
  mutate(transcript_id = if_else(pairType=="novelPair", pairs_id, transcript_id)) %>% 
  separate_wider_delim(pairs_id, ":", names = c("gene", "promoter", "pas")) %>% 
  mutate(gene_id = gene, promoter_id = paste(gene_id, promoter, sep = ":"),
         tes_id=paste(gene_id, pas, sep = ":")) %>% 
  left_join(db, by="gene_id") %>% 
  mutate(promoter_type=case_when(promoter == distal_tss ~ "distal",
                                 promoter == proximal_tss ~ "proximal",
                                 TRUE ~ 'intermediate'),
        utr_type=case_when(pas==proximal_tes ~ "proximal",
                           pas==distal_tes ~ "distal",
                           TRUE ~ "intermediate"))
  
# gene counts: TSS and PAS
status_tss_tes <- info %>% group_by(gene_id) %>% 
  summarise(tss=length(unique(promoter)), tes=length(unique(pas))) %>% 
  mutate(category=case_when(tss==1 & tes==1 ~ "TSS-PAS",
                            tss>1 & tes==1 ~ "ATSS-PAS",
                            tss==1 & tes>1 ~ "TSS-APA",
                            tss>1 & tes>1 ~ "ATSS-APA"))
status_tss_tes
counts_tss_tes <- status_tss_tes %>% group_by(category) %>% summarise(n=n()) %>% 
  pivot_wider(names_from = category, values_from = n) %>% 
  mutate(celltype=celltype, .before = 1)
counts_tss_tes %>% readr::write_delim(glue::glue("{outdir}/{celltype}.genes_tss_tes.counts.tsv"))

# (TSS and Dominance Promoter) vs (PAS: proximal and distal)
domP <- readr::read_table(file.path(indir, "total/promoter_dominance.LATER.result.txt")) %>% 
  filter(tss.status == "ATSS", apa.status == "APA") %>%
  group_by(promoter_type, utr_type) %>% summarise(n=n()) %>% 
  mutate(category="domP", utr_type=if_else(utr_type=="other", "intermediate", utr_type))
domP
single_tss <- status_tss_tes %>% filter(tss==1) %>% 
  left_join(info, by="gene_id") %>% distinct(gene_id, promoter_type, utr_type) %>% 
  group_by(promoter_type, utr_type) %>% summarise(n=n()) %>% mutate(category="singleTSS")
single_tss
genes.single_domP_stats <- bind_rows(single_tss, domP) %>% 
  mutate(category=factor(category, levels = c("singleTSS","domP"))) %>%  
  pivot_longer(cols=c(promoter_type, utr_type), 
               names_to = c("type"), values_to = c("location")) %>% 
  group_by(category, type, location) %>% summarise(counts = sum(n)) %>% 
  pivot_wider(id_cols = c(category, type), names_from = location, values_from = counts) %>% 
  mutate(celltype=celltype, .before = 1) %>% mutate_if(is.numeric,coalesce,0) 
genes.single_domP_stats %>% 
  readr::write_delim(glue::glue("{outdir}/{celltype}.genes_tss_domP.location.counts.tsv"))

p <- bind_rows(single_tss, domP) %>% 
  mutate(category=factor(category, levels = c("singleTSS","domP"))) %>%  
  pivot_longer(cols=c(promoter_type, utr_type), 
               names_to = c("type"), values_to = c("location")) %>% 
  ggplot(aes(x=category, y = n, fill = location)) + 
  geom_bar(stat='identity', position="fill",width = .7) + 
  ggsci::scale_fill_aaas(alpha = 0.7) + facet_grid(~type) + 
  labs(x = "", y = "", title = glue::glue("{celltype}: single TSS vs DomP")) + 
  theme(aspect.ratio=2,axis.ticks = element_blank(),
        panel.grid=element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(0.6, "lines")) 
print(p)

# (PAS and Dominance PAS) vs (TSS: proximal and distal)
domPAS <- readr::read_table(file.path(indir, "total/pas_dominance.LATER.result.txt")) %>% 
  filter(tss.status == "ATSS", apa.status == "APA") %>%
  group_by(promoter_type, utr_type) %>% summarise(n=n()) %>% 
  mutate(category="domPAS", utr_type=if_else(utr_type=="other", "intermediate", utr_type))
domPAS
single_tes <- status_tss_tes %>% filter(tes==1) %>% 
  left_join(info, by="gene_id") %>% distinct(gene_id, promoter_type, utr_type) %>% 
  group_by(promoter_type, utr_type) %>% summarise(n=n()) %>% mutate(category="singleTES")
single_tes
genes.single_domPAS_stats <- bind_rows(single_tes, domPAS) %>% 
  mutate(category=factor(category, levels = c("singleTES","domPAS"))) %>%  
  pivot_longer(cols=c(promoter_type, utr_type), 
               names_to = c("type"), values_to = c("location")) %>% 
  group_by(category, type, location) %>% summarise(counts = sum(n)) %>% 
  pivot_wider(id_cols = c(category, type), names_from = location, values_from = counts) %>% 
  mutate(celltype=celltype, .before = 1) %>% mutate_if(is.numeric,coalesce,0) 
genes.single_domPAS_stats %>% 
  readr::write_delim(glue::glue("{outdir}/{celltype}.genes_tes_domPAS.location.counts.tsv"))

p <- bind_rows(single_tes, domPAS) %>% 
  mutate(category=factor(category, levels = c("singleTES","domPAS"))) %>% 
  pivot_longer(cols=c(promoter_type, utr_type), 
               names_to = c("type"), values_to = c("location")) %>% 
  ggplot(aes(x=category, y = n, fill = location)) + 
  geom_bar(stat='identity', position="fill",width = .7) + 
  ggsci::scale_fill_aaas(alpha = 0.7) + facet_grid(~type) + 
  labs(x = "", y = "", title = glue::glue("{celltype}: single TES vs DomPAS")) + 
  theme(aspect.ratio=2,axis.ticks = element_blank(),
        panel.grid=element_blank(), axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.spacing = unit(0.6, "lines")) 
print(p)

#######################################################
# Significant Links
# domP for TSS-PAS link
# exon-link for AS link
#######################################################
genes.domP <- readr::read_table(file.path(indir, "total/promoter_dominance.LATER.result.txt")) %>% 
  distinct(gene_id) %>% mutate(significant=1, significant_level = "DomP")
genes.domPAS <- readr::read_table(file.path(indir, "total/pas_dominance.LATER.result.txt")) %>% 
  distinct(gene_id) %>% mutate(significant=1, significant_level = "DomPAS")
genes.tss_exon_links <- readr::read_table(file.path(indir, "diu/promoter_diu.celltype.tsv")) %>% 
  filter(DIU==1) %>% mutate(significant=1, significant_level = "exon_links") %>% 
  select(-q.value, -DIU)
genes.exon_pas_links <- readr::read_table(file.path(indir, "diu/pas_diu.celltype.tsv")) %>% 
  filter(DIU==1) %>% mutate(significant=1, significant_level = "exon_links") %>% 
  select(-q.value, -DIU)

genes.atss <- filter(status_tss_tes, tss > 1) %>% select(gene_id) %>% 
  mutate(category="ATSS-AS genes", ATSS=1) %>% 
  left_join(genes.tss_exon_links, by="gene_id") %>% 
  mutate_if(is.numeric,coalesce,0)
genes.apa <- filter(status_tss_tes, tes > 1) %>% select(gene_id) %>% 
  mutate(category="AS-APA genes", APA=1) %>% 
  left_join(genes.exon_pas_links, by="gene_id") %>% 
  mutate_if(is.numeric,coalesce,0)
genes.atss_apa <- filter(status_tss_tes, tss > 1, tes > 1) %>% select(gene_id) %>% 
  mutate(category="ATSS-APA genes", ATSSAPA=1) %>% 
  left_join(genes.domP, by="gene_id") %>% 
  mutate_if(is.numeric,coalesce,0)
# genes.atss_apa_1 <- filter(status_tss_tes, tss > 1, tes > 1) %>% select(gene_id) %>% 
#   mutate(category="ATSS-APA genes ", ATSSAPA=1) %>% 
#   left_join(genes.domPAS, by="gene_id") %>% 
#   mutate_if(is.numeric,coalesce,0)

# significant links in atss, apa genes
genes.stat <- bind_rows(genes.atss, genes.apa, genes.atss_apa) %>% 
  mutate_if(is.numeric,coalesce,0) %>% 
  replace_na(list(significant_level = "NO"))
p1 <- genes.stat %>%
  mutate(category=factor(category, levels = c("ATSS-APA genes", "AS-APA genes", "ATSS-AS genes"))) %>% 
  group_by(category, significant_level) %>% 
  summarise(n=n(), significant = sum(significant)) %>% 
  ggplot(aes(x=category, y = n, fill = significant_level)) + 
  geom_bar(stat='identity', position="fill",width = .7) + 
  ggsci::scale_fill_aaas(alpha = 0.7) + 
  labs(x = "", y = "Proportion of genes", 
       title = glue::glue("{celltype}: Sigificant Links")) + 
  theme(aspect.ratio=1.2,axis.ticks = element_blank(),
        panel.grid=element_blank(), axis.text.x = element_blank(),
        panel.spacing = unit(0.01, "lines")) + coord_flip()
print(p1)

# mutual significant links
genes.tss_pas_link <- select(genes.domP, gene_id) %>% mutate(category="TSS-PAS link") %>% 
  left_join(mutate(genes.tss_exon_links, significant_level = "TSS-exon"), by = "gene_id")
genes.tss_pas_link
genes.tss_pas_link_2 <- select(genes.domP, gene_id) %>% mutate(category="TSS-PAS link ") %>% 
  left_join(mutate(genes.exon_pas_links, significant_level = "exon-PAS"), by = "gene_id")
genes.tss_pas_link_2
genes.tss_exon_pas <- select(genes.tss_exon_links, gene_id) %>% mutate(category="TSS-exon link") %>% 
  left_join(mutate(genes.domP, significant_level = "Dom P"), by = "gene_id")
genes.tss_exon_pas
genes.exon_pas_link <- select(genes.exon_pas_links, gene_id) %>% mutate(category="exon-PAS link") %>% 
  left_join(mutate(genes.domP, significant_level = "Dom P"), by = "gene_id")
genes.exon_pas_link
genes.stat_2 <- bind_rows(genes.tss_pas_link, genes.tss_pas_link_2, genes.tss_exon_pas, genes.exon_pas_link) %>% 
  mutate_if(is.numeric,coalesce,0) %>% 
  replace_na(list(significant_level = "NO"))
p2 <- genes.stat_2 %>%
  mutate(category=factor(category, levels = c("exon-PAS link", "TSS-exon link", "TSS-PAS link ", "TSS-PAS link")),
         significant_level=factor(significant_level, levels = c("Dom P", "TSS-exon", "exon-PAS", "NO"))) %>% 
  group_by(category, significant_level) %>% 
  summarise(n=n(), significant = sum(significant)) %>% 
  ggplot(aes(x=category, y = n, fill = significant_level)) + 
  geom_bar(stat='identity', position="fill",width = .7) + 
  ggsci::scale_fill_aaas(alpha = 0.7) + 
  labs(x = "", y = "Proportion of genes", 
       title = glue::glue("{celltype}: Mutual Sigificant Links")) + 
  theme(aspect.ratio=1.2,axis.ticks = element_blank(),
        panel.grid=element_blank(), axis.text.x = element_blank(),
        panel.spacing = unit(0.01, "lines")) + coord_flip()

print(p2)
#######################################################
# Differential Isoform Usage at stage
#######################################################
genes.diu <- readr::read_table(file.path(indir, "diu/exon_diu.stage.tsv")) %>% 
  filter(DIU==1) %>% select(gene_id, DIU) 
genes.diu


#######################################################
# Intersection
#######################################################
genes.link <- mutate(genes.stat_2,  category = stringr::str_split_i(category, ' ', 1), n = 1) %>% 
  distinct(gene_id, category, n) %>% 
  pivot_wider(id_cols = gene_id, names_from = category, values_from = n) 
genes.total <- select(genes.stat, gene_id, ATSS, APA, ATSSAPA) %>% 
    full_join(genes.link, by="gene_id") %>% 
    left_join(genes.diu, by="gene_id") %>% 
    mutate_if(is.numeric,coalesce,0) %>% 
    rename(ATSS_APA=ATSSAPA, TSS_PAS_link = `TSS-PAS`, 
           TSS_exon_link = `TSS-exon`, exon_PAS_link = `exon-PAS`)
genes.total %>% readr::write_delim(glue::glue("{outdir}/{celltype}.genes_status.upsetR.tsv"))
p <- upset(as.data.frame(genes.total), order.by = "freq", empty.intersections = "on",
      sets = c("ATSS_APA", "TSS_PAS_link", "TSS_exon_link", "exon_PAS_link", "DIU"),
      queries = list(list(query = intersects, params = list("TSS_PAS_link",  "ATSS_APA"), active = T),
                     list(query = intersects, params = list("TSS_exon_link",  "exon_PAS_link"), active = T)))
print(p)
dev.off()


