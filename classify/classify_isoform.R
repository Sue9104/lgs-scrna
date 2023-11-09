library(plyranges)
library(stringr)

gtf <- "~/projects/seq3/20230530/gtf/RIF.bambu_collapse.gtf"

bambu_gtf <- read_gff(gtf) 

bambu_gtf %>% 
  arrange(gene_id, transcript_id) %>% 
  mutate(transcript_id = paste0(stringr::str_split_i(gene_id, '\\.',1), '.', transcript_id)) %>% 
  write_gff2("~/projects/seq3/20230530/6_as_events/pigeon_classify/RIF.bambu_collapse.retrans.gtf")

bambu <- readRDS("~/projects/seq3/20230530/3_collapse/3_2_bambu_by_CB/default/RIF.bamsboo.discovery.rds")

obj.genes <- readRDS("~/projects/seq3/20230530/4_scrna/counts/7_final_cluster/RIF.genes.cluster.rds")
obj.genes@meta.data %>% filter(orig.ident == "Pro3") %>% lengths


