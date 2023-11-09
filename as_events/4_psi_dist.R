library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

indir <- "~/projects/seq3/20230530/6_as_events/pool_genes/"

celltypes <- c("Fib1", "Fib2", "Fib3")
categorys <- c("A5","A3","AF","AL","MX","RI","SE")
avg_psi.tib <- c()
# aspect.ratio = 1000/300
for (category in categorys){
  # category <- "AL"
  # prefix: all, qc, hvt, diff, diff2, diff3, diff2_up, diff2_down, diff3_up
  prefix <- "qc"
  trans_avail <- read.table(glue::glue("{indir}/../tpm/RIF.celltype_tpm.{prefix}.tsv"), header = T) %>% rownames()
  events <- read.table(glue::glue("{indir}/RIF.bambu_collapse.as_events_{category}_variable_10.ioe"), header = T) %>%   
    filter(alternative_transcripts %in% trans_avail) %>% 
    .$event_id
  
  avg_psi <- read.table(glue::glue("{indir}/RIF.{category}.celltype_{prefix}.psi")) %>% 
    tibble::rownames_to_column("event") %>% 
    ## filter only trans-related events
    # filter(event %in% events) %>%
    tidyr::gather("celltype", "psi", -event) %>% 
    filter(celltype %in% celltypes) %>% 
    mutate(category = category) %>% drop_na()
  avg_psi.tib <- bind_rows(avg_psi, avg_psi.tib)
}
groups = list(c("Fib1", "Fib2"),c("Fib1", "Fib3"),c("Fib2", "Fib3"))
avg_psi.tib %>% 
  ggplot(aes(y=psi, x=celltype, color = celltype)) + 
  geom_boxplot() + facet_grid(~category) +
  ggsignif::geom_signif(comparisons = groups, textsize = 2, color ="black",
                        test = "wilcox.test", 
                        # test.args = "greater",
                        map_signif_level = T,
                        tip_length = 0, step_increase = 0.1) +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(fill = NA),
        # aspect.ratio = 12/4,
        axis.title = element_text(size=10),
        strip.text = element_text(size=10, face = "bold"))

avg_psi.tib %>% filter(category == "A3") %>% 
  # filter(tran == "ENSG00000117650.13_and_ENSG00000224763.1;AL:chr1:211658657:211658712-211667106:211662772:211663652-211667106:-") %>% head()
  # pivot_wider(id_cols = c(tran), names_from = celltype, values_from = psi) 
  spread("celltype", "psi")  
  filter(Fib2 < 0.5 )




# not done yet, 20231031
sample_psi.tib <- c()
for(celltype in celltypes){
  # celltype <- "Fib2"
  for (category in categorys){
    # category <- "AL"
    # sample_psi <- read.table(glue::glue("{indir}/{celltype}.{category}_event.psi")) %>% 
    #   tibble::rownames_to_column("tran") %>% 
    #   tidyr::gather("sample", "psi", -tran) %>% 
    #   mutate(category = category, celltype = celltype) %>% drop_na()
    # sample_psi.tib <- bind_rows(sample_psi, sample_psi.tib)
  }
}





