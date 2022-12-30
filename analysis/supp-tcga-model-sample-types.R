load("../modeling/data/train_data_extra_info_v3.RData")

library(dplyr)
data = UCSCXenaShiny::tcga_clinical %>% 
  filter(sample %in% data_gene_load$sample) %>% 
  select(sample, type)
  

library(ggplot2)
p = ggplot(data, aes(x = forcats::fct_infreq(type))) +
  geom_bar() +
  labs(x = NULL, y = "Patient count") +
  ggthemes::theme_base() +
  ggpubr::rotate_x_text(30)
ggsave(filename = "plots/TCGA_modeling_sample_type_dist.pdf", p, width = 8, height = 4)
