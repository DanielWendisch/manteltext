
library(tidyverse)
library(readxl)
library(Seurat)
library(janitor)
setwd("C:/Users/danne/R_projects/manteltext_git/scripts")

path_data_file <- "C:\\Users\\danne\\raw_data\\schulte_schrepping_et_al\\"
path_output_file <- "..\\output\\"
path_raw_data <- "C:/Users/danne/raw_data/"
seurat_c1_full_FG <- readRDS(paste0(path_raw_data,
                                    "schulte_schrepping_et_al/seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"))
Idents(seurat_c1_full_FG) <- "RNA_snn_res.0.4"

new_cluster_names<-
  c(
    `0`="Classical Monocytes",
    `7`="HLA DRhi CD83 Monocytes",
    `12`="CD163hi Monocytes 1_d7",
    `6`="HLA DRlo S100A Monocytes",
    `8`="Non classical Monocytes",
    `11`="Neutrophils",  
    `13`="Immature Neutrophils",
    `17`="mDCs",
    `18`="pDCs",
    `1`="CD4 T",
    `3`="CD4 T",
    `10`="CD4 T",
    `2`="CD8 T",
    `15`="CD8 T",
    `9`="CD8 T",
    `4`="NK",
    `5`="B",
    `19`="B",
    `20`="B",
    `16`="Plasmablasts",
    `14`="Megakaryocyte",
    `21`="mixed",
    `22`="undefined"
  )
#make converter named vector
new_cluster_names_clean <- make_clean_names(new_cluster_names)
names(new_cluster_names_clean) <- names(new_cluster_names)
#create new metadata variable for cluster names
seurat_c1_full_FG$cluster_clean <- new_cluster_names_clean[seurat_c1_full_FG$RNA_snn_res.0.4 |> as.character()]

meta_data_all <- tibble(cluster=seurat_c1_full_FG$cluster_clean,
                        donor=seurat_c1_full_FG$donor,
                        purification= seurat_c1_full_FG$purification,
                        condition_per_sample=seurat_c1_full_FG$group_per_sample,
                        days_post_symptom_onset=seurat_c1_full_FG$days_after_onset,
                        disease_stage=seurat_c1_full_FG$disease_stage)


meta_data_condition_per_sample_donor_day <- meta_data_all %>%
  mutate(condition_per_sample_donor_day=paste0(
    condition_per_sample,"__0_",donor, "_",days_post_symptom_onset)) %>% 
  #combine t-cell and b-cell subclusters
  mutate(cluster=str_replace_all(cluster,
                                 c( "cd8_t_3"="CD8_T_cells",
                                    "cd8_t_2"="CD8_T_cells",
                                    "cd8_t"="CD8_T_cells",
                                    "m_d_cs"= "mDCs",
                                    "cd4_t_3"="CD4_T_cells",
                                    "cd4_t_2"="CD4_T_cells",
                                    "cd4_t"="CD4_T_cells",
                                    "b_3"="B_cells",
                                    "b_2"="B_cells",
                                    "plasmablasts"="B_cells",
                                    "b"="B_cells",
                                    "p_d_cs"="pDCs",
                                    "hla_d_rhi_cd83_monocytes"="HLA_DRhi_CD83_Monocytes",
                                    "classical_monocytes" ="Classical_Monocytes",
                                    "cd163hi_monocytes_1_d7"="CD163hi_Monocytes_donor1_day1",
                                    "non_classical_monocytes"="Non_classical_Monocytes",
                                    "neutrophils" ="Neutrophils" ,
                                    "hla_d_rlo_s100a_monocytes"="HLA_DRlo_s100a_Monocytes",
                                    "megakaryocyte"="Megakaryocyte",
                                    "nk"="NK"
                                    
                                 )
  )) %>% 
  group_by(cluster,condition_per_sample_donor_day) %>% 
  summarise(n=n())%>% 
  ungroup()  %>% 
  pivot_wider(names_from = "cluster",values_from = n) %>%
  ungroup() %>% replace(is.na(.), 0) %>%
  rename(sample=condition_per_sample_donor_day) %>% 
  select(-mixed,-undefined)

meta_data_condition_per_sample_donor_day_grouping <- meta_data_condition_per_sample_donor_day %>% 
  select(sample) %>% 
  separate(col=sample, sep="__", into = c("condition","donor_day"),remove = FALSE) %>% 
  select(sample,condition) 

meta_data_condition_per_sample_donor_day %>% write_csv(paste0(path_output_file,"meta_data_condition_per_sample_donor_day.csv"))
meta_data_condition_per_sample_donor_day_grouping %>% write_csv(paste0(path_output_file,"meta_data_condition_per_sample_donor_day_grouping.csv"))


join_monos <- . %>% mutate(Monocytes=non_Classical_Monocytes+
                            HLA_DRlo_s100a_Monocytes+
                            HLA_DRhi_CD83_Monocytes+
                            CD163hi_Monocytes_donor1_day1+
                            Classical_Monocytes) %>% 
  select(-non_Classical_Monocytes,
         -HLA_DRlo_s100a_Monocytes,
         -HLA_DRhi_CD83_Monocytes,
         -CD163hi_Monocytes_donor1_day1,
         -Classical_Monocytes)


meta_data_condition_per_sample_donor_day %>%join_monos %>%
  write_csv(paste0(path_output_file,"meta_data_condition_per_sample_donor_day_monocytes_as_one.csv"))


















