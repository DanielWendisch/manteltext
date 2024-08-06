library(tidyverse)
library(Seurat)
library(janitor)

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

janitor_to_nice <- c( "cd8_t_3"="CD8_T_cells",
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
                      "nk"="NK",
                      "mixed"="mixed",
                      "undefined" ="undefined",
                      "immature_neutrophils"="immature_neutrophils"
                      
)


path_data_file <- "C:\\Users\\danne\\raw_data\\schulte_schrepping_et_al\\"
path_output_file <- ".\\output\\"
seurat_c1_full_FG <- readRDS("C:/Users/danne/raw_data/schulte_schrepping_et_al/seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds")
Idents(seurat_c1_full_FG) <- "RNA_snn_res.0.4"

new_cluster_names_clean <- make_clean_names(new_cluster_names)
names(new_cluster_names_clean) <- names(new_cluster_names)
#create new metadata variable for cluster names
seurat_c1_full_FG$cluster_clean <- new_cluster_names_clean[seurat_c1_full_FG$RNA_snn_res.0.4 |> as.character()]
seurat_c1_full_FG$cluster_combined <- janitor_to_nice[seurat_c1_full_FG$cluster_clean |> as.character()]


seurat_monocytes <-subset(x = seurat_c1_full_FG, idents = c("0","6","7","12","8"#, "11"#,"17"
))
seurat_monocytes |> write_rds(paste0(path_data_file,
                                     "schulte_schrepping_monocytes_from_seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"))                    

#monocyte&neutrophil subset

seurat_monocytes_neutrophils <-subset(x = seurat_c1_full_FG, idents = c("0","6","7","12","8", "11","13"))

seurat_monocytes_neutrophils |> write_rds(paste0(path_data_file,
                                                 "schulte_schrepping_monocytes_neutrophils_from_seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"))
