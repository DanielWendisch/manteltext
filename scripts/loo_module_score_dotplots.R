
library(effsize)
print("start leave one out module score scropt")


###
### only recalulate if files doesn't exist yet
if(!"scores_loo_corr.rds" %in% list.files(here::here("intermediate_data"))){
  print("not in abolute path, recalculate scores_loo_corr")
  print(getwd())

if (!("scores_loo_corr.rds" %in% (list.files("../intermediate_data/")))) {
  #construct list of genestes leaving out one gene each
  print("not in realtive path, recalculate scores_loo_corr")

Reyf_MCluster1 <- list()
Reyf_MCluster1[["full"]] <- genesets_fibrosis$Reyf_MCluster1

for (i in seq_along(Reyf_MCluster1[["full"]])) {
  Reyf_MCluster1[[genesets_fibrosis$Reyf_MCluster1[i]]] <- genesets_fibrosis$Reyf_MCluster1[-i]
  
}

Ayau_IPFeM <- list()
Ayau_IPFeM[["full"]] <- genesets_fibrosis$Ayau_IPFeM

for (i in seq_along(Ayau_IPFeM[["full"]])) {
  Ayau_IPFeM[[genesets_fibrosis$Ayau_IPFeM[i]]] <- genesets_fibrosis$Ayau_IPFeM[-i]
  
}

Reyf_Mfibrosis<- list()
Reyf_Mfibrosis[["full"]] <- genesets_fibrosis$Reyf_Mfibrosis

for (i in seq_along(Reyf_Mfibrosis[["full"]])) {
  Reyf_Mfibrosis[[genesets_fibrosis$Reyf_Mfibrosis[i]]] <- genesets_fibrosis$Reyf_Mfibrosis[-i]
  
}

Mors_1_SPP_fib <- list()
Mors_1_SPP_fib[["full"]] <- genesets_fibrosis$`Mors_1SPP1+MIPF`

for (i in seq_along(Mors_1_SPP_fib[["full"]])) {
  Mors_1_SPP_fib[[genesets_fibrosis$`Mors_1SPP1+MIPF`[i]]] <- genesets_fibrosis$`Mors_1SPP1+MIPF`[-i]
  
}

#takes loooookg!
#compute all module score scores of the leave one out genesets
#lo_corr means "leave one out corrrected"
BAL_macrophages <- add_named_module_score_ctrl_corrected(BAL_macrophages,
                                                         score_id = "Reyf_MCluster1_lo_corr",
                                                         genesets_list_named = lapply(Reyf_MCluster1, function(x){ensemblizer(x)}) )

BAL_macrophages <- add_named_module_score_ctrl_corrected(BAL_macrophages,
                                                         score_id = "Reyf_fibr_lo_corr",
                                                         genesets_list_named = lapply(Reyf_Mfibrosis, function(x){ensemblizer(x)}) )

BAL_macrophages <- add_named_module_score_ctrl_corrected(BAL_macrophages,
                                                         score_id = "Ayau_IPFeMp_lo_corr",
                                                         genesets_list_named = lapply(Ayau_IPFeM, function(x){ensemblizer(x)}) )

BAL_macrophages <- add_named_module_score_ctrl_corrected(BAL_macrophages,
                                                         score_id = "Mors_1_SPP_fib_lo_corr",
                                                         genesets_list_named = lapply(Mors_1_SPP_fib, function(x){ensemblizer(x)}) )

#extract scores from SEurat object
scores_reyf_clus1_loo_corr <- BAL_macrophages@meta.intermediate_data[, which(
  stringr::str_detect(string = names(BAL_macrophages@meta.intermediate_data), pattern = "Reyf_MCluster1_lo_corr")
)]

scores_loo_corr <- BAL_macrophages@meta.intermediate_data[, which(
  stringr::str_detect(string = names(BAL_macrophages@meta.intermediate_data), pattern = "lo_corr")
)]



saveRDS(scores_loo_corr,"intermediate_data/scores_loo_corr.rds")
}
}else("no calculation, loading files")

####calculate reduction of abolute scores, probably not so intersting because change of other clusters unkown




















# reyf_clus1_loo_corr<- donor_clus_avg_scores_reyf_clus1_loo_corr %>% filter(Cluster=="CD163LGMNM") %>% group_by(Geneset) %>% summarise(avg_score=mean(Score), stdev=sd(Score)) %>% arrange(avg_score)
# ref_ <- reyf_clus1_loo_corr %>% filter(str_detect(Geneset,"full")) %>% pull(avg_score)
# 
# reyf_clus1_loo_corr %<>% mutate(avg_score_diff=avg_score-ref_)
# reyf_clus1_loo_corr %>% mutate(Geneset=str_remove_all(Geneset,pattern="Reyf_Clus1_w_o_")) %>% ggplot(aes(fct_reorder(Geneset,avg_score_diff),avg_score_diff))+ geom_col()+theme(axis.text.x=element_text(angle = 90))  +coord_flip()
# reyf_clus1_loo_corr

scores_loo_corr <- read_rds(here::here("intermediate_data", "scores_loo_corr.rds"))

# #calculate mean score
# donor_clus_avg_scores_reyf_clus1_loo_corr <- scores_reyf_clus1_loo_corr %>%  mutate(
#             Cluster= as_factor(BAL_macrophages@meta.intermediate_data$Cluster),
#             donor=as_factor(BAL_macrophages@meta.intermediate_data$patient)) %>%
#     pivot_longer(cols = contains("Reyf_MCluster1_lo_corr"),names_to = "Geneset", values_to = "Score") %>% 
#   group_by(Cluster,donor, Geneset) %>% summarise(Score=mean(Score))%>% 
#   mutate(Cluster= factor(Cluster, levels=cluster_names))


donor_clus_avg_scores_loo_cor <- scores_loo_corr %>%  mutate(
  Cluster= as_factor(BAL_macrophages@meta.data$Cluster),
  donor=as_factor(BAL_macrophages@meta.data$patient)) %>%
  pivot_longer(cols = contains("lo_corr"),names_to = "gene_lo", values_to = "Score") 

# separate(col= gene_lo, into = c("gene_lo", "Geneset"),sep="_"#,extra = "merge" )
##this was to slow
gene_lo_slit <- donor_clus_avg_scores_loo_cor%>%  pull(gene_lo) %>% str_split_fixed(string=., pattern = "_", n=2 ) %>% as_tibble()

colnames(gene_lo_slit) <- c("gene_lo", "Geneset")
donor_clus_avg_scores_loo_cor <- donor_clus_avg_scores_loo_cor %>% 
  select(-gene_lo) %>% 
  bind_cols(gene_lo_slit) 

# group_by(Cluster,donor,gene_lo, Geneset) %>% summarise(Score=mean(Score))%>% 
#   mutate(Cluster= factor(Cluster, levels=cluster_names))




cohens_d_lo <- donor_clus_avg_scores_loo_cor %>%
  mutate(Cluster= fct_other(Cluster,keep = "CD163LGMNM")) %>% 
  group_by(gene_lo,Geneset) %>% 
  summarise(cohens=(cohen.d(Score[Cluster=="CD163LGMNM"],Score[Cluster!="CD163LGMNM"])$estimate), mean_dif=mean(Score[Cluster=="CD163LGMNM"])-mean(Score[Cluster!="CD163LGMNM"]) ) %>%
  ungroup() %>% 
  mutate(rank_conhens=rank(cohens), rank_meandiff=rank(mean_dif)) #%>% mutate(gene=str_remove_all(Geneset, "_Reyf_MCluster1_lo_corr") )

cohens_d_lo %>% arrange(gene_lo)
#write_rds(cohens_d_lo, "../output/cohens_d_lo.rds")

lo_corr_dot_plot_corr<- function(geneset, n_=16){
  x <- cohens_d_lo %>% filter(Geneset==geneset) %>% arrange(mean_dif) %>% slice_min(n = n_,order_by =mean_dif ) %>%  pull(gene_lo) 
  x <- x[!str_detect(x,"full")]
  y <- dot_plot_cluster(genes_sorted=x) +coord_flip()
  return(y)}
  



