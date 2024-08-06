library(effsize)

#compute averages of pre-calculated module scores over genesets, clusters and donors
donor_clus_avg_scores <- scores_fibrosis %>%  mutate(
  Cluster= as_factor(BAL_macrophages@meta.data$Cluster),
  donor=as_factor(BAL_macrophages@meta.data$patient)) %>%
  pivot_longer(cols = contains("score_fibrosis"),names_to = "Geneset", values_to = "Score") %>% 
  group_by(Cluster,donor, Geneset) %>% summarise(Score=mean(Score))%>% 
  mutate(Cluster= factor(Cluster, levels=cluster_names),
         Geneset=str_replace_all(Geneset, pattern= "_score_fibrosis", replacement = "") %>% factor())

p_avg_scores_wilcox <- wilcox_cluster_to_rest(data=donor_clus_avg_scores)

p_avg_scores_t <- t_test_cluster_to_rest(donor_clus_avg_scores)%>% arrange(p.value)

p_avg_scores_t %>% arrange(p.value)



#calculate cohnes d for each cluster average score of donors  in comparison to all other cluster-donor average scores
cohens_d_module_scores_fibrosis_donor <- tibble(geneset=character(),cluster=character(),cohens_d=numeric())
data <- donor_clus_avg_scores
for (gs in levels(data$Geneset)) {
  #print(gs)
  for (clus in levels(data$Cluster) ) { #data %>% filter(Geneset==gs) %>% pull(Cluster) %>% unique()
    #print(clus)
    ctrl <- data%>% filter(Geneset==gs, Cluster!=clus) %>% pull(Score)
    treat <- data%>% filter(Geneset==gs, Cluster==clus) %>% pull(Score)
    d <- cohen.d(treat, ctrl)
    cohens_d_module_scores_fibrosis_donor %<>%  add_row(geneset=gs,cluster=clus,cohens_d=d$estimate)
  }
}

two_highest <- donor_clus_avg_scores %>% group_by(Geneset, Cluster) %>% summarise(mean_score=mean(Score)) %>% top_n(wt=mean_score,2) %>% ungroup() %>%  mutate(rank="highest")%>% clean_names()
n_highest <- (two_highest$rank=="highest") %>% sum()



dat <- p_avg_scores_t %>% clean_names() %>% left_join(cohens_d_module_scores_fibrosis_donor,by = c("geneset", "cluster")) %>% arrange(-p_adjust) %>% left_join(two_highest,by = c("geneset", "cluster"))%>% 
  mutate(p_adjust_highest=ifelse(rank=="highest", p_value*n_highest, NULL)) %>% mutate(sig=if_else(p_adjust_highest<0.05,"sig","not sig.")) %>% arrange(p_adjust_highest) %>% mutate(cluster= factor(cluster, levels = names(color.cluster)))

dat <- dat%>% 
  filter(geneset!="Adam_ProfibroticM") %>% mutate(geneset=factor(geneset, levels =  c("Ayau_FABP4+Mp","Ayau_IPFeMp","Reyf_MCluster0" ,"Reyf_MCluster1" ,"Reyf_Mfibrosis", "Reyf_MCluster2", "Reyf_MCluster3","Mors_0FABP4+M" ,   "Mors_6FCN1+M"  ,"Mors_1SPP1+M","Mors_1SPP1+MIPF","overlap_set","intersect_set")   
  ), Cluster= factor(cluster, levels = cluster_names))

