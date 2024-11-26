

BAL_macrophages <- load_clean_bal_seurat()
gene_ensb <- BAL_macrophages@misc$features$SYMBOL
names(gene_ensb) <- BAL_macrophages@misc$features$ENSEMBL
features <- BAL_macrophages@misc$features
cluster_names <- c("Monocytes" , "MonoM" , "CD163LGMNM" ,    "AM1"  ,  "AM2"  ,   "ProlifAM")

color.cluster <- c(
  "Monocytes"     = "#009E73",
  "MonoM"       = "#E69F00",
  "CD163LGMNM" = "#D55E00",
  "AM1"         = "#A020F0",
  "AM2"         = "#0072B2",
  "ProlifAM"   = "#56B4E9"
)
BAL_macrophages@misc$colors$Cluster <- color.cluster

# Retrieve gene set dictionary



path_data_file <- "C:\\Users\\danne\\raw_data\\"
dictionary <- load_clean_dictionary(path_ = paste0(path_data_file,"wendisch_et_al_final\\oli_github\\genesets.xlsx"))
genesets_covid <- make_genesets(dictionary = dictionary, key = "COVID-19")
dict_covid <- make_dict(dictionary = dictionary, key = "COVID-19")

genesets_fibrosis<- make_genesets(dictionary = dictionary, key = "IPF")
dict_fibrosis <- make_dict(dictionary = dictionary, key="IPF")
###
overlap_set <- c("LILRB4","CTSZ", "APOE", "MMP9", "CCL2","CHIT1","CD14", "FNIP2", "SPP1", "LGMN","CD84","SDC2", "GPNMB")
intersect_set <- c( "SPP1", "LGMN","CD84","SDC2", "GPNMB")
genesets_fibrosis[["overlap_set"]] <- overlap_set
genesets_fibrosis[["intersect_set"]] <- intersect_set

#genesets_fibrosis$reyfman_fibrosis_cleaned <- genesets_fibrosis$Reyf_Mfibrosis[!genesets_fibrosis$Reyf_Mfibrosis %>% str_detect(pattern = "IG|RP11")]
BAL_macrophages$ClusternonASCII <- BAL_macrophages$Cluster
BAL_macrophages$Cluster <- BAL_macrophages$Cluster %>% iconv(from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% str_replace_all( "[[:punct:]]", "") %>% str_remove_all(" ")
Idents(BAL_macrophages) <- "Cluster"
BAL_macrophages@assays$SCoV2@key <- "scov2_"

BAL_macrophages <- ScaleData(BAL_macrophages)
########

############
BAL_macrophages <- add_named_module_score( 
  seurat_object   = BAL_macrophages, 
  genesets_list_named = genesets_fibrosis %>% lapply(function(x) ensemblizer(x)),
  score_id = "score_fibrosis"
)
scores_fibrosis <- BAL_macrophages@meta.data[, which(
  stringr::str_detect(string = names(BAL_macrophages@meta.data), pattern = "score_fibrosis")
)]

data_fibrosis <- extract_score_data(seurat_obj = BAL_macrophages, scores=scores_fibrosis, score_id = "score_fibrosis")%>% 
  mutate(Cluster= factor(Cluster, levels=cluster_names),
         Geneset=str_replace_all(Geneset, pattern= "_score_fibrosis", replacement = "") %>% factor())
ps <- wilcox_cluster_to_rest(data = data_fibrosis, geneset_var = "Geneset", cluster_var="Cluster")

hla_d<- BAL_macrophages@assays$RNA@data %>% 
  as_tibble(rownames = "ensg") %>% 
  mutate(gene=genifiy(ensg)) %>% relocate(gene) %>% 
  filter(str_detect(gene,"HLA-D")) %>% summarise(across(where(is.numeric), sum))
BAL_macrophages$HLA_d <- as.numeric(hla_d[1,])


color.cluster <- c(
  "Monocytes"     = "#009E73",
  "MonoM"       = "#E69F00",
  "CD163LGMNM" = "#D55E00",
  "AM1"         = "#A020F0",
  "AM2"         = "#0072B2",
  "ProlifAM"   = "#56B4E9"
)
BAL_macrophages@misc$colors$Cluster <- color.cluster
BAL_macrophages <- ProjectDim(BAL_macrophages)

clean_hgnc <- . %>% ensemblizer_2() %>% genifiy_2() %>% discard(is.na)
s1 <- genesets_fibrosis$Reyf_MCluster1 %>% clean_hgnc
s2 <- genesets_fibrosis$`Mors_1SPP1+MIPF`%>% clean_hgnc
s3 <- genesets_fibrosis$Ayau_IPFeMp%>% clean_hgnc
all_intersect <- Reduce(intersect, list(s1,s2,s3 ))
common_genes <- c(
  intersect(s1,s2),
  intersect(s2,s3),
  intersect(s1,s3)) %>% unique()
