
#################convenience to translate olis gene idententifiers
ensemblizer <- function(gene){
  BAL_macrophages@misc$features %>% 
    filter(SYMBOL%in%gene) %>%
    pull(ENSEMBL) %>% 
    return()
}

ensemblizer_2 <- function(gene){
  ens <- BAL_macrophages@misc$features$ENSEMBL
  names(ens) <- BAL_macrophages@misc$features$SYMBOL
    return(ens[gene])
}
#######################
genifiy<- function(gene){
  BAL_macrophages@misc$features %>% 
    filter(ENSEMBL%in%gene) %>%
    pull(SYMBOL) %>% 
    return()
}

genifiy_2<- function(gene){
  symb <- BAL_macrophages@misc$features$SYMBOL
  names(symb) <- BAL_macrophages@misc$features$ENSEMBL
  return(symb[gene])
}
load_clean_bal_seurat <- function(path_="C:/Users/danne/raw_data/wendisch_et_al_final/oli_github/BAL-macrophages.Rds"){
  BAL_macrophages <- readRDS(path_)
  #get rid of special characeters and greek letters that harm downstream functions
  BAL_macrophages$ClusternonASCII <- BAL_macrophages$Cluster
  BAL_macrophages$Cluster <- BAL_macrophages$Cluster %>% iconv(from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% str_replace_all( "[[:punct:]]", "") %>% str_remove_all(" ")
  
  BAL_macrophages@assays$SCoV2@key <- "scov2_"
  BAL_macrophages@assays$RNA@key <- "rna_"
  
  Idents(BAL_macrophages) <- "Cluster"
  levels(Idents(BAL_macrophages)) <- c("Monocytes" , "MonoM", "CD163LGMNM" ,"AM1",  "AM2", "ProlifAM" )
  return(BAL_macrophages)}
###############################
load_clean_dictionary <- function(path_="C:\\Users\\danne\\raw_data\\wendisch_et_al_final\\oli_github\\genesets.xlsx"){
  
  file <- path_
    dictionary <- readxl::read_excel(file)
  names(dictionary) <- c("ref", "term", "disease", "gene")
  dictionary$term <- dictionary$term %>% iconv(from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% str_replace_all( "[[:punct:]]", "") %>% str_remove_all(" ")
  dictionary <- dictionary %>% mutate(term=paste0(substr(ref,1,4),"_",term))
  return(dictionary)
  }
##############list of charcter vectors of genenames of given genesets
make_genesets <- function(key ="COVID-19", dictionary, n = 50){
  genesets <- list()
  for (i in unique(dictionary$term[dictionary$disease == key])) {
    genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
  }
  return(genesets)
}
############make list of dictionary dataframes
make_dict <- function(key, dictionary_, n = 50){
  dict <- list()
  for (i in unique(dictionary_$term[dictionary$disease == key])) {
    dict[[i]] <- head(dictionary_[dictionary$term == i, ], n)
  }
  return(dict)
}
########################
seurat_plot_renamer <- function(ggplot_, named_vec=gene_ensb){
  y_names<- ggplot_[[1]] %>% rownames()
  str_length_1 <-y_names[1] %>% nchar()
  str_length_2 <-y_names[length(y_names)] %>% nchar()
  
  y_groups <- y_names %>% substr(1,str_length_1)
  y_names_appendix <- y_names %>% substr(str_length_1+1,str_length_2)
  y_groups_new <- named_vec[y_groups]
  y_names_new <- paste0(y_groups_new,y_names_appendix)
  rownames(ggplot_[[1]]) <- y_names_new
  levels(ggplot_[[1]]$features.plot) <- named_vec[levels(ggplot_[[1]]$features.plot)]
  return(ggplot_)
}
###############
add_module_scores_custom <- function(){
  genesets <- genesets_list
  idsets <- list()
  # change to ensb id
  for (i in names(genesets)) {
    print(i)
    idsets[[i]] <- na.omit(
      features$ENSEMBL[match(genesets[[i]], features$SYMBOL)]
    )
  }
  #remove empty lsit elements
  idsets <- idsets[which(sapply(idsets, length) > 0)]
  
  #remove metadata columns that
  object@meta.data <- object@meta.data[, which(
    !stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
  )]
  
  object <- Seurat::AddModuleScore(
    object   = object, 
    features = idsets, 
    name     = "geneset",
    assay    = "RNA", 
    seed     = 1993
  )
  object
  
  return(object)
  
}
###########################
add_module_scores_custom_corrrected <- function(){
  genesets <- genesets_list
  idsets <- list()
  # change to ensb id
  for (i in names(genesets)) {
    print(i)
    idsets[[i]] <- na.omit(
      features$ENSEMBL[match(genesets[[i]], features$SYMBOL)]
    )
  }
  #remove empty lsit elements
  idsets <- idsets[which(sapply(idsets, length) > 0)]
  
  #remove metadata columns that
  object@meta.data <- object@meta.data[, which(
    !stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
  )]
  
  object <- AddModuleScore_ctr_corrected(
    object   = object, 
    features = idsets, 
    name     = "geneset",
    assay    = "RNA", 
    seed     = 1993
  )
  object
  
  return(object)
  
}
#################################


#####
###this is the functional one
add_named_module_score <- function(
    seurat_object, 
    assay="RNA",
    genesets_list_named,
    score_id,
    seed_=1993){
  
  seurat_object <- Seurat::AddModuleScore(
  object   = seurat_object, 
  features = genesets_list_named, 
  name     = "ertz1234",
  assay    = assay, 
  seed     = seed_
)

meta_feats <-  names(seurat_object@meta.data)
j <- 1
for (i in seq_along(meta_feats)){
  if(str_detect(meta_feats[i], pattern = "ertz1234")){
    meta_feats[i] <-  names(genesets_list_named)[j] %>% paste(score_id, sep="_")
    j <- j+1}
}
names(seurat_object@meta.data) <- meta_feats

return(seurat_object)
}
#########
add_named_module_score_ctrl_corrected <- function(
    seurat_object, 
    assay="RNA",
    genesets_list_named,
    score_id,
    seed_=1993){
  
  seurat_object <- AddModuleScore_ctr_corrected(
    object   = seurat_object, 
    features = genesets_list_named, 
    name     = "ertz1234",
    assay    = assay, 
    seed     = seed_
  )
  
  meta_feats <-  names(seurat_object@meta.data)
  j <- 1
  for (i in seq_along(meta_feats)){
    if(str_detect(meta_feats[i], pattern = "ertz1234")){
      meta_feats[i] <-  names(genesets_list_named)[j] %>% paste(score_id, sep="_")
      j <- j+1}
  }
  names(seurat_object@meta.data) <- meta_feats
  
  return(seurat_object)
}

#################
module_score_umaps <- function(object, 
                               #genesets_list, 
                               genesets_levels="",

                               limits = c(-0.5, 1),
                               pattern,
                               pt_size=0.25) {
  
  #copy metadata into separte dataframe
  scores <- object@meta.data[, which(
    stringr::str_detect(string = names(object@meta.data), pattern )
  )] #%>% as_tibble()
  #names(scores) <- names(idsets) ################
  ##############
  ############# still sonderzeichen
  
  scores$x <- object@reductions$umap@cell.embeddings[, 1]
  scores$y <- object@reductions$umap@cell.embeddings[, 2]
  
  plot_dat <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
  
  if (genesets_levels!=""){
    plot_dat$Geneset <- factor(
      plot_dat$Geneset, genesets_levels)}
  
  # Set colorscale limits
  
  plot_dat$Score[plot_dat$Score > max(limits)] <- max(limits)
  plot_dat$Score[plot_dat$Score < min(limits)] <- min(limits)
  ##############??
  plot_dat <- plot_dat[order(plot_dat$Score), ]
  
  # Plot
  p1 <- ggplot2::ggplot(
    data = plot_dat,
    mapping = ggplot2::aes(
      x   = x,
      y   = y,
      col = Score
    )
  ) +
    ggplot2::geom_point(size = 0.25) +
    ggplot2::facet_wrap(~Geneset, nrow = 4)+
    ggplot2::scale_color_distiller(
      palette = "RdYlBu", direction = -1, limits = limits
    ) +
    ggplot2::theme_void(base_size = 20) +
    ggplot2::guides(
      col = ggplot2::guide_colorbar(
        barwidth = 1, barheight = 15, 
        frame.colour = "black", ticks.colour = "black"
      )
    ) +
    ggplot2::coord_fixed() + 
    ggplot2::annotate(
      geom  = "segment", 
      size  = 0.25,
      x     = min(plot_dat$x)*1.1, 
      xend  = seq(min(plot_dat$x), max(plot_dat$x), length.out = 100)[25], 
      y     = min(plot_dat$y)*1.1, 
      yend  = min(plot_dat$y)*1.1, 
      arrow = grid::arrow(
        angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
      )
    ) +
    ggplot2::annotate(
      geom  = "segment",
      size  = 0.25,
      x     = min(plot_dat$x)*1.1, 
      xend  = min(plot_dat$x)*1.1, 
      y     = min(plot_dat$y)*1.1, 
      yend  = seq(min(plot_dat$y), max(plot_dat$y), length.out = 100)[25], 
      arrow = grid::arrow(
        angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
      )
    )
  return(p1)
}
#######

module_score_one_umap <- function(object, 
                               #genesets_list, 
                               genesets_levels="",
                               
                               limits = c(-0.5, 1),
                               pattern,
                               pt_size=0.25) {
  
  #copy metadata into separte dataframe
  scores <- tibble(
    object@meta.data[, which(
    stringr::str_detect(string = names(object@meta.data), pattern )
  )] )
  
  colnames(scores)[1] <- pattern
  #names(scores) <- names(idsets) ################
  ##############
  ############# still sonderzeichen
  
  scores$x <- object@reductions$umap@cell.embeddings[, 1]
  scores$y <- object@reductions$umap@cell.embeddings[, 2]
  
  plot_dat <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
  
  if (genesets_levels!=""){
    plot_dat$Geneset <- factor(
      plot_dat$Geneset, genesets_levels)}
  
  # Set colorscale limits
  
  plot_dat$Score[plot_dat$Score > max(limits)] <- max(limits)
  plot_dat$Score[plot_dat$Score < min(limits)] <- min(limits)
  ##############??
  plot_dat <- plot_dat[order(plot_dat$Score), ]
  
  # Plot
  p1 <- ggplot2::ggplot(
    data = plot_dat,
    mapping = ggplot2::aes(
      x   = x,
      y   = y,
      col = Score
    )
  ) +
    ggplot2::geom_point(size = pt_size) +
    ggplot2::facet_wrap(~Geneset, nrow = 4)+
    ggplot2::scale_color_distiller(
      palette = "RdYlBu", direction = -1, limits = limits
    ) +
    ggplot2::theme_void(base_size = 20) +
    ggplot2::guides(
      col = ggplot2::guide_colorbar(
        barwidth = 1, barheight = 15, 
        frame.colour = "black", ticks.colour = "black"
      )
    ) +
    ggplot2::coord_fixed() + 
    ggplot2::annotate(
      geom  = "segment", 
      size  = 0.25,
      x     = min(plot_dat$x)*1.1, 
      xend  = seq(min(plot_dat$x), max(plot_dat$x), length.out = 100)[25], 
      y     = min(plot_dat$y)*1.1, 
      yend  = min(plot_dat$y)*1.1, 
      arrow = grid::arrow(
        angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
      )
    ) +
    ggplot2::annotate(
      geom  = "segment",
      size  = 0.25,
      x     = min(plot_dat$x)*1.1, 
      xend  = min(plot_dat$x)*1.1, 
      y     = min(plot_dat$y)*1.1, 
      yend  = seq(min(plot_dat$y), max(plot_dat$y), length.out = 100)[25], 
      arrow = grid::arrow(
        angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
      )
    )
  return(p1)
}



#####################
score_violin_theme <- function(limits, p_score){list(
  ggplot2::geom_point(
    position = "jitter"
  ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) ,
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) ,
  ggplot2::facet_wrap(~Geneset, nrow = 4) ,
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) ,
  ggplot2::scale_fill_manual(values = color.cluster) ,
  ggplot2::theme_classic(base_size = 15) ,
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) , 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) ,
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      #barwidth = 8, barheight = 1,
      frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) ,
  ggplot2::labs(col = NULL, x = NULL) ,
  ggnewscale::new_scale("fill") ,
  ggplot2::geom_label(
    data = p_score,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) ,
  viridis::scale_fill_viridis(direction = -1, option = "A")
)}
#######################

#####################
score_violin_theme_patchwork <- function(limits, p_score){list(
  ggplot2::geom_point(
    position = "jitter",size=0.1
  ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 1
  ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) ,
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) ,
  #ggplot2::facet_wrap(~Geneset, nrow = 4) ,
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", 
    direction = -1, limits = limits
  ) ,
  ggplot2::scale_fill_manual(values = color.cluster) ,
  ggplot2::theme_classic(#base_size = 1
                         ) ,
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) , 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) ,
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      #barwidth = 8, barheight = 1,
      frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) ,
  ggplot2::labs(col = NULL, x = NULL) ,
  ggnewscale::new_scale("fill") ,
  ggplot2::geom_label(
    data = p_score,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 2,
    color = "white"
  ) ,
  viridis::scale_fill_viridis(direction = -1, option = "A")
)}
#######################
score_violin_theme_patchwork_no_points <- function(limits, p_score){list(
  # ggplot2::geom_point(
  #   position = "jitter",size=0.1
  # ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 1
  ) ,
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) ,
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) ,
  #ggplot2::facet_wrap(~Geneset, nrow = 4) ,
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", 
    direction = -1, limits = limits
  ) ,
  ggplot2::scale_fill_manual(values = color.cluster) ,
  ggplot2::theme_classic(#base_size = 1
  ) ,
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) , 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) ,
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      #barwidth = 8, barheight = 1,
      frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) ,
  ggplot2::labs(col = NULL, x = NULL) ,
  ggnewscale::new_scale("fill") ,
  ggplot2::geom_label(
    data = p_score,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 2,
    color = "white"
  ) ,
  viridis::scale_fill_viridis(direction = -1, option = "A")
)}
#################


wilcox_cluster_to_rest <- function(data, geneset_var="Geneset", cluster_var="Cluster"){
  # Calculate p-values
  ps <- unique(data[, c(geneset_var,cluster_var)])
  ps$p.value <- NA
  # DW calculate p_values of Cluster versus rest for all clusters
  for (i in levels(data[[geneset_var]])) {
    #print(i)
    for (ii in levels(data[[cluster_var]])) {
      #print(ii)
      ps$p.value[ps[[geneset_var]] == i & ps[[cluster_var]] == ii] <- wilcox.test(
        x = data$Score[data[[geneset_var]] == i & data[[cluster_var]] == ii],
        y = data$Score[data[[geneset_var]] == i], 
        alternative = "greater"
      )$p.value
    }
  }
  ps$p.adjust <- p.adjust(ps$p.value)
  ps$label <- round(-log10(ps$p.adjust))
  
  ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])
  ps %>% arrange(p.value)
  return(ps)
}
###############
t_test_cluster_to_rest <- function(data, geneset_var="Geneset", cluster_var="Cluster"){
  # Calculate p-values
  ps <- unique(data[, c(geneset_var,cluster_var)])
  ps$p.value <- NA
  # DW calculate p_values of Cluster versus rest for all clusters
  for (i in levels(data[[geneset_var]])) {
    #print(i)
    for (ii in levels(data[[cluster_var]])) {
      #print(ii)
      ps$p.value[ps[[geneset_var]] == i & ps[[cluster_var]] == ii] <- t.test(
        x = data$Score[data[[geneset_var]] == i & data[[cluster_var]] == ii],
        y = data$Score[data[[geneset_var]] == i], 
        alternative = "greater"
      )$p.value
    }
  }
  ps$p.adjust <- p.adjust(ps$p.value)
  ps$label <- round(-log10(ps$p.adjust))
  
  ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])
  ps %>% arrange(p.value)
  return(ps)
}
############
cohens_d_cluster_to_rest <- function(data, geneset_var="Geneset", cluster_var="Cluster"){
  # Calculate p-values
  ps <- unique(data[, c(geneset_var,cluster_var)])
  ps$cohens_d <- NA
  # DW calculate p_values of Cluster versus rest for all clusters
  for (i in levels(data[[geneset_var]])) {
    #print(i)
    for (ii in levels(data[[cluster_var]])) {
      #print(ii)
      ps$cohens_d[ps[[geneset_var]] == i & ps[[cluster_var]] == ii] <- wilcox.test(
        x = data$Score[data[[geneset_var]] == i & data[[cluster_var]] == ii],
        y = data$Score[data[[geneset_var]] == i], 
        alternative = "greater"
      )$cohens_d
    }
  }
  ps$p.adjust <- p.adjust(ps$p.value)
  ps$label <- round(-log10(ps$p.adjust))
  
  ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])
  ps %>% arrange(p.value)
  return(ps)
}


#####################
extract_score_data <-  function(seurat_obj=BAL_macrophages, scores_, limits=c(-0.5,1), score_id) { scores_ %>%
    mutate( x= seurat_obj@reductions$umap@cell.embeddings[, 1],
            y=seurat_obj@reductions$umap@cell.embeddings[, 2],
            Cluster= as_factor(seurat_obj@meta.data$Cluster)) %>% 
    pivot_longer(cols = contains(score_id),names_to = "Geneset", values_to = "Score") %>%
    mutate(Geneset=as_factor(Geneset)) %>% 
    group_by(Geneset) %>% 
    mutate(ms=mean(Score)) %>% 
    ungroup() %>% 
    mutate(col=if_else(Score>max(limits), 
                       max(limits),
                       if_else(Score<min(limits),
                               min(limits),
                               Score
                       )
    )
    ) 
  
}
###############
# extract_score_data <-  function(seurat_obj=BAL_macrophages, scores_, limits=c(-0.5,1), score_id) { scores_ %>%
#     mutate( x= seurat_obj@reductions$umap@cell.embeddings[, 1],
#             y=seurat_obj@reductions$umap@cell.embeddings[, 2],
#             Cluster= as_factor(seurat_obj@meta.data$Cluster),
#             donor=as_factor(seurat_obj@meta.data$patient)) %>% 
#     pivot_longer(cols = contains(score_id),names_to = "Geneset", values_to = "Score") %>%
#     mutate(Geneset=as_factor(Geneset)) %>% 
#     group_by(Geneset) %>% 
#     mutate(ms=mean(Score)) %>% 
#     ungroup() %>% 
#     mutate(col=if_else(Score>max(limits), 
#                        max(limits),
#                        if_else(Score<min(limits),
#                                min(limits),
#                                Score
#                        )
#     )
#     ) 
#   
# }
################
pca_loadings_plot <- function(object=BAL_macrophages, feature, dims_ = c(1:2), pc="PC_1"){
  p1 <- FeaturePlot(object, reduction = "pca", features=ensemblizer(feature), dims_)+ scale_color_viridis_c() + xlim(-0.5,0.3) + ylim(-0.5,0.25) +ggtitle(feature)
  colnames(p1$data)[4] <- feature
  p2 <- p1$data %>% ggplot(aes_string(pc, feature )) + geom_point(aes(color=ident),alpha=0.3) +geom_smooth(color="black", size=2) + ggtitle(feature) + theme(legend.position = "none")
  return(p2+theme_1)
}
pca_loadings_plot_clean <- function(object=BAL_macrophages, feature, dims_ = c(1:2), pc="PC_1"){#same as pca_loadings_plot but with name gene name conversion as error handling with genes with "-"
  p1 <- FeaturePlot(object, reduction = "pca", features=ensemblizer(feature), dims_)+ scale_color_viridis_c() + xlim(-0.5,0.3) + ylim(-0.5,0.25) +ggtitle(feature)
  colnames(p1$data)[4] <- feature
  p2 <- p1$data %>% ggplot(aes_string(pc, janitor::make_clean_names(feature) )) + geom_point(aes(color=ident),alpha=0.3) +geom_smooth(color="black", size=2) + ggtitle(feature) + theme(legend.position = "none")
  return(p2+theme_1)
}

pca_loadings_plot_90 <- function(object=BAL_macrophages, feature, dims_ = c(1:2), pc="PC_1"){
  p1 <- FeaturePlot(object, reduction = "pca", features=ensemblizer(feature), dims_)+ scale_color_viridis_c() + xlim(-0.5,0.3) + ylim(-0.5,0.25) +ggtitle(feature)
  colnames(p1$data)[4] <- feature
  p1$data[[feature]] <- p1$data[[feature]]*-1
  p2 <- p1$data %>% mutate()%>% ggplot(aes_string( pc,feature )) + geom_point(aes(color=ident),alpha=0.3) +geom_smooth(color="black", size=2) + ggtitle(feature) + theme(legend.position = "none")#+coord_flip()
  return(p2+coord_flip()+ theme_2)
}

pca_loadings_plot_90_titled <- function(object=BAL_macrophages, feature, dims_ = c(1:2), pc="PC_1"){
  p1 <- FeaturePlot(object, reduction = "pca", features=ensemblizer(feature), dims_)+ scale_color_viridis_c() + xlim(-0.5,0.3) + ylim(-0.5,0.25) +ggtitle(feature)
  colnames(p1$data)[4] <- feature
  p1$data[[feature]] <- p1$data[[feature]]*-1
  p2 <- p1$data %>% mutate()%>% ggplot(aes_string( pc,feature )) + geom_point(aes(color=ident),alpha=0.3) +geom_smooth(color="black", size=2) + ggtitle(feature) + theme(legend.position = "none")#+coord_flip()
  return(p2+coord_flip()+ theme_2_titeled)
}
##############
pca_loadings_plot_meta_data <- function(object=BAL_macrophages, feature, dims_ = c(1:2), pc="PC_1"){
  p1 <- FeaturePlot(object, reduction = "pca", features=feature, dims_)+
    scale_color_viridis_c() + xlim(-0.5,0.3) + ylim(-0.5,0.25) +ggtitle(feature)
  colnames(p1$data)[4] <- feature
  p2 <- p1$data %>% ggplot(aes_string(pc, feature )) + geom_point(aes(color=ident),alpha=0.3) +
    geom_smooth(color="black", size=2) + ggtitle(feature) + theme(legend.position = "none")
  return(p2+theme_1)
}
##########################
theme_1 <-   theme(text = element_text(size=6),
                   axis.text = element_text(size=6),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   plot.title = element_blank(),legend.background = element_blank(), 
                   legend.key = element_blank(), panel.background =  element_rect(fill = "white", colour = "white"),
                   panel.grid.major = element_line(colour = "grey"),
                   panel.border = element_blank(), strip.background = element_blank(), 
                   plot.background = element_blank() ,
                   plot.margin = unit(c(0, 0, 0, 0), "cm")#removes spacing around plots
)

theme_2 <- theme(text = element_text(size=6),
                 axis.text = element_text(size=6),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 plot.title = element_blank(),
                 legend.background = element_blank(), 
                 legend.key = element_blank(), 
                 panel.background =  element_rect(fill = "white", colour = "white"),
                 panel.grid.major = element_line(colour = "grey"),
                 panel.border = element_blank(), strip.background = element_blank(), 
                 plot.background = element_blank(), 
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank(),
                 plot.margin = unit(c(0, 0, 0, 0), "cm")#removes spacing around plots
)

theme_2_titeled <- theme(text = element_text(size=6),
                         axis.text = element_text(size=6),
                         axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 #plot.title = element_blank(),
                 legend.background = element_blank(), 
                 legend.key = element_blank(), 
                 panel.background =  element_rect(fill = "white", colour = "white"),
                 panel.grid.major = element_line(colour = "grey"),
                 panel.border = element_blank(), strip.background = element_blank(), 
                 plot.background = element_blank(), 
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank(),
                 plot.margin = unit(c(0, 0, 0, 0), "cm")#removes spacing around plots
)
##################


my_hclust_BuildCLusterTree <- function (object, assay = NULL, features = NULL, dims = NULL, 
          reduction = "pca", graph = NULL, slot = "data", 
          reorder = FALSE, reorder.numeric = FALSE, verbose = TRUE) 
{
  if (!PackageCheck("ape", error = FALSE)) {
    stop(cluster.ape, call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = graph)) {
    idents <- levels(x = object)
    nclusters <- length(x = idents)
    data.dist <- matrix(data = numeric(length = 1L), nrow = nclusters, 
                        ncol = nclusters, dimnames = list(idents, idents))
    graph <- object[[graph]]
    cxi <- CellsByIdentities(object = object)
    cpairs <- na.omit(object = unique(x = t(x = apply(X = expand.grid(1:nclusters, 
                                                                      1:nclusters)[, c(2, 1)], MARGIN = 1, FUN = function(x) {
                                                                        if (length(x = x) == length(x = unique(x = x))) {
                                                                          return(sort(x = x))
                                                                        }
                                                                        return(c(NA, NA))
                                                                      }))))
    if (verbose) {
      pb <- txtProgressBar(style = 3, file = stderr())
    }
    for (i in 1:nrow(x = cpairs)) {
      i1 <- cpairs[i, ][1]
      i2 <- cpairs[i, ][2]
      graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
      d <- mean(x = graph.sub)
      if (is.na(x = d)) {
        d <- 0
      }
      data.dist[i1, i2] <- d
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i/nrow(x = cpairs))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(x = data.dist)
  }
  else if (!is.null(x = dims)) {
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    embeddings <- Embeddings(object = object, reduction = reduction)[, 
                                                                     dims]
    data.dims <- my.lapply(X = levels(x = object), FUN = function(x) {
      cells <- WhichCells(object = object, idents = x)
      if (length(x = cells) == 1) {
        cells <- c(cells, cells)
      }
      temp <- colMeans(x = embeddings[cells, ])
    })
    data.dims <- do.call(what = "cbind", args = data.dims)
    colnames(x = data.dims) <- levels(x = object)
    data.dist <- dist(x = t(x = data.dims))
  }
  else {
    features <- features %||% VariableFeatures(object = object)
    features <- intersect(x = features, y = rownames(x = object))
    data.avg <- AverageExpression(object = object, assays = assay, 
                                  features = features, slot = slot, verbose = verbose)[[1]]
    data.dist <- dist(x = t(x = data.avg[features, ]))
  }
  x = hclust(d = data.dist)

  return(x)
}
###################
AddModuleScore_ctr_corrected <- function (object, features, pool = NULL, nbin = 24, ctrl = 100, 
          k = FALSE, assay = NULL, name = "Cluster", seed = 1, 
          search = FALSE, ...) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                         i))
    }
    cluster.length <- length(x = features)
  }
  else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning("The following features are not present in the object: ", 
                paste(missing.features, collapse = ", "), 
                ifelse(test = search, yes = ", attempting to find updated synonyms", 
                       no = ", not searching for symbol synonyms"), 
                call. = FALSE, immediate. = TRUE)
        if (search) {
          tryCatch(expr = {
            updated.features <- UpdateSymbolList(symbols = missing.features, 
                                                 ...)
            names(x = updated.features) <- missing.features
            for (miss in names(x = updated.features)) {
              index <- which(x == miss)
              x[index] <- updated.features[miss]
            }
          }, error = function(...) {
            warning("Could not reach HGNC's gene names database", 
                    call. = FALSE, immediate. = TRUE)
          })
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning("The following features are still not present in the object: ", 
                    paste(missing.features, collapse = ", "), 
                    call. = FALSE, immediate. = TRUE)
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  # if (!all(LengthCheck(values = features))) {
  #   warning(paste("Could not find enough features in the object from the following feature lists:", 
  #                 paste(names(x = which(x = !LengthCheck(values = features)))), 
  #                 "Attempting to match case..."))
  #   features <- lapply(X = features.old, FUN = CaseMatch, 
  #                      match = rownames(x = object))
  # }
  
  # if (!all(LengthCheck(values = features))) {
  #   stop(paste("The following feature lists do not have enough features present in the object:", 
  #              paste(names(x = which(x = !LengthCheck(values = features)))), 
  #              "exiting..."))
  # }
  print("length check distabled")
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ### ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
    ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}
##############
#dotplot

dot_plot_cluster <- function(genes_sorted){
  ids_1 <- ensemblizer_2(genes_sorted)
  ids <- ids_1[!is.na(ids_1)]
  if(length(setdiff(ids_1,ids))>0){
  paste("gene names....", paste(genes_sorted[is.na(ids_1)]), "....converted to NA and removed from geneset") %>%
    print()}
data <- tidyr::as_tibble(t(as.matrix(BAL_macrophages@assays$RNA@data[ids, ])))
colnames(data) <- names(ids)

# Add metadata
data$Cluster <- factor(BAL_macrophages$Cluster, levels = cluster_names)

# Tidy data & compute summaries
data <- tidyr::gather(data, "Gene", "Expression", -Cluster)
data$Gene <- factor(data$Gene, levels = genes_sorted)
data$N <- rep(1, length(data$Cluster))
data$Expressed <- data$Expression > 0
data <- dplyr::group_by(data, Cluster, Gene)
data <- dplyr::summarize(
  .data        = data,
  "Expression" = mean(Expression),
  "Cells"      = sum(N),
  "Expressed"  = sum(Expressed)
)
data$PCT <- data$Expressed / data$Cells * 100

# Scale data & set limits
data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(
  data, zscore = scale(Expression)
)
limits <- c(-2, 2)
data$zscore[data$zscore > max(limits)] <- max(limits)
data$zscore[data$zscore < min(limits)] <- min(limits)

# Plot
ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x    = Gene,
    y    = Cluster,
    col  = zscore,
    size = PCT
  )
) +
  ggplot2::geom_point() +
  ###
  ggplot2::  scale_size_continuous(range = c(0, 3))+
  
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
     # barwidth  = 18, barheight = 1,
      ticks = FALSE, frame.colour = "black",
      title.position = "top", title.hjust = 0.5
    ),
    size  = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  ) +
  ggplot2::labs(
    col  = "z-score", size = "% expressed"
  ) +
   ggplot2::theme_classic(#base_size = 25 +
   )+
  ggplot2::theme(text = element_text(size = 6),
                 axis.text = element_text(size = 6),
    legend.position = "top",
    axis.text.x     = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, face = "italic"
    ),
    axis.title      = ggplot2::element_blank()
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdBu", limits = limits
  ) #+coord_flip()
}
######gene_tf_cor_heatmap

gene_tf_cor_heatmap <- function(genes=c("CD14","S100A8",  "APOC1", "HLA-DRB1", "HLA-DQA1", "C1QA", "LGMN", "CD163"),
                                gene_matrix=BAL_macrophages@assays$RNA@data){
  library(readxl)
  library(janitor)
  tf <- read_excel("data/DatabaseExtract_v_1.01.xlsx")
  tf_genes <- tf %>% clean_names() %>%  filter(is_tf=="Yes") %>% pull(hgnc_symbol) %>% unique()
  #matrix_mod<-BAL_macrophages_small_1@assays$RNA@data %>% as.matrix()
  rele_genes <- c(tf_genes,genes) %>% ensemblizer_2()
  rele_genes <- rele_genes[!is.na(rele_genes)]
  matrix_mod<-gene_matrix[ rownames(gene_matrix) %in% rele_genes ,]  %>% as.matrix()
  #matrix_mod <- matrix_mod[rowSums(matrix_mod)>5,]
  rownames(matrix_mod) <- rownames(matrix_mod) %>% genifiy_2()
  cor_list <- list()
  
  for (gene_name in genes) {
    #gene_name <- "CSF1R"
    gene<-as.numeric(matrix_mod[gene_name,])
    correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
    correlations_tbl <- tibble(gene=names(correlations), cor=correlations, cor_gene=gene_name)
    
    cor_list[[gene_name]] <- correlations_tbl
    
  }
  
  cor_tbl <- cor_list %>% purrr::reduce(bind_rows) #%>% %>% mutate(cor_gene=fct_reorder(cor_gene,.x = genes))
  
  cor_tbl <- cor_list %>% purrr::reduce(bind_rows) %>%  mutate(cor_gene=fct_relevel(cor_gene,genes))
  #write_rds(cor_tbl,"cor_tbl.rds")
  
  top_cor_genes <- cor_tbl %>% filter(gene %in% tf_genes) %>% group_by(cor_gene) %>% slice_max(cor,n=3) %>% pull(gene) %>% unique()
  
  cor_tbl <- cor_tbl%>%filter(gene%in%top_cor_genes)%>% group_by(cor_gene) %>% arrange(cor,.by_group = TRUE) %>% mutate(gene=fct_inorder(gene))
  #write_rds(cor_tbl,"cor_tbl.rds")
  cor_hp <- cor_tbl %>%filter(gene%in%top_cor_genes) %>%
    group_by(cor_gene)%>% ungroup() %>% 
    ggplot(aes(gene,cor_gene, fill=cor)) +geom_tile()+ theme(axis.text.x=element_text(angle = 90))+scale_fill_distiller(
      palette = "RdBu"
    )
  return(cor_hp)

}

#############


