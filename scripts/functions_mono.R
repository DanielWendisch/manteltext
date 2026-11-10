ensemblizer_seur <- function(gene, seurat_obj=mono){
  ens <- seurat_obj@misc$features$ENSEMBL
  names(ens) <- seurat_obj@misc$features$SYMBOL
  return(ens[gene])
}


genifiy_seur<- function(gene, seurat_obj=mono){
  symb <- seurat_obj@misc$features$SYMBOL
  names(symb) <- seurat_obj@misc$features$ENSEMBL
  return(symb[gene])
}


add_named_module_score <- 
function(
    seurat_object, 
    assay="RNA",
    genesets_list_named,
    score_id,
    seed_=1993,
    nbin_=25){
  
  seurat_object <- Seurat::AddModuleScore(
    object   = seurat_object, 
    features = genesets_list_named, 
    name     = "ertz1234",
    assay    = assay, 
    seed     = seed_,
    nbin = nbin_
      
    
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

########
dot_plot_condition_seur <- function(genes_sorted, seur_obj, grouping_){
  ids_1 <- ensemblizer_seur(genes_sorted)
  ids <- ids_1[!is.na(ids_1)]
  if(length(setdiff(ids_1,ids))>0){
    paste("gene names....", paste(genes_sorted[is.na(ids_1)]), "....converted to NA and removed from geneset") %>%
      print()}
  data <- tidyr::as_tibble(t(as.matrix(seur_obj@assays$RNA@data[ids, ])))
  colnames(data) <- names(ids)
  
  # Add metadata
  data$condition <- factor(seur_obj$condition#, levels = condition_names
                       )
  
  # Tidy data & compute summaries
  data <- tidyr::gather(data, "Gene", "Expression", -condition)
  data$Gene <- factor(data$Gene, levels = genes_sorted)
  data$N <- rep(1, length(data$condition))
  data$Expressed <- data$Expression > 0
  data <- dplyr::group_by(data, condition, Gene)
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
      y    = condition,
      col  = zscore,
      size = PCT
    )
  ) +
    ggplot2::geom_point() +
    ###
    ggplot2::  scale_size_continuous(range = c(0, 3))+
    
    # ggplot2::guides(
    #   color = ggplot2::guide_colorbar(
    #     # barwidth  = 18, barheight = 1,
    #     ticks = FALSE, frame.colour = "black",
    #     title.position = "top", title.hjust = 0.5
    #   ),
    #   size  = ggplot2::guide_legend(
    #     title.position = "top", title.hjust = 0.5
    #   )
    # ) +
    # ggplot2::labs(
    #   col  = "z-score", size = "% expressed"
    # ) +
    # ggplot2::theme_classic(#base_size = 25 +
    # )+
    # ggplot2::theme(text = element_text(size = 6),
    #                axis.text = element_text(size = 6),
    #                legend.position = "top",
    #                axis.text.x     = ggplot2::element_text(
    #                  angle = 45, hjust = 1, vjust = 1, face = "italic"
    #                ),
    #                axis.title      = ggplot2::element_blank()
    # ) +
    # ggplot2::scale_color_distiller(
    #   palette = "RdBu", limits = limits
    # ) 
  coord_flip()
}

#####intitalize
dictionary <- load_clean_dictionary(path_ = here::here("..", "..", "raw_data", "wendisch_et_al_final", "oli_github", "genesets.xlsx"))
genesets_covid <- make_genesets(dictionary = dictionary, key = "COVID-19")
dict_covid <- make_dict(dictionary = dictionary, key = "COVID-19")

genesets_fibrosis<- make_genesets(dictionary = dictionary, key = "IPF")
dict_fibrosis <- make_dict(dictionary = dictionary, key="IPF")

overlap_set <- c("LILRB4","CTSZ", "APOE", "MMP9", "CCL2","CHIT1","CD14", "FNIP2", "SPP1", "LGMN","CD84","SDC2", "GPNMB")
intersect_set <- c( "SPP1", "LGMN","CD84","SDC2", "GPNMB")
genesets_fibrosis[["overlap_set"]] <- overlap_set
genesets_fibrosis[["intersect_set"]] <- intersect_set




