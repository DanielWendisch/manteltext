








 c("NRP1") %in% (BAL_macrophages@assays$RNA@var.features %>% genifiy())








 # .downsampled seurat oject

# ..Module scores downsampled 
## ..1. UMAP embedding with cell-based scores
```{r}

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_COVID-19-geneset-score-umap_small_1.png", width = 10, height = 6
)

```
## ..2. Violin plots as cluster-based population summaries

```{r}
# Re-shape data
data <- scores
data$Cluster <- object_small_1@meta.data$Cluster
data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(11,1,8,12,5,6,2,13,14,7,3,15,4,16,9,10)]
)
data$col <- data$Score

# Set colorscale limits
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)

data <- dplyr::mutate(dplyr::group_by(data, Geneset), ms = mean(Score))
data %>% head()
```

```{r}
i <- levels(data$Geneset)[1]
ii <- unique(data$Cluster)[2]

data[data$Geneset == i & data$Cluster == ii,] %>% head()
```


```{r}
# Calculate p-values
ps <- unique(data[, c("Geneset", "Cluster")])
ps$p.value <- NA
# DW calculate pairwise 
for (i in levels(data$Geneset)) {
  print(i)
  for (ii in levels(data$Cluster)) {
    print(ii)
    ps$p.value[ps$Geneset == i & ps$Cluster == ii] <- wilcox.test(
      x = data$Score[data$Geneset == i & data$Cluster == ii],
      y = data$Score[data$Geneset == i], 
      alternative = "greater"
    )$p.value
  }
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))

ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])
ps %>% arrange(p.value)
```


```{r}
# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter"
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) +
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) +
  ggplot2::facet_wrap(~Geneset, nrow = 4) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::scale_fill_manual(values = color.cluster) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 8, barheight = 1, frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_COVID-19-geneset-score-violins_small_1.png", width = 12, height = 8
)
```

```{r}
# Select disease and gene set size
key <- "IPF"
n <- 50

dict <- list()
genesets <- list()
for (i in unique(dictionary$term[dictionary$disease == key])) {
  dict[[i]] <- head(dictionary[dictionary$term == i, ], n)
  genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
}
dict <- dplyr::bind_rows(dict)
```


# ..Module scores 
## ..1. UMAP embedding with cell-based scores


```{r}
idsets <- list()
for (i in names(genesets)) {
  print(i)
  idsets[[i]] <- na.omit(
    features$ENSEMBL[match(genesets[[i]], features$SYMBOL)]
  )
}
idsets <- idsets[which(sapply(idsets, length) > 0)]

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

scores <- object@meta.data[, which(
  stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
)]
names(scores) <- names(idsets)

scores$x <- object@reductions$umap@cell.embeddings[, 1]
scores$y <- object@reductions$umap@cell.embeddings[, 2]

data <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
data$Geneset <- factor(data$Geneset, levels = unique(data$Geneset))
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
)

# Set colorscale limits
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)

data <- data[order(data$Score), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(size = 0.25) +
  ggplot2::facet_wrap(~Geneset, nrow = 3)+
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
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  )

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_IPF-geneset-score-umap.png", width = 12, height = 6
)
```
## ..2. Violin plots as cluster-based population summaries

```{r}


data <- scores
data$Cluster <- object@meta.data$Cluster

data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)

data$Geneset <- factor(data$Geneset, unique(data$Geneset))
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
)

data$col <- data$Score
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)

data <- dplyr::mutate(dplyr::group_by(data, Geneset), ms = mean(Score))

# Calculate p-values
ps <- unique(data[, c("Geneset", "Cluster")])
ps$p.value <- NA
for (i in levels(data$Geneset)) {
  print(i)
  for (ii in levels(data$Cluster)) {
    print(ii)
    ps$p.value[ps$Geneset == i & ps$Cluster == ii] <- wilcox.test(
      x = data$Score[data$Geneset == i & data$Cluster == ii],
      y = data$Score[data$Geneset == i], 
      alternative = "greater"
    )$p.value
  }
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))

ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter"
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) +
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) +
  ggplot2::facet_wrap(~Geneset, nrow = 3) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::scale_fill_manual(values = color.cluster) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 8, barheight = 1, frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_IPF-geneset-score-violins.png", width = 12, height = 6
)

# ------------------------------------------------------------------------------
# Save dataset to disk

#saveRDS(object, "BAL-macrophages_small_1.Rds")

# end of document
################################################################################
```








```{r}
# Select disease and gene set size
key <- "IPF"
n <- 50

dict <- list()
genesets <- list()
for (i in unique(dictionary$term[dictionary$disease == key])) {
  dict[[i]] <- head(dictionary[dictionary$term == i, ], n)
  genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
}
dict <- dplyr::bind_rows(dict)


```
## downsampled Seurat object



```{r}
#DW I guess this is the genes object that oli wanted to use for the code beneath
cells <- colnames(object)[which(
  object@meta.data$Celltype == "Macrophages"
)]
genes <- list()
for (assay in names(object@assays)) {
  print(assay)
  genes[[assay]] <- names(which(
    Matrix::rowSums(
      slot(subset(object, cells = cells)@assays[[assay]], "counts")
    ) > 0
  ))
}
genes <- as.character(unlist(genes))
```





## ..Module scores 
## ..1. UMAP embedding with cell-based scores


```{r}
idsets <- list()
for (i in names(genesets)) {
  print(i)
  idsets[[i]] <- na.omit(
    features$ENSEMBL[match(genesets[[i]], features$SYMBOL)]
  )
}
idsets <- idsets[which(sapply(idsets, length) > 0)]

object_small_1@meta.data <- object_small_1@meta.data[, which(
  !stringr::str_detect(string = names(object_small_1@meta.data), pattern = "geneset")
)]

object_small_1 <- Seurat::AddModuleScore(
  object   = object_small_1, 
  features = idsets, 
  name     = "geneset",
  assay    = "RNA", 
  seed     = 1993
)

scores <- object_small_1@meta.data[, which(
  stringr::str_detect(string = names(object_small_1@meta.data), pattern = "geneset")
)]
names(scores) <- names(idsets)

scores$x <- object_small_1@reductions$umap@cell.embeddings[, 1]
scores$y <- object_small_1@reductions$umap@cell.embeddings[, 2]

data <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
data$Geneset <- factor(data$Geneset, levels = unique(data$Geneset))
# data$Geneset <- factor(
#   data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
# )

# Set colorscale limits
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)

data <- data[order(data$Score), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(size = 0.25) +
  ggplot2::facet_wrap(~Geneset, nrow = 3)+
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
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  )

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_IPF-geneset-score-umap_small_1.png", width = 12, height = 6
)
```
## ..2. Violin plots as cluster-based population summaries

```{r}
data <- scores
data$Cluster <- object_small_1@meta.data$Cluster

data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)

data$Geneset <- factor(data$Geneset, unique(data$Geneset))
# data$Geneset <- factor(
#   data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
# )

data$col <- data$Score
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)

data <- dplyr::mutate(dplyr::group_by(data, Geneset), ms = mean(Score))
data$Cluster %<>% factor(c("Monocytes" , "MonoM", "CD163LGMNM" ,"AM1",  "AM2", "ProlifAM" ))

# Calculate p-values
ps <- unique(data[, c("Geneset", "Cluster")])
ps$p.value <- NA
for (i in levels(data$Geneset)) {
  print(i)
  for (ii in levels(data$Cluster)) {
    print(ii)
    ps$p.value[ps$Geneset == i & ps$Cluster == ii] <- wilcox.test(
      x = data$Score[data$Geneset == i & data$Cluster == ii],
      y = data$Score[data$Geneset == i], 
      alternative = "greater"
    )$p.value
  }
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))

ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter"
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) +
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) +
  ggplot2::facet_wrap(~Geneset, nrow = 3) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::scale_fill_manual(values = color.cluster) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 8, barheight = 1, frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# Save plot
ggplot2::ggsave(
  "BAL-macrophages_IPF-geneset-score-violins_small_1.png", width = 12, height = 6
)

# ------------------------------------------------------------------------------
# Save dataset to disk

#saveRDS(object_small_1, "BAL-macrophages_small_1.Rds")

# end of document
################################################################################




layout_dings <- "
AAABCD
AAABCD
AAABCD
EEEF#G
HHHF#G
"
layout_dings <-
  "
AAAAAABBCCDD
AAAAAABBCCDD
AAAAAABBCCDD
AAAAAABBCCDD
AAAAAABBCCDD
AAAAAABBCCDD
EEEEEEFFFGGG
EEEEEEFFFGGG
GGGGGGFFFGGG
GGGGGGFFFGGG
IIIIIIFFFGGG
IIIIIIFFFGGG
"

pc1_pc2_titeled+#A
  CD163_pc2+#B
  LGMN_pc2+#C
  LYZ_pc2 +#D
  HLA_d_pca_1+#E
  plot_spacer()+#F
  S100A8_pca_1+#G
  
  plot_layout(design = layout_dings)

pc1_pc2_titeled+#A
  CD163_pc2+#B
  LGMN_pc2+#C
  LYZ_pc2 +#D
  HLA_d_pca_1+#E
  PC_1_lollypop+#F
  PC_2_lollypop+#G
  S100A8_pca_1+#H
  
  plot_layout(design = layout_dings, guides = "collect")

layout_dings <- "
AAABCD
AAABCD
AAABCD
EEEFFF
GGGFFF
"
a2 <- pc1_pc2_titeled+#A
  CD163_pc2+#B
  LGMN_pc2+#C
  LYZ_pc2 +#D
  HLA_d_pca_1+#E
  x+#F
  S100A8_pca_1+#G
  
  plot_layout(design = layout_dings, guides = "collect") +plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 8))

a1/a2 +
  plot_layout(heights =   c(1, 2))
ggsave("try_3.png", width = 19, height = 20, units = "cm")

```