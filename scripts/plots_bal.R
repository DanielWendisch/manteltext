library(patchwork)
library(cowplot)
library(tidyverse)
library(Seurat)
library(magrittr)
library(ggpubr)
library(gridGraphics)
library(ggdendro)

source("functions_bal_macrophages.R")
source("initialize_bal_data.R")

########




layout_dings <- "
AAABCD
AAABCD
AAABCD
EEEFFF
GGGFFF
"

pc1_pc2_titeled+#A
  CD163_pc2+#B
  LGMN_pc2+#C
  LYZ_pc2 +#D
  HLA_d_pca_1+#E
  plot_spacer()+#F
  S100A8_pca_1+#G
  
  plot_layout(design = layout_dings)

############





colnames(dings)


library(Seurat)
library(tidyverse)
library(janitor)
library(magrittr)
library(effsize)

#source()
#source_rmd( "BAL_effect_size.Rmd")



cohens_d_module_scores_fibrosis <- tibble(geneset=character(),cluster=character(),cohens_d=numeric())

data <- data_fibrosis
for (gs in levels(data$Geneset)) {
  print(gs)
  for (clus in levels(data$Cluster) ) { #data %>% filter(Geneset==gs) %>% pull(Cluster) %>% unique()
    #print(clus)
    ctrl <- data%>% filter(Geneset==gs, Cluster!=clus) %>% pull(Score)
    treat <- data%>% filter(Geneset==gs, Cluster==clus) %>% pull(Score)
    d <- cohen.d(treat, ctrl)
    cohens_d_module_scores_fibrosis %<>%  add_row(geneset=gs,cluster=clus,cohens_d=d$estimate)
  }
}
cohens_d_module_scores_fibrosis

eff_size_plot <- cohens_d_module_scores_fibrosis %>% ggplot(aes(x=geneset, y=cluster, fill=cohens_d)) + geom_tile()+theme(axis.text.x=element_text(angle = 90))+scale_fill_viridis_c()
eff_size_plot

p_ayaub_fib_score_point <- ggplot2::ggplot(
  data = data_fibrosis%>% filter(Geneset=="Ayau_IPFeMp"),
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +score_violin_theme(limits = limits, p_score=ps%>% filter(Geneset=="Ayau_IPFeMp"))


p_ayaub_fib_score_point <- ggplot2::ggplot(
  data = data_fibrosis%>% filter(Geneset=="Ayau_IPFeMp"),
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +score_violin_theme(limits = limits, p_score=ps%>% filter(Geneset=="Ayau_IPFeMp"))


BAL_macrophages$Cluster %<>% factor(levels = c("Monocytes" , "MonoM", "CD163LGMNM" ,"AM1",  "AM2", "ProlifAM" ) )
Idents(BAL_macrophages) <- "Cluster"
p_ayaub_fib_score_dot <- DotPlot(BAL_macrophages, features = genesets_fibrosis$Ayau_IPFeMp %>% ensemblizer()) + scale_y_discrete(guide = guide_axis(angle = 90)) +coord_flip() 
p_ayaub_fib_score_dot %<>% seurat_plot_renamer()

ayaub_score_plus_features_plot <- list(p_ayaub_fib_score_point,p_ayaub_fib_score_dot)
plot_grid(plotlist=ayaub_score_plus_features_plot, ncol=1, align='v')


p_ayaub_fib_score_point/p_ayaub_fib_score_dot


