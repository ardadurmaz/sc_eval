library(uwot)
library(Rfast)
library(ggsci)

setwd('~/sc_eval/data/cellranger_data/scimpute_res/sctransform_res')

all.labels <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1',
                'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3')
common.colors <- setNames(ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all.labels))],
                          all.labels)

## UMAP ##
for(d in c('alectinib', 'lorlatinib', 'crizotinib')){
  message('Processing Drug: ', d)
  if(d == 'alectinib'){
    in_file <- 'alec/cellranger_alectinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_alectinib_scimpute_sctransform_umap'
  }else if(d == 'lorlatinib'){
    in_file <- 'lor/cellranger_lorlatinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_lorlatinib_scimpute_sctransform_umap'
  }else{
    in_file <- 'criz/cellranger_crizotinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_crizotinib_scimpute_sctransform_umap'
  }
  sc.data <- readRDS(in_file)
  
  ## UMAP ##
  umap.res <- umap(as.matrix(t(sc.data$norm_mat)),
                   n_neighbors = 30,
                   n_components = 3,
                   n_epochs = 1000,
                   scale = TRUE,
                   init = 'spectral',
                   min_dist = 0.01,
                   spread = 1,
                   pca = NULL,
                   n_threads = 20)
  
  ## HTML Output ##
  local.colors <- common.colors[match(sc.data$labels, table = names(common.colors))]
  plot.data.umap <- data.frame('Dimension.1' = umap.res[,1],
                               'Dimension.2' = umap.res[,2],
                               'Dimension.3' = umap.res[,3],
                               'Label' = sc.data$labels,
                               'Id' = colnames(sc.data$norm_mat),
                               'Colors' = local.colors,
                               stringsAsFactors = FALSE)
  write.table(plot.data.umap,
              row.names = FALSE,
              col.names = TRUE,
              sep = '\t',
              append = FALSE,
              quote = FALSE,
              file = paste0(out_file, '.txt'))
  
  require(plotly)
  plot_ly(plot.data.umap,
          x = ~Dimension.1,
          y = ~Dimension.2,
          z = ~Dimension.3,
          text = plot.data.umap$Label,
          hoverinfo = 'text',
          hoverlabel = list(bgcolor = I(local.colors)),
          marker = list(size = 2,
                        color = I(local.colors))) %>%
    add_markers() %>%
    htmlwidgets::saveWidget(file = paste0(out_file, '.html'), 
                            selfcontained = TRUE)
}


## Dhaka ##
for(d in c('alectinib', 'lorlatinib', 'crizotinib')){
  message('Processing Drug: ', d)
  if(d == 'alectinib'){
    in_file <- 'alec/dhaka_res/cellranger_alectinib_scimpute_sctransform_dhaka.txt'
    meta_file <- 'alec/cellranger_alectinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_alectinib_scimpute_sctransform_dhaka_meta'
  }else if(d == 'lorlatinib'){
    in_file <- 'lor/dhaka_res/cellranger_lorlatinib_scimpute_sctransform_dhaka.txt'
    meta_file <- 'lor/cellranger_lorlatinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_lorlatinib_scimpute_sctransform_dhaka_meta'
  }else{
    in_file <- 'criz/dhaka_res/cellranger_crizotinib_scimpute_sctransform_dhaka.txt'
    meta_file <- 'criz/cellranger_crizotinib_scimpute_sctransform.rds'
    out_file <- 'cellranger_crizotinib_scimpute_sctransform_dhaka_meta'
  }
  embedding <- read.table(in_file, header = FALSE, sep = ' ')
  meta.data <- readRDS(meta_file)

  ## HTML Output ##
  local.colors <- common.colors[match(meta.data$labels, table = names(common.colors))]
  plot.data <- data.frame('Dimension.1' = embedding[,1],
                          'Dimension.2' = embedding[,2],
                          'Dimension.3' = embedding[,3],
                          'Label' = meta.data$labels,
                          'Id' = colnames(meta.data$norm_mat),
                          'Colors' = local.colors,
                          stringsAsFactors = FALSE)
  write.table(plot.data,
              row.names = FALSE,
              col.names = TRUE,
              sep = '\t',
              append = FALSE,
              quote = FALSE,
              file = paste0(out_file, '.txt'))
  
  require(plotly)
  plot_ly(plot.data,
          x = ~Dimension.1,
          y = ~Dimension.2,
          z = ~Dimension.3,
          text = plot.data$Label,
          hoverinfo = 'text',
          hoverlabel = list(bgcolor = I(local.colors)),
          marker = list(size = 2,
                        color = I(local.colors))) %>%
    add_markers() %>%
    htmlwidgets::saveWidget(file = paste0(out_file, '.html'), 
                            selfcontained = TRUE)
}
