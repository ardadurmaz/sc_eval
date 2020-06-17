library(uwot)
library(Rfast)
library(ggsci)

all_labels <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 7, 'Dark2'), all_labels)

setwd('~/sc_eval/data/gse87375/')

for(base_name in c('raw', 'scimpute')){
  for(type in c('a', 'b')){
    for(norm_meth in c('deconv', 'sctransform')){
      sc_data <- readRDS(sprintf('%s_res/%s_res/gse87375_%s_%s_%s.rds', base_name, norm_meth, type, base_name, norm_meth))
      
      ## Write Exp ##
      write.table(round(as.matrix(t(sc_data$norm_mat)), digits = 3),
                  file = sprintf('%s_res/%s_res/dhaka_res/gse87375_%s_%s_%s_mat.txt', base_name, norm_meth, type, base_name, norm_meth),
                  sep = '\t',
                  col.names = FALSE,
                  row.names = FALSE,
                  append = FALSE,
                  quote = FALSE)
      
      
      ## UMAP ##
      umap.res <- umap(as.matrix(t(sc_data$norm_mat)),
                       n_neighbors = 10,
                       n_components = 3,
                       n_epochs = 1000,
                       scale = TRUE,
                       init = 'spectral',
                       min_dist = 0.01,
                       spread = 1,
                       pca = NULL,
                       n_threads = 20)
      
      ## HTML Output ##
      plot.data.umap <- data.frame('Dimension.1' = umap.res[,1],
                                   'Dimension.2' = umap.res[,2],
                                   'Dimension.3' = umap.res[,3],
                                   'Label' = sc_data$labels$time,
                                   'Id' = sc_data$labels$id,
                                   'Colors' = common_colors[match(sc_data$labels$time, table = names(common_colors))],
                                   stringsAsFactors = FALSE)
      write.table(plot.data.umap,
                  row.names = FALSE,
                  col.names = TRUE,
                  sep = '\t',
                  append = FALSE,
                  quote = FALSE,
                  file = sprintf('%s_res/%s_res/umap_res/gse87375_%s_%s_%s_umap.txt', base_name, norm_meth, type, base_name, norm_meth))
      
      require(plotly)
      plot_ly(plot.data.umap,
              x = ~Dimension.1,
              y = ~Dimension.2,
              z = ~Dimension.3,
              text = plot.data.umap$Label,
              hoverinfo = 'text',
              hoverlabel = list(bgcolor = I(plot.data.umap$Colors)),
              marker = list(size = 3,
                            color = I(plot.data.umap$Colors))) %>%
        add_markers() %>%
        htmlwidgets::saveWidget(file = sprintf('gse87375_%s_%s_%s_umap.html', type, base_name, norm_meth), 
                                selfcontained = TRUE)
      
    }
  }
}