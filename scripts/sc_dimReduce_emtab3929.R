library(uwot)
library(Rfast)
library(ggsci)

all_labels <- c('embryonic day 3', 'embryonic day 4', 'embryonic day 5', 'embryonic day 6', 'embryonic day 7')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 5, 'Set2')[as.numeric(as.factor(all_labels))], all_labels)

setwd('~/sc_eval/data/emtab3929_data/')

for(data_name in c('emtab3929')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('sctransform', 'deconv')){
      setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res', data_name, base_name, norm_name))
      sc_data <- readRDS(sprintf('emtab3929_%s_%s.rds', base_name, norm_name))
      
      ## Write Exp ##
      write.table(round(as.matrix(t(sc_data$norm_mat)), digits = 3),
                  file = sprintf('dhaka_res/emtab3929_%s_%s_mat.txt', base_name, norm_name),
                  sep = '\t',
                  col.names = FALSE,
                  row.names = FALSE,
                  append = FALSE,
                  quote = FALSE)
      
      ## UMAP ##
      umap_res <- umap(round(as.matrix(t(sc_data$norm_mat))),
                       n_neighbors = 20,
                       n_components = 3,
                       n_epochs = 1000,
                       scale = TRUE,
                       init = 'spectral',
                       min_dist = 0.01,
                       spread = 1,
                       pca = NULL,
                       n_threads = 20)
      
      ## HTML Output ##
      plot_data_umap <- data.frame('Dimension.1' = umap_res[,1],
                                   'Dimension.2' = umap_res[,2],
                                   'Dimension.3' = umap_res[,3],
                                   'Label' = sc_data$labels$Label,
                                   'Id' = sc_data$labels$Id,
                                   'Colors' = common_colors[match(sc_data$labels$Label, table = names(common_colors))],
                                   stringsAsFactors = FALSE)
      write.table(plot_data_umap,
                  row.names = FALSE,
                  col.names = TRUE,
                  sep = '\t',
                  append = FALSE,
                  quote = FALSE,
                  file = sprintf('umap_res/emtab3929_%s_%s_umap.txt', base_name, norm_name))
      
      require(plotly)
      plot_ly(plot_data_umap,
              x = ~Dimension.1,
              y = ~Dimension.2,
              z = ~Dimension.3,
              text = plot_data_umap$Label,
              hoverinfo = 'text',
              hoverlabel = list(bgcolor = I(plot_data_umap$Colors)),
              marker = list(size = 3,
                            color = I(plot_data_umap$Colors))) %>%
        add_markers() %>%
        htmlwidgets::saveWidget(file = sprintf('emtab3929_%s_%s_umap.html', base_name, norm_name), selfcontained = TRUE)
    }
  }
}
