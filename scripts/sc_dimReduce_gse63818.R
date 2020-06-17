library(uwot)
library(Rfast)
library(ggsci)

setwd('~/sc_eval/data/gse63818/')
all_data <- readRDS('gse63818_norm.rds')
common.colors <- ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all_data$labels$time))]

## UMAP ##
umap.res <- umap(as.matrix(t(all_data$norm_mat)),
                 n_neighbors = 5,
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
                             'Label' = all_data$labels$time,
                             'Id' = colnames(all_data$norm_mat),
                             'Colors' = common.colors,
                             stringsAsFactors = FALSE)
write.table(plot.data.umap,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t',
            append = FALSE,
            quote = FALSE,
            file = 'umap_res/GSE63818_UMAP.txt')

require(plotly)
plot_ly(plot.data.umap,
        x = ~Dimension.1,
        y = ~Dimension.2,
        z = ~Dimension.3,
        text = plot.data.umap$Label,
        hoverinfo = 'text',
        hoverlabel = list(bgcolor = I(common.colors)),
        marker = list(size = 2,
                      color = I(common.colors))) %>%
  add_markers() %>%
  htmlwidgets::saveWidget(file = 'GSE63818_UMAP.html', 
                          selfcontained = TRUE)

## Dhaka ##
embedding <- read.table('~/sc_eval/data/gse63818/dhaka_res/GSE63818_Dhaka.txt.txt', header = FALSE, sep = '')
plot.data <- data.frame('Dimension.1' = embedding[,1],
                        'Dimension.2' = embedding[,2],
                        'Dimension.3' = embedding[,3],
                        'Label' = all_data$labels$time,
                        'Id' = colnames(all_data$norm_mat),
                        'Colors' = common.colors,
                        stringsAsFactors = FALSE)
write.table(plot.data, file = '~/sc_eval/data/gse63818/dhaka_res/GSE63818_DHAKA_META.txt', col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE, append = FALSE)
