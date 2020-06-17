library(Matrix)
library(mclust)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

title_sub <- 'Alectinib'


clust_res <- readRDS(sprintf('~/sc_eval/data/%s_Clusters.rds', title_sub))

## ARI ##
ari_res <- sapply(clust_res, simplify = FALSE, function(x){
  local_ari_res <- matrix(NA, ncol=length(x), nrow=length(x))
  for(i in 1:(nrow(local_ari_res)-1)){
    for(j in (i+1):ncol(local_ari_res)){
      local_ari_res[i,j] <- adjustedRandIndex(x[[i]], x[[j]])
    }
  }
  return(local_ari_res)
})

## Plot ##
require(ComplexHeatmap)
require(circlize)
require(viridis)

for(i in 1:length(ari_res)){
  local_res <- ari_res[[i]]
  colnames(local_res) <- names(clust_res[[1]])
  rownames(local_res) <- names(clust_res[[1]])
  
  col_fun <- colorRamp2(breaks=seq(0,1,0.001), colors=viridis(length(seq(0,1,0.001))))
  
  h <- Heatmap(local_res, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_column_names = TRUE, 
               show_row_names = TRUE,
               column_names_gp = gpar(fontsize=10),
               row_names_gp = gpar(fontsize=10),
               column_names_rot = 60,
               column_title=sprintf('%s S%d', title_sub, i),
               column_title_gp = gpar(fontsize=16, fontface='bold'),
               rect_gp = gpar(col = "white", lwd = 2),
               column_names_side = 'top',
               col = col_fun, 
               na_col = 'white',
               heatmap_legend_param = list(title='Adjusted Rand Index', title_position='leftcenter-rot', legend_height=unit(6,'cm')))
  pdf(sprintf('~/sc_eval/plots/ClusterComparison_%s_S%d.pdf', title_sub, i), width = 12, height = 10)
  draw(h)
  dev.off()
}


imp_res <- readRDS(sprintf('~/sc_eval/data/%s_Importance.rds', title_sub))

## Within workflow ##
comp_res <- sapply(1:length(imp_res), simplify = FALSE, function(x.idx){
  x <- imp_res[[x.idx]]
  local_comp_res <- sapply(c(25, 50, 100, 250, 500, 1000), simplify = FALSE, function(thr){
    comp_mat <- matrix(NA, ncol=7, nrow=7)
    for(i in 1:6){
      for(j in (i+1):7){
        feat.x <- names(sort(x[[i]], decreasing = TRUE))[1:thr]
        feat.y <- names(sort(x[[j]], decreasing = TRUE))[1:thr]
        comp_mat[i,j] <- length(intersect(feat.x, feat.y))/length(union(feat.x, feat.y))
      }
    }
    return(na.omit(as.vector(comp_mat)))
  })
  local_comp_res <- do.call('rbind', sapply(1:length(local_comp_res), simplify = FALSE, function(i){
    return(data.frame('Jacc'=local_comp_res[[i]],
                      'Method'=names(imp_res)[x.idx],
                      'Thr'=c(25,50,100,250,500,1000)[i]))
  }))
  return(local_comp_res)
})
comp_res <- do.call('rbind', comp_res)
comp_res$Thr <- paste0('Thr-', comp_res$Thr)
comp_res$Thr <- factor(comp_res$Thr, levels = paste0('Thr-', c(25,50,100,250,500,1000)))

comp_res$Data <- ifelse(grepl(pattern = 'raw', comp_res$Method), 'Raw', 'ScImpute')
comp_res$Normalization <- ifelse(grepl(pattern = 'deconv', comp_res$Method), 'Deconvolution',
                                 ifelse(grepl(pattern = 'sctransform', comp_res$Method), 'ScTransform',
                                        ifelse(grepl(pattern = 'dca', comp_res$Method), 'DCA', 'Unknown')))
comp_res$Dimension <- ifelse(grepl(pattern = 'paga_umap', comp_res$Method), 'Paga_UMAP',
                             ifelse(grepl(pattern = 'umap', comp_res$Method), 'UMAP',
                                    ifelse(grepl(pattern = 'dhaka', comp_res$Method), 'Dhaka',
                                           ifelse(grepl(pattern = 'dm', comp_res$Method), 'DM',
                                                  ifelse(grepl(pattern = 'tsne', comp_res$Method), 't-SNE', 'Unknown')))))

require(ggplot2)
p <- ggplot(comp_res, aes(x=Thr, y=Jacc, fill=Data)) +
  geom_boxplot() +
  theme_minimal() +
  ylab('Jaccard Index') +
  xlab('') +
  labs(fill='') +
  ggtitle(sprintf('%s', title_sub)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(size=14, face='bold', hjust = 0.5),
        strip.text = element_text(size=12),
        legend.key.size = unit(0.5, 'cm')) +
  facet_grid(Normalization~Dimension)
ggsave(p, filename = sprintf('~/sc_eval/plots/%s_FeatureComparison.pdf', title_sub), width = 8, height = 6)
