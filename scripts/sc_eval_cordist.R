library(Matrix)
library(ggplot2)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

title_sub <- 'GSE107863'

## Within workflow ##
all_res <- list()
for(d_type in c('raw', 'scimpute')){
  for(n_type in c('deconv', 'sctransform', 'dca')){
    for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      
      traject_res <- sapply(1:7, simplify = FALSE, function(i){
        title <- paste0(title_sub, sprintf('_S%d', i))
        tag <- sprintf('~/sc_eval/data/%s/%s_%s_%s_%s_Slingshot.rds', title, title, d_type, n_type, r_type)
        traj_res <- readRDS(tag)[[1]]
        return(subset(traj_res, traj_res$Type == 'CB'))
      })
      comp_res <- sapply(1:(length(traject_res)-1), simplify = FALSE, function(i){
        res.1 <- traject_res[[i]]
        ncurv.1 <- length(grep(pattern='curve', colnames(res.1)))
        local_comp_res <- sapply((i+1):length(traject_res), simplify = FALSE, function(j){
          res.2 <- traject_res[[j]]
          ncurv.2 <- length(grep(pattern='curve', colnames(res.2)))
          
          rank_cor <- sapply(1:ncurv.1, simplify = FALSE, function(k){
            temp.1 <- na.omit(setNames(res.1[,which(colnames(res.1) == paste0('curve', k))], res.1$Id))
            local_rank_cor <- sapply(1:ncurv.2, simplify = FALSE, function(z){
              temp.2 <- na.omit(setNames(res.2[,which(colnames(res.2) == paste0('curve', z))], res.2$Id))
              common.ids <- intersect(names(temp.1), names(temp.2))
              if(length(common.ids) < 10)
                return(NA)
              
              temp.1.filt <- temp.1[match(common.ids, table=names(temp.1))]
              temp.2.filt <- temp.2[match(common.ids, table=names(temp.2))]
              abs(cor(temp.1.filt, temp.2.filt, method='spearman'))
            })
            return(unlist(local_rank_cor))
          })
          return(na.omit(unlist(rank_cor)))
        })
        return(local_comp_res)
      })
      comp_res <- unlist(comp_res, recursive = TRUE)
      all_res[[length(all_res)+1]] <- data.frame('Cor' = comp_res,
                                                 'Method' = sprintf("%s_%s_%s", d_type, n_type, r_type))
    }
  }
}
comp_res <- do.call('rbind', all_res)

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
p <- ggplot(comp_res, aes(x=Data, y=Cor, fill=Normalization)) +
  geom_boxplot() +
  theme_minimal() +
  ylab('Rank Correlation') +
  xlab('') +
  labs(fill='') +
  ggtitle(sprintf('%s', title_sub)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(size=14, face='bold', hjust = 0.5),
        strip.text = element_text(size=12),
        legend.key.size = unit(0.5, 'cm')) +
  facet_grid(~Dimension)
ggsave(p, filename = sprintf('~/sc_eval/plots/%s_SlingshotCorrelation.pdf', title_sub), width = 8, height = 4)
