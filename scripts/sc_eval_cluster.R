library(Matrix)
library(parallel)

setwd('~/sc_eval/data/')

title_sub <- 'Neurodegeneration'
CombinedData <- readRDS(sprintf("%s.rds", title_sub))
rowSums(sapply(CombinedData[1], function(x){dim(counts(x, 'raw'))}))

cl <- makeCluster(12)
clust_res <- parSapply(cl, CombinedData[1:12], simplify = FALSE, function(local_sceval){
  library(Matrix)

  suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
  suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
  
  wd <- sprintf('~/sc_eval/data/%s', local_sceval@title)
  setwd(wd)
  
  ## Dimension Reduction ##
  local_clust_res <- list()
  for(d_type in c('raw', 'scimpute', 'drimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca_nb', 'dca_zinb')){
      for(r_type in c('umap', 'tsne', 'vae', 'dm', 'paga_umap')){
        if((d_type == 'scimpute' || d_type == 'drimpute') && (n_type == 'dca_nb' || n_type == 'dca_zinb'))
          next
        
        tag <- sprintf('%s_%s_%s', d_type, n_type, r_type)
        dim_res <- readRDS(sprintf('%s_%s.rds', local_sceval@title, tag))
        rownames(dim_res) <- paste0('cell.', 1:nrow(dim_res))
        
        # Check NA
        if(sum(apply(dim_res, 2, function(x){any(is.na(x))})) > 0)
          next
        
        ## Cluster ##
        write.table(dim_res, file = 'ReducedDimension.txt', col.names = FALSE, row.names = FALSE, sep = '\t')
        system(sprintf('python3 /home/arda/sc_eval/scripts/sc_eval_cluster.py %s', wd))
        local.clust <- read.table('LeidenClusters.txt', header=FALSE)
        local.clust <- setNames(local.clust$V1, rownames(dim_res))
        local_clust_res[[length(local_clust_res) + 1]] <- local.clust
        names(local_clust_res)[length(local_clust_res)] <- sprintf('%s_%s_%s', d_type, n_type, r_type)
      }
    }
  }
  return(local_clust_res)
})
stopCluster(cl)
names(clust_res) <- sapply(CombinedData, function(x){x@title})

title_sub <- 'Alectinib'
saveRDS(clust_res, file = sprintf('~/sc_eval/data/%s_Clusters.rds', title_sub))
quit(save='no')


# Plot Dim Reduced
library(Matrix)
library(ggplot2)
library(data.table)
library(dplyr)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

setwd('~/sc_eval/data/v2/')

title_sub <- 'TKI_Treatment'
CombinedData <- readRDS(sprintf("%s.rds", title_sub))

## TKI Treatment ##
uniq.order <- c('H3122-1', 'H3122-2', 
                'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor-1', 'H3122-erLor',
                'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz-3', 'H3122-erCz')
uniq.labels <- c('Naive-1', 'Naive-2', 
                 'Alec., 4h', 'Alec., 48h', 'Alec., 3w', 'erAlec., 1', 'erAlec., 2',
                 'Lor., 48h', 'Lor., 3w', 'erLor., 1', 'erLor., 2',
                 'Criz., 48h', 'Criz., 3w', 'erCriz., 1', 'erCriz., 2')

## Neurodegeneration ##
uniq.labels <- c('Week-0', 'Week-1', 'Week-2', 'Week-6')
uniq.order <- c('Week-0', 'Week-1', 'Week-2', 'Week-6')


## E2 Treatment ##
uniq.labels <- c('0h', '3h', '6h', '12h')
uniq.order <- c('0h', '3h', '6h', '12h')

## Pancreatic Differentiation ##
uniq.labels <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')
uniq.order <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')

#scale_color_manual(values=uniq.colors, labels=setNames(c(, uniq.labels[1:7])) +
#scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', ), uniq.labels[c(1,2,8,9,10,11)])) +
#scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', ), uniq.labels[c(1,2,12,13,14,15)])) +
#scale_color_manual(values=uniq.colors) +


#uniq.colors <- setNames(ggsci::pal_jco()(4), uniq.labels)

#uniq.colors <- setNames(ggsci::pal_jco()(4), uniq.labels)
temp_idx <- 12
clust_res <- readRDS('~/sc_eval/data/v2/Alectinib_Clusters.rds')
local_sceval <- CombinedData[[temp_idx]]

wd <- sprintf('~/sc_eval/data/v2/%s', local_sceval@title)
setwd(wd)
local_annot <- annot(local_sceval)

tag <- sprintf('%s_%s_%s', 'raw', 'dca_nb', 'paga_umap')
dim_res <- readRDS(sprintf('%s_%s.rds', local_sceval@title, tag))
plot_ft <- data.table('id' = local_annot$Id,
                      'Label' = local_annot$Label,
                      'Color' = local_annot$Color,
                      'Dim.1' = dim_res[,1],
                      'Dim.2' = dim_res[,2])
plot_ft$LabelFT <- uniq.labels[match(plot_ft$Label, table=uniq.order)]
temp <- plot_ft %>%
  group_by(Label, LabelFT, Color) %>%
  summarise(Label=Label, Color=Color) %>%
  distinct() %>%
  arrange(match(Label, uniq.order))

p <- ggplot(plot_ft, aes(x=Dim.1, y=Dim.2, fill=I(Color))) +
  geom_point(shape=21, size=1) +
  theme_classic() +
  labs(color='') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  theme(axis.title = element_text(face='bold', size=12)) +
  scale_fill_identity(guide='legend',
                      breaks=temp$Color,
                      labels=temp$LabelFT) +
  guides(fill=guide_legend(title='', override.aes = list(size=2)))

ggsave(p, filename = '~/sc_eval/plots/cluster_dim/Alectinib_SampleUMAPPAGA.pdf', width = 8, height = 6)

local_clust_res <- clust_res[[temp_idx]]
data.1 <- local_clust_res[[15]]
data.2 <- local_clust_res[[3]]
annot_ft <- annot(local_sceval)
plot_ft <- data.table('id' = annot_ft$Id,
                      'Label' = annot_ft$Label,
                      'Cluster1' = paste0('Cluster.', data.1),
                      'Cluster2' = paste0('Cluster.', data.2),
                      'Dim.1' = dim_res[,1],
                      'Dim.2' = dim_res[,2])
plot_ft$LabelFT <- uniq.labels[match(plot_ft$Label, table=uniq.order)]
p <- ggplot(plot_ft, aes(x=Dim.1, y=Dim.2, fill=Cluster1)) +
  geom_point(shape=21, size=1) +
  theme_classic() +
  labs(color='') +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  theme(axis.title = element_text(face='bold', size=12))

require(ComplexHeatmap)
require(circlize)

mat_1 <- plot_ft %>%
  group_by(LabelFT, Cluster1) %>%
  summarise(Count=n())
mat_2 <- plot_ft %>%
  group_by(LabelFT, Cluster2) %>%
  summarise(Count=n())

mat_1 <- reshape2::dcast(mat_1, LabelFT~Cluster1, value.var = 'Count', fill = 0)
mat_2 <- reshape2::dcast(mat_2, LabelFT~Cluster2, value.var = 'Count', fill = 0)
mat_1 <- data.matrix(mat_1[match(uniq.labels[1:7], table=mat_1$LabelFT),-1])
mat_2 <- data.matrix(mat_2[match(uniq.labels[1:7], table=mat_2$LabelFT),-1])
rownames(mat_1) <- uniq.labels[1:7]
rownames(mat_2) <- uniq.labels[1:7]

col_fun <- colorRamp2(breaks=seq(0,1,0.001), colors=colorRampPalette(c('white', 'firebrick'))(length(seq(0,1,0.001))))
rowNorm <- function(x){
  return(as.matrix(Matrix::Diagonal(x=Matrix::rowSums(x)^-1) %*% x))
}


h_1 <- Heatmap(rowNorm(mat_1),
               col = col_fun,
               cell_fun = function(j, i, x, y, width, height, fill){
                 if(round(rowNorm(mat_1)[i, j], digits=2) > 0){
                   grid.text(sprintf("%.2f", round(rowNorm(mat_1)[i, j], digits = 2)), x, y, gp = gpar(fontsize = 10))  
                 }else{
                   grid.text('', x, y, gp = gpar(fontsize = 10))
                 }
               },
               cluster_rows = FALSE,
               cluster_columns = TRUE,
               show_column_names = FALSE,
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize=12, fontface='bold'),
               show_heatmap_legend = FALSE,
               show_column_dend = FALSE,
               na_col = 'white',
               column_title = 'Workflow 1',
               rect_gp = gpar(col='black', lty=2),
               column_title_gp = gpar(fontsize=16, fontface='bold'))

cairo_pdf('~/sc_eval/plots/cluster_dim/Workflow1.pdf', width = 12, height = 3)
draw(h_1)
dev.off()

h_2 <- Heatmap(rowNorm(mat_2),
               col = col_fun,
               cell_fun = function(j, i, x, y, width, height, fill){
                 if(round(rowNorm(mat_2)[i, j], digits=2) > 0){
                   grid.text(sprintf("%.2f", round(rowNorm(mat_2)[i, j], digits = 2)), x, y, gp = gpar(fontsize = 10))
                 }else{
                   grid.text('', x, y, gp = gpar(fontsize = 10))
                 }
               },
               cluster_rows = FALSE,
               cluster_columns = TRUE,
               show_column_names = FALSE,
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize=12, fontface='bold'),
               show_heatmap_legend = FALSE,
               show_column_dend = FALSE,
               na_col = 'white',
               column_title = 'Workflow 2',
               rect_gp = gpar(col='black', lty=2),
               column_title_gp = gpar(fontsize=16, fontface='bold'))

cairo_pdf('~/sc_eval/plots/cluster_dim/Workflow2.pdf', width = 12, height = 3)
draw(h_2)
dev.off()
