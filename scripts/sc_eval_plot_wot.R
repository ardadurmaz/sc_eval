library(Matrix)
library(ggplot2)
library(dbscan)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)


get_nn <- function(x=NULL, k=15){
  local.d <- as.matrix(dist(x))
  diag(local.d) <- NA
  local.e <- do.call('rbind', sapply(1:nrow(local.d), simplify = FALSE, function(i){
    local.c <- colnames(local.d)[order(local.d[i,], decreasing = TRUE)[1:k]]
    return(data.frame('cell.a' = rep(rownames(local.d)[i], length(local.c)),
                      'cell.b' = local.c))
  }))
  return(local.e)
}

get_clust <- function(edge.ft=NULL){
  require(igraph)
  local.net <- graph_from_data_frame(edge.ft, directed = FALSE)
  clust.res <- cluster_louvain(local.net)
  return(membership(clust.res))
}

CombinedData <- readRDS('data/EGEOD103334_Data.rds')

title_sub <- 'EGEOD103334'
for(d_type in c('raw', 'scimpute')){
  for(n_type in c('deconv', 'sctransform', 'dca')){
    for(r_type in c('paga_umap')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      ## Read Data ##
      cat(sprintf('Processing: T:%s D:%s N:%s\n', title_sub, d_type, n_type))
      dim_res_c <- sapply(1:7, simplify = FALSE, function(i){
        title <- sprintf('%s_S%d', title_sub, i)
        dim_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_%s_%s_%s.rds', title, title, d_type, n_type, r_type))  
        
        ## Expression ##
        tag.expr <- sprintf('%s_%s', d_type, n_type)
        if(tag.expr == 'raw_deconv'){
          expr_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_RawDeconv.rds', title, title))  
        }else if(tag.expr == 'raw_sctransform'){
          expr_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_RawScTransform.rds', title, title))  
        }else if(tag.expr == 'raw_dca'){
          expr_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_RawDCA.rds', title, title))  
        }else if(tag.expr == 'scimpute_deconv'){
          expr_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_ScImputeDeconv.rds', title, title))  
        }else if(tag.expr == 'scimpute_sctransform'){
          expr_res <- readRDS(sprintf('~/sc_eval/data/%s/%s_ScImputeScTransform.rds', title, title))  
        }
        rownames(dim_res) <- colnames(expr_res)
        return(dim_res)
      })
      
      ## Choose Cell Sets ##
      local.nn <- get_nn(x=dim_res_c[[7]])
      local.clust <- get_clust(local.nn)
      
      ## Get Annotation ##
      title_idx <- which(sapply(CombinedData, function(x){x@title}) == sprintf('%s_S7', title_sub))
      annot_ft <- annot(CombinedData[[title_idx]])
      local.clust.ft <- data.frame('id' = names(local.clust),
                                   'cluster' = paste0('cluster.', local.clust),
                                   'day' = annot_ft$Days[match(names(local.clust), table=annot_ft$Id)])
      cell_sets <- sapply(unique(local.clust.ft$cluster), simplify = FALSE, function(l.c){
        local.local.clust.ft <- subset(local.clust.ft, local.clust.ft$cluster == l.c)
        local_cell_sets <- sapply(unique(local.local.clust.ft$day), simplify = FALSE, function(l.d){
          temp <- subset(local.local.clust.ft, local.local.clust.ft$day == l.d)
          if(nrow(temp) < 5)
            return(NULL)
          
          return(data.frame('cluster' = sprintf('%s_D%f', l.c, l.d),
                            'cells' = paste(as.character(temp$id), collapse = ',')))
        })
        local_cell_sets <- do.call('rbind', local_cell_sets)
        return(local_cell_sets)
      })
      cell_sets <- do.call('rbind', cell_sets)
      wot_res <- sapply(1:7, simplify = FALSE, function(s){
        local_wot_res <- sapply(1:nrow(cell_sets), function(i){
          local_tag <- cell_sets$cluster[i]
          local_cells <- cell_sets$cells[i]
          local_cells <- gsub(pattern = ',', replacement = '\t', local_cells)
          cat(sprintf('%s\t%s', local_tag, local_cells), file = '~/sc_eval/data/CellSet.gmt', sep = '\n', append=FALSE)
          title <- sprintf('%s_S%d', title_sub, s)
          comm <- sprintf('/usr/bin/python3 /home/arda/sc_eval/scripts/sc_eval_wot_traject.py /home/arda/sc_eval/data/%s/%s_%s_%s_WOT %s',
                          title, title, d_type, n_type, gsub(pattern = '^D', replacement = '', unlist(strsplit(as.character(local_tag), split = '_'))[2]))
          system(comm)
          tr_res <- read.table('~/sc_eval/data/TrajectoryWOT.txt', header=FALSE, sep = '\t')
          return(tr_res$V1)
        })
        return(local_wot_res)
      })
      for(i in 1:7){
        rownames(wot_res[[i]]) <- rownames(dim_res_c[[i]])
      }
      saveRDS(wot_res, file = sprintf('~/sc_eval/data/%s_%s_%s_WOT_results.rds', title_sub, d_type, n_type))
    }
  }
}
quit(save='no')

comp_res_all <- list()
for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'EGEOD103334', 'GSE107863', 'GSE107858')){
  comp_res <- list()
  for(d_type in c('raw', 'scimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      wot_res <- readRDS(sprintf('~/sc_eval/data/%s_%s_%s_WOT_results.rds', title_sub, d_type, n_type))
      cor_res <- sapply(1:ncol(wot_res[[1]]), simplify = TRUE, function(i){
        val_res <- sapply(wot_res, simplify = FALSE, function(j){
          j[,i]
        })
        common_cells <- Reduce(intersect, sapply(val_res, simplify = FALSE, names))
        val_mat <- do.call('cbind', sapply(val_res, simplify = FALSE, function(x){
          return(x[match(common_cells, table=names(x))])
        }))
        cor_val <- cor(val_mat, method = 'spearman')
        diag(cor_val) <- NA
        return(mean(as.vector(cor_val), na.rm=TRUE))
      })
      comp_res[[length(comp_res)+1]] <- cor_res
      names(comp_res)[length(comp_res)] <- sprintf('%s_%s', d_type, n_type)
    }
  }
  comp_res_ft <- do.call('rbind',
                         sapply(1:length(comp_res), simplify = FALSE, function(i){
                           return(data.frame('spearman' = comp_res[[i]],
                                             'name' = names(comp_res)[i]))
                         }))
  comp_res_ft$data <- rep(title_sub, nrow(comp_res_ft))
  comp_res_all[[length(comp_res_all) + 1]] <- comp_res_ft
}
comp_res_all <- do.call('rbind', comp_res_all)


require(ggplot2)
comp_res_all$data <- factor(comp_res_all$data,
                            levels = c('Alectinib', 'Lorlatinib', 'Crizotinib', 'EGEOD103334', 'GSE107863', 'GSE107858'))
p <- ggplot(comp_res_all, aes(x=name, y=spearman, fill=data)) +
  geom_boxplot(outlier.size = 1) +
  geom_point(aes(color=data), shape=20, alpha=0.4, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0)) +
  theme_classic() +
  labs(fill='', color='') +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(face='bold', size=10, angle = 30, hjust = 1),
        axis.text.y = element_text(face='bold', size=10),
        plot.margin = unit(c(0.5,0.5,0.5,1), 'cm'),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8)) +
  scale_x_discrete(labels=c('raw_deconv'= 'Raw Deconvolution', 'raw_sctransform'='Raw ScTransform', 'raw_dca'='Raw DCA', 'scimpute_deconv'='ScImpute Deconvolution', 'scimpute_sctransform'='ScImpute ScTransform'))
ggsave(p, filename = '~/WOT_Comparison.pdf', device = cairo_pdf, width = 6, height = 5)

CombinedData <- readRDS('data/CombinedData.rds')
local_sceval <- CombinedData[[7]]

norm_mat <- readRDS('data/Alectinib_S7/Alectinib_S7_RawDeconv.rds')
dim_res <- readRDS('data/Alectinib_S7/Alectinib_S7_raw_deconv_umap.rds')
rownames(dim_res) <- colnames(norm_mat)
annot_ft <- local_sceval@annot


plot_ft_time <- data.frame('Dim.1' = dim_res[,1],
                           'Dim.2' = dim_res[,2],
                           'Label' = annot_ft$Label[match(colnames(norm_mat), table=annot_ft$Id)])
plot_ft_time$Label <- factor(plot_ft_time$Label,
                             levels = c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec'))
uniq.colors <- setNames(ggsci::pal_simpsons()(7), unique(plot_ft_time$Label))

p.time <- ggplot(plot_ft_time, aes(x=Dim.1, y=Dim.2, color=Label)) +
  geom_point(shape=20) +
  theme_classic() +
  scale_color_manual(values=uniq.colors,
                     labels=setNames(c('Naive-1', 'Naive-2', 'Alec., 4h', 'Alec., 48h', 'Alec., 3w', 'erAlec., 1', 'erAlec., 2'),
                                     unique(plot_ft_time$Label)))
ggsave(p.time, filename = '~/sc_eval/plots/Alectinib_S7_RawDeconv_Label.pdf', width = 7, height = 5)

clust_res <- kmeans(dim_res, centers = 8)
clust_res <- split(names(clust_res$cluster), as.numeric(clust_res$cluster))
names(clust_res) <- paste0('K', names(clust_res))
for(i in c(1,8)){
  local.cells <- clust_res[[i]]
  local.cells <- intersect(local.cells, annot_ft$Id[annot_ft$Days == 21])
  local.s <- sprintf('%s\t-\t%s', names(clust_res)[i], paste(local.cells, collapse = '\t'))
  cat(local.s, file = '~/Alectinib_S7_Clusters.gmt', append=TRUE, sep = "\n")
}
plot_ft <- data.frame('Dim.1' = dim_res[,1],
                      'Dim.2' = dim_res[,2],
                      'Cluster' = paste0('K-', clust_res$cluster))
p.clust <- ggplot(plot_ft, aes(x=Dim.1, y=Dim.2, color=Cluster)) +
  geom_point() +
  theme_classic()

require(cowplot)
p <- plot_grid(p.time, p.clust)
save_plot(filename = '~/Alectinib_S7_Test.pdf', p, nrow = 1, base_height = 4, base_width = 12)

traject_dist <- read.table('~/Alectinib_S7_TestTrajectories.tsv', header=FALSE, sep = '\t')
rownames(traject_dist) <- colnames(norm_mat)

plot_ft_time$K1 <- traject_dist$V1
plot_ft_time$K8 <- traject_dist$V2
plot_ft_time$Day <- annot_ft$Days[match(rownames(plot_ft_time), table=annot_ft$Id)]
max_val <- max(plot_ft_time$K1)
day.mapping <- setNames(c(0, 0.15, 2, 21, 112), c('Naive', '4h', '48h', '3w', 'Resist'))

p.list <- sapply(unique(plot_ft_time$Day), simplify = FALSE, function(t){
  local_plot_ft_time <- plot_ft_time
  local_plot_ft_time$K1[local_plot_ft_time$Day != t] <- NA
  local.alpha <- ifelse(!is.na(local_plot_ft_time$K1), 1, 0.4)
  local_plot_ft_time$NormVal <- local_plot_ft_time$K1/max(local_plot_ft_time$K1, na.rm = TRUE)
  p.tr <- ggplot(local_plot_ft_time, aes(x=Dim.1, y=Dim.2, color=NormVal)) +
    geom_point(alpha=local.alpha, shape=20) +
    labs(color='Association\nProbability') +
    ggtitle(sprintf('Day: %s', names(day.mapping)[which(day.mapping == t)])) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face='bold')) +
    scale_color_viridis_c()
  ggsave(p.tr, filename = sprintf('~/Alectinib_K1_WOT_%.2f.pdf', t), width = 7, height = 5)
  return(p.tr)
})
combined.plot <- plot_grid(plotlist = p.list, ncol=3)
save_plot(filename = '~/Alectinib_S7_K8.pdf', plot = combined.plot, ncol = 3, base_height = 10, base_width = 8)
