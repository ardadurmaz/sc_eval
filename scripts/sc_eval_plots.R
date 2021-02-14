library(Matrix)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(data.table)
library(dplyr)
library(tidyr)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))

## Calculate Cluster Comparison ##
set.seed(42)
data_tag <- 'Alectinib'
clust_res <- readRDS(sprintf('data/v2/%s_Clusters.rds', data_tag))
clust_stat <- sapply(clust_res, simplify = FALSE, function(x){
  comp_res <- matrix(NA, ncol=length(x), nrow=length(x), dimnames=list(names(x), names(x)))
  for(i in 1:(nrow(comp_res)-1)){
    clust_a <- x[[i]]
    for(j in i:ncol(comp_res)){
      clust_b <- x[[j]]
      comp_res[i,j] <- adjustedRandIndex(clust_a, clust_b)
    }
  }
  return(comp_res)
})
saveRDS(clust_stat, file = sprintf('%s_Cluster_Comparisons.rds', data_tag))


## Plot Cluster Comparison Boxplot ##
combined_ari <- sapply(c('Alectinib', 'Lorlatinib', 'Crizotinib', 'E2_Treatment_MCF7', 'E2_Treatment_T47D', 'Neurodegeneration', 'PancreaticMaturation_A', 'PancreaticMaturation_B'), simplify = FALSE, function(data_tag){
  cluster_comp <- readRDS(sprintf('data/v2/%s_Cluster_Comparisons.rds', data_tag))
  cluster_comp <- sapply(1:length(cluster_comp), simplify = FALSE, function(x.idx){
    x <- cluster_comp[[x.idx]]
    x_ft <- na.omit(reshape2::melt(x))
    colnames(x_ft) <- c('Method.a', 'Method.b', 'ARI')
    x_ft$Data <- rep(paste0('data.', x.idx), nrow(x_ft))
    return(x_ft)
  })
  cluster_comp <- rbindlist(cluster_comp)
  cluster_comp$DataTag <- rep(data_tag, nrow(cluster_comp))
  return(cluster_comp)
})
combined_ari <- rbindlist(combined_ari) %>%
  distinct() %>%
  filter(Method.a != Method.b)

combined_ari$DataTag <- factor(combined_ari$DataTag,
                               levels = c('Alectinib', 'Lorlatinib', 'Crizotinib',
                                          'Neurodegeneration', 'PancreaticMaturation_B', 'PancreaticMaturation_A',
                                          'E2_Treatment_T47D', 'E2_Treatment_MCF7'))
combined_ari$Data <- factor(combined_ari$Data,
                            levels = paste0('data.', 1:12))

## Linear Regression 
tki <- readRDS('~/sc_eval/data/v2/TKI_Treatment.rds')
neur <- readRDS('~/sc_eval/data/v2/Neurodegeneration.rds')
mcf7 <- readRDS('~/sc_eval/data/v2/E2_Treatment_MCF7.rds')
t47d <- readRDS('~/sc_eval/data/v2/E2_Treatment_T47D.rds')
panc.a <- readRDS('~/sc_eval/data/v2/PancreaticMaturation_A.rds')
panc.b <- readRDS('~/sc_eval/data/v2/PancreaticMaturation_B.rds')

tki_depth <- sapply(tki, function(x){
  ncol(counts(x,'raw'))
})
tki_depth <- data.table('data' = rep(paste0('data.', 1:12), 3),
                        'depth' = tki_depth,
                        'tag' = c(rep('Alectinib', 12), rep('Lorlatinib', 12), rep('Crizotinib', 12)))
neur_depth <- data.table('data' = paste0('data.', 1:12),
                         'depth' = sapply(neur, function(x){ncol(counts(x,'raw'))}),
                         'tag'= rep('Neurodegeneration', 12))
mcf7_depth <- data.table('data' = paste0('data.', 1:12),
                         'depth' = sapply(mcf7, function(x){ncol(counts(x,'raw'))}),
                         'tag'= rep('E2_Treatment_MCF7', 12))
t47d_depth <- data.table('data' = paste0('data.', 1:12),
                         'depth' = sapply(t47d, function(x){ncol(counts(x,'raw'))}),
                         'tag'= rep('E2_Treatment_T47D', 12))
panca_depth <- data.table('data' = paste0('data.', 1:12),
                         'depth' = sapply(panc.a, function(x){ncol(counts(x,'raw'))}),
                         'tag'= rep('PancreaticMaturation_A', 12))
pancb_depth <- data.table('data' = paste0('data.', 1:12),
                          'depth' = sapply(panc.b, function(x){ncol(counts(x,'raw'))}),
                          'tag'= rep('PancreaticMaturation_B', 12))
dept_data <- rbindlist(list(tki_depth, neur_depth, mcf7_depth, t47d_depth, panca_depth, pancb_depth))
ari_stats <- data.table('data' = combined_ari$Data,
                        'ARI' = combined_ari$ARI,
                        'tag' = combined_ari$DataTag)
ari_stats <- ari_stats %>%
  left_join(dept_data, by=c('data', 'tag'))
lm_res <- glm(ARI~log10(depth), data=ari_stats, family = 'gaussian')



p <- ggplot(combined_ari, aes(x=DataTag, y=ARI)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.01, binwidth=1/80) + 
  #geom_jitter(position=position_jitter(0.2), size=0.1, alpha=0.1) +
  ylab('Adjusted Rand Index') +
  xlab('') +
  theme_classic() +
  annotate(geom='text', x=7.5, y=0.75, label=sprintf('%s: %.3f', expr('beta'), 0.033), parse=TRUE, fontface=2, size=4) +
  #facet_wrap(~Data) +
  theme(axis.text.x = element_text(face='bold', size=10, angle = 45, hjust = 1),
        axis.title.y = element_text(face='bold', size=12)) +
  scale_x_discrete(labels=setNames(c('TKI Treatment Alectinib',
                                   'TKI Treatment Lorlatinib',
                                   'TKI Treatment Crizotinib',
                                   'Neurodegeneration Model',
                                   'Pancreatic Differentiation B',
                                   'Pancreatic Differentiation A',
                                   'E2 Treatment T47D',
                                   'E2 Treatment MCF7'),
                                   c('Alectinib', 'Lorlatinib', 'Crizotinib',
                                     'Neurodegeneration', 'PancreaticMaturation_B', 'PancreaticMaturation_A',
                                     'E2_Treatment_T47D', 'E2_Treatment_MCF7')))

ggsave(p, filename = '~/sc_eval/plots/ClusterBoxplot.pdf', width = 8, height = 6)

## Plot Cluster Comparisons Heatmap ##
for(data_tag in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'E2_Treatment_MCF7', 'E2_Treatment_T47D', 'Neurodegeneration', 'PancreaticMaturation_A', 'PancreaticMaturation_B')){
  cluster_comp <- readRDS(sprintf('data/v2/%s_Cluster_Comparisons.rds', data_tag))
  cluster_comp <- sapply(1:length(cluster_comp), simplify = FALSE, function(x.idx){
    x <- cluster_comp[[x.idx]]
    x_ft <- na.omit(reshape2::melt(x))
    colnames(x_ft) <- c('Method.a', 'Method.b', 'ARI')
    x_ft$Data <- rep(paste0('data.', x.idx), nrow(x_ft))
    return(x_ft)
  })
  cluster_comp <- rbindlist(cluster_comp)
  cluster_stats <- cluster_comp %>%
    group_by(Method.a, Method.b) %>%
    summarise(med=median(ARI),
              sd=sd(ARI))
  cluster_stats <- rbindlist(list(cluster_stats,
                                  data.table('Method.a'=cluster_stats$Method.b,
                                             'Method.b'=cluster_stats$Method.a,
                                             'med'=cluster_stats$med,
                                             'sd'=cluster_stats$sd))) %>%
    distinct()
  cluster_stats$Method.a <- factor(cluster_stats$Method.a,
                                   levels = sort(union(cluster_stats$Method.a, cluster_stats$Method.b)))
  cluster_stats$Method.b <- factor(cluster_stats$Method.b,
                                   levels = sort(union(cluster_stats$Method.a, cluster_stats$Method.b)))
  
  med_ft <- reshape2::dcast(data=cluster_stats, Method.a~Method.b, value.var = 'med')
  med_mat <- data.matrix(med_ft[,-1])
  rownames(med_mat) <- med_ft$Method.a
  med_mat <- forceSymmetric(med_mat)
  
  sd_ft <- reshape2::dcast(data=cluster_stats, Method.a~Method.b, value.var = 'sd')
  sd_mat <- data.matrix(sd_ft[,-1])
  rownames(sd_mat) <- sd_ft$Method.a
  
  combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
  combined_stats[upper.tri(combined_stats)] <- med_mat[upper.tri(med_mat)]
  combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
  sd_max <- max(combined_stats[lower.tri(combined_stats)])
  sd_min <- 0
  ari_col_fun <- colorRamp2(breaks=seq(0,1,0.001), colors=viridis(length(seq(0,1,0.001))))
  sd_col_fun <- colorRamp2(breaks=seq(0,sd_max,0.001), colors=magma(length(seq(0,sd_max,0.001))))
  
  cell_fun <- function(j, i, x, y, width, height, fill){
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = NA))
    if(i == j) {
      grid.text('', x = x, y = y)
    } else if(i > j) {
      grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
                  gp = gpar(fill = sd_col_fun(combined_stats[i, j]), col = NA))
    } else {
      grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = "grey", fill = ari_col_fun(combined_stats[i, j])))
    }
  }
  
  ## Annotations ##
  Data <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){s[1]})
  Normalization <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
  DimensionRed <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[4:length(s)], collapse = '_'), paste(s[3:length(s)], collapse = '_'))})
  
  data_image_paths <- ifelse(Data == 'raw', '',
                             ifelse(Data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
  norm_image_paths <- ifelse(Normalization == 'deconv', 'plots/deconv_icon.png',
                             ifelse(Normalization == 'sctransform', 'plots/sctransform_icon.png',
                                    ifelse(Normalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
  red_image_paths <- ifelse(DimensionRed == 'umap', 'plots/umap_icon.png',
                            ifelse(DimensionRed == 'paga_umap', 'plots/paga_umap_icon.png',
                                   ifelse(DimensionRed == 'vae', 'plots/vae_icon.png',
                                          ifelse(DimensionRed == 'tsne', 'plots/tsne_icon.png', 'plots/dm_icon.png'))))
  
  ha = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                         Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                         DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                         height = unit(1.2, 'cm'),
                         show_annotation_name = FALSE,
                         show_legend = FALSE)
  
  hr = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                         Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                         DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                         width = unit(1.2, 'cm'),
                         show_annotation_name = FALSE,
                         show_legend = FALSE, 
                         which = 'row')
  
  
  h <- Heatmap(combined_stats,
               cell_fun = cell_fun,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE,
               show_row_names = FALSE,
               column_names_gp = gpar(fontsize=10),
               row_names_gp = gpar(fontsize=10),
               column_names_rot = 60,
               show_heatmap_legend = FALSE,
               top_annotation = ha,
               left_annotation = hr,
               #column_title=sprintf('%s', data_tag),
               #column_title_gp = gpar(fontsize=16, fontface='bold'),
               rect_gp = gpar(type = "none"),
               column_names_side = 'top',
               na_col = 'white')
  
  ari_leg <- Legend(col_fun = ari_col_fun, 
                    title='Adjusted Rand Index', 
                    title_position='leftcenter-rot', 
                    legend_height=unit(6,'cm'),
                    at = c(0, 0.25, 0.5, 0.75, 1))
  
  sd_leg <- Legend(col_fun = sd_col_fun,
                   title = 'Variation',
                   title_position = 'leftcenter-rot',
                   legend_height=unit(6,'cm'))
  
  pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))
  
  pdf(sprintf('plots/ClusterComparison_%s.pdf', data_tag), width = 10, height = 8)
  draw(h, annotation_legend_list=pd)
  dev.off()
}

## Plot WOT ##
combined_wot <- list()
for(title_sub in c('Alectinib', 'Crizotinib', 'Lorlatinib', 'GSE103334', 'GSE87375_A', 'GSE87375_B', 'GSE107858', 'GSE107863')){
  for(i in c(1,2,3,4)){
    for(j in c(1,2,3)){
      wot_cor <- readRDS(sprintf("~/sc_eval/data/v2/%s_MADS%d_LevS%d/%s_MADS%d_LevS%d_Correlation_WOT.rds", title_sub, i, j, title_sub, i, j))
      wot_cor_ft <- sapply(wot_cor, simplify = FALSE, function(x){
        x[lower.tri(x)] <- NA
        return(reshape2::melt(x))
      })
      wot_cor_ft <- na.omit(do.call('rbind', wot_cor_ft))
      colnames(wot_cor_ft) <- c('Method.A', 'Method.B', 'Cor')
      wot_cor_ft <- wot_cor_ft[wot_cor_ft$Method.A != wot_cor_ft$Method.B,]
      wot_cor_ft <- rbind(data.frame('Method.A' = wot_cor_ft$Method.A,
                                     'Method.B' = wot_cor_ft$Method.B,
                                     'Cor' = wot_cor_ft$Cor,
                                     stringsAsFactors = FALSE),
                          data.frame('Method.A' = wot_cor_ft$Method.B,
                                     'Method.B' = wot_cor_ft$Method.A,
                                     'Cor' = wot_cor_ft$Cor,
                                     stringsAsFactors = FALSE))
      wot_cor_ft$MAD <- rep(paste0('MADS', i), nrow(wot_cor_ft))
      wot_cor_ft$Lev <- rep(paste0('LevS', j), nrow(wot_cor_ft))
      wot_cor_ft$Data <- rep(title_sub, nrow(wot_cor_ft))
      combined_wot[[length(combined_wot)+1]] <- wot_cor_ft
      
      # p <- ggplot(wot_cor_ft, aes(x=Method.A, y=Cor, fill=Method.B)) +
      #   geom_boxplot(position = position_dodge2()) +
      #   xlab('Normalization Method A') +
      #   ylab('Spearman Correlation') +
      #   geom_hline(yintercept = 0, linetype='dashed') +
      #   ylim(-1,1) +
      #   theme_classic() +
      #   theme(axis.text = element_text(size=8),
      #         axis.title = element_text(size=12, face='bold'),
      #         axis.text.x = element_text(angle = 60, hjust = 1)) +
      #   scale_fill_brewer(palette='Dark2') +
      #   guides(fill=guide_legend(title='Normalization Method B',
      #                            title.position = 'left',
      #                            title.theme = element_text(size=12, face='bold', angle = 90)))
      # ggsave(plot=p, filename=sprintf('~/sc_eval/plots/%s_MADS%d_LevS%d_WOT_Correlation.pdf', title_sub, i, j), width = 8, height = 6)
    }
  }
}
combined_wot <- rbindlist(combined_wot)
combined_wot <- rbindlist(list(combined_wot,
                               data.table('Method.A' = combined_wot$Method.B,
                                          'Method.B' = combined_wot$Method.A,
                                          'Cor' = combined_wot$Cor,
                                          'MAD' = combined_wot$MAD,
                                          'Lev' = combined_wot$Lev,
                                          'Data' = combined_wot$Data)))
combined_wot$Data <- factor(combined_wot$Data, levels = c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE103334', 'GSE87375_A', 'GSE87375_B', 'GSE107858', 'GSE107863'))
man_colors <- setNames(RColorBrewer::brewer.pal(n=8, 'Dark2'), levels(combined_wot$Data))
all.methods <- c('raw_deconv', 'raw_sctransform', 'raw_dca_nb', 'raw_dca_zinb', 'scimpute_deconv', 'scimpute_sctransform', 'drimpute_deconv', 'drimpute_sctransform')

## Global Boxplot ##
custom_labels <- function(s_in){
  sapply(s_in, function(s){
    if(s == 'raw_deconv'){
      return('Raw')
    }else if(s == 'raw_sctransform'){
      return('Raw')
    }else if(s == 'raw_dca_nb'){
      return('Raw')
    }else if(s == 'raw_dca_zinb'){
      return('Raw')
    }else if(s == 'scimpute_deconv'){
      return('ScImpute')
    }else if(s == 'scimpute_sctransform'){
      return('ScImpute')
    }else if(s == 'drimpute_deconv'){
      return('DrImpute')
    }else if(s == 'drimpute_sctransform'){
      return('DrImpute')
    }
  })
}

p <- ggplot(combined_wot, aes(x=Data, y=Cor, fill=Data)) +
  geom_boxplot(position = position_dodge2()) +
  theme_classic() +
  theme(axis.text = element_text(face='bold', size=10),
        strip.text = element_text(face='bold', size=10),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text.align = 0.5,
        axis.title= element_text(face='bold', size=12),
        panel.spacing = unit(1.3, "lines"),
        legend.text = element_text(size=10)) +
  ylim(0,1) +
  ylab('Spearman Correlation') +
  xlab('') +
  facet_grid(Method.A~Method.B, labeller = labeller(Method.B=custom_labels,Method.A=custom_labels)) +
  scale_fill_manual(values = man_colors,
                    breaks = names(man_colors),
                    labels = c(expression('NSCLC Alectinib\nTreatment'), 
                               expression('NSCLC Lorlatinib\nTreatment'), 
                               expression('NSCLC Crizotinib\nTreatment'),
                               'Neurodegeneration',
                               expression(paste('Maturation ',alpha)),
                               expression(paste('Maturation ',beta)),
                               expression('MCF7 E2\nTreatment'),
                               expression('T47D E2\nTreatment'))) +
  guides(fill=guide_legend(title=''))
ggsave(p, filename = sprintf('plots/wot/Boxplot_WOT_Stats.pdf'), width = 18, height = 16)

# Compare ScTransform vs Deconvolution given different imputations #
get_uniq_comb <- function(s=NULL,k=NULL){
  s <- s[grepl(pattern=k, s)]
  local_l <- t(combn(s, 2))
  return(data.table('Method.A'=local_l[,1],
                    'Method.B'=local_l[,2]))
}

require(cowplot)
filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'sctransform')
combined_wot_sct <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'deconv')
combined_wot_dec <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
combined_wot_sct$Label <- paste0(custom_labels(combined_wot_sct$Method.A), ",", custom_labels(combined_wot_sct$Method.B))
combined_wot_sct$Type <- 'ScTransform'
combined_wot_dec$Label <- paste0(custom_labels(combined_wot_dec$Method.A), ",", custom_labels(combined_wot_dec$Method.B))
combined_wot_dec$Type <- 'Deconvolution'
temp <- rbindlist(list(combined_wot_sct, combined_wot_dec))

p <- ggplot(temp, aes(x=Data, y=Cor, fill=Data)) +
  geom_boxplot(position = position_dodge2()) +
  theme_minimal() +
  theme(axis.text = element_text(face='bold', size=10),
        strip.text.x = element_text(face='bold', size=12),
        strip.text.y = element_text(face='bold', size=12),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text.align = 0.5,
        axis.title= element_text(face='bold', size=12),
        panel.spacing = unit(1.3, "lines"),
        legend.text = element_text(size=10)) +
  ylim(0,1) +
  ylab('Spearman Correlation') +
  xlab('') +
  facet_grid(Type~Label) +
  scale_fill_manual(values = man_colors,
                    breaks = names(man_colors),
                    labels = c(expression('TKI Treatment Alectinib\n'), 
                               expression('TKI Treatment Lorlatinib\n'), 
                               expression('TKI Treatment Crizotinib\n'),
                               'Neurodegeneration',
                               expression(paste('Pancreatic Maturation ',alpha)),
                               expression(paste('Pancratic Maturation ',beta)),
                               expression('E2 Treatment\nMCF7'),
                               expression('E2 Treatment\nT47D'))) +
  guides(fill=guide_legend(title=''))
ggsave(p, filename = sprintf('plots/wot/WOT_Norm.pdf'), width = 10, height = 5)

# Compare DrImpute, ScImpute and Raw given different normalizations #
custom_labels <- function(s_in){
  sapply(s_in, function(s){
    if(s == 'raw_deconv'){
      return('Deconvolution')
    }else if(s == 'raw_sctransform'){
      return('ScTransform')
    }else if(s == 'raw_dca_nb'){
      return('DCA_NB')
    }else if(s == 'raw_dca_zinb'){
      return('DCA_ZINB')
    }else if(s == 'scimpute_deconv'){
      return('Deconvolution')
    }else if(s == 'scimpute_sctransform'){
      return('ScTransform')
    }else if(s == 'drimpute_deconv'){
      return('Deconvolution')
    }else if(s == 'drimpute_sctransform'){
      return('ScTransform')
    }
  })
}

filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'raw')
filt_list <- filt_list[!grepl(pattern='dca', filt_list$Method.A) & !grepl(pattern='dca', filt_list$Method.B),]
combined_wot_raw <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'drimpute')
combined_wot_drimpute <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'scimpute')
combined_wot_scimpute <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
combined_wot_raw$Type <- 'Raw'
combined_wot_drimpute$Type <- 'DrImpute'
combined_wot_scimpute$Type <- 'ScImpute'
combined_wot_raw$Label <- paste0(custom_labels(combined_wot_raw$Method.A), ",", custom_labels(combined_wot_raw$Method.B))
combined_wot_drimpute$Label <- paste0(custom_labels(combined_wot_drimpute$Method.A), ",", custom_labels(combined_wot_drimpute$Method.B))
combined_wot_scimpute$Label <- paste0(custom_labels(combined_wot_scimpute$Method.A), ",", custom_labels(combined_wot_scimpute$Method.B))

temp <- rbindlist(list(combined_wot_raw, combined_wot_drimpute, combined_wot_scimpute))

p <- ggplot(temp, aes(x=Data, y=Cor, fill=Data)) +
  geom_boxplot(position = position_dodge2()) +
  theme_minimal() +
  theme(axis.text = element_text(face='bold', size=10),
        strip.text.x = element_text(face='bold', size=12),
        strip.text.y = element_text(face='bold', size=12),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text.align = 0.5,
        axis.title= element_text(face='bold', size=12),
        panel.spacing = unit(1.3, "lines"),
        legend.text = element_text(size=10)) +
  ylim(0,1) +
  ylab('Spearman Correlation') +
  xlab('') +
  facet_grid(Type~Label) +
  scale_fill_manual(values = man_colors,
                    breaks = names(man_colors),
                    labels = c(expression('TKI Treatment Alectinib\n'), 
                               expression('TKI Treatment Lorlatinib\n'), 
                               expression('TKI Treatment Crizotinib\n'),
                               'Neurodegeneration',
                               expression(paste('Pancreatic Maturation ',alpha)),
                               expression(paste('Pancratic Maturation ',beta)),
                               expression('E2 Treatment\nMCF7'),
                               expression('E2 Treatment\nT47D'))) +
  guides(fill=guide_legend(title=''))
ggsave(p, filename = sprintf('plots/wot/WOT_Impute.pdf'), width = 7, height = 6)

# Plot Effect of preprocessing #
get_uniq_comb <- function(s=NULL,k=NULL){
  s <- s[grepl(pattern=k, s)]
  local_l <- t(combn(s, 2))
  return(data.table('Method.A'=local_l[,1],
                    'Method.B'=local_l[,2]))
}

require(cowplot)
filt_list <- get_uniq_comb(s=sort(union(combined_wot$Method.A, combined_wot$Method.B)), k = 'sctransform')
combined_wot_sct <- subset(combined_wot, combined_wot$Method.A %in% filt_list$Method.A & combined_wot$Method.B %in% filt_list$Method.B)
combined_wot_sct <- combined_wot_sct %>%
  group_by(Method.A, Method.B, MAD, Lev, Data) %>%
  summarise(CorMed = mean(Cor))
combined_wot_sct$Label <- paste0(custom_labels(combined_wot_sct$Method.A), ",", custom_labels(combined_wot_sct$Method.B))
combined_wot_sct$Type <- 'ScTransform'
combined_wot_sct$MADN <- as.numeric(gsub(pattern='^MADS', replacement = '', combined_wot_sct$MAD))

level_labs <- function(s){
  if(s == 'LevS1'){
    return("1.0%")
  }else if(s == 'LevS2'){
    return("5.0%")
  }else if(s == 'LevS3'){
    return("10.0%")
  }
}

p <- ggplot(combined_wot_sct, aes(y=CorMed, x=MADN, color=Data)) +
  geom_line() +
  theme_minimal() +
  theme(axis.text = element_text(face='bold', size=10),
        strip.text.x = element_text(face='bold', size=12),
        strip.text.y = element_text(face='bold', size=12),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text.align = 0.5,
        axis.title= element_text(face='bold', size=12),
        panel.spacing = unit(1.3, "lines"),
        legend.text = element_text(size=10)) +
  ylim(0,1) +
  ylab('Spearman Correlation') +
  xlab('Cell Level Filter') +
  facet_grid(Label~Lev, labeller = labeller(Lev=setNames(c('1.0%', '5.0%', '10.0%'),c('LevS1', 'LevS2', 'LevS3')))) +
  scale_color_manual(values = man_colors,
                    breaks = names(man_colors),
                    labels = c(expression('TKI Treatment Alectinib\n'), 
                               expression('TKI Treatment Lorlatinib\n'), 
                               expression('TKI Treatment Crizotinib\n'),
                               'Neurodegeneration',
                               expression(paste('Pancreatic Maturation ',alpha)),
                               expression(paste('Pancratic Maturation ',beta)),
                               expression('E2 Treatment\nMCF7'),
                               expression('E2 Treatment\nT47D'))) +
  scale_x_reverse() +
  guides(color=guide_legend(title=''))
ggsave(p, filename = 'plots/wot/WOT_Filter.pdf', width = 10, height = 6)

#### WOT Heatmap ####
combined_wot <- list()
for(title_sub in c('Alectinib', 'Crizotinib', 'Lorlatinib', 'GSE103334', 'GSE87375_A', 'GSE87375_B', 'GSE107858', 'GSE107863')){
  for(i in c(1,2,3,4)){
    for(j in c(1,2,3)){
      wot_cor <- readRDS(sprintf("~/sc_eval/data/v2/%s_MADS%d_LevS%d/%s_MADS%d_LevS%d_Correlation_WOT.rds", title_sub, i, j, title_sub, i, j))
      wot_cor_ft <- sapply(wot_cor, simplify = FALSE, function(x){
        x[lower.tri(x)] <- NA
        return(reshape2::melt(x))
      })
      wot_cor_ft <- na.omit(do.call('rbind', wot_cor_ft))
      colnames(wot_cor_ft) <- c('Method.A', 'Method.B', 'Cor')
      wot_cor_ft <- wot_cor_ft[wot_cor_ft$Method.A != wot_cor_ft$Method.B,]
      wot_cor_ft <- rbind(data.frame('Method.A' = wot_cor_ft$Method.A,
                                     'Method.B' = wot_cor_ft$Method.B,
                                     'Cor' = wot_cor_ft$Cor,
                                     stringsAsFactors = FALSE),
                          data.frame('Method.A' = wot_cor_ft$Method.B,
                                     'Method.B' = wot_cor_ft$Method.A,
                                     'Cor' = wot_cor_ft$Cor,
                                     stringsAsFactors = FALSE))
      wot_cor_ft$MAD <- rep(paste0('MADS', i), nrow(wot_cor_ft))
      wot_cor_ft$Lev <- rep(paste0('LevS', j), nrow(wot_cor_ft))
      wot_cor_ft$Data <- rep(title_sub, nrow(wot_cor_ft))
      combined_wot[[length(combined_wot)+1]] <- wot_cor_ft
    }
  }
}
combined_wot <- rbindlist(combined_wot)
combined_wot <- rbindlist(list(combined_wot,
                               data.table('Method.A' = combined_wot$Method.B,
                                          'Method.B' = combined_wot$Method.A,
                                          'Cor' = combined_wot$Cor,
                                          'MAD' = combined_wot$MAD,
                                          'Lev' = combined_wot$Lev,
                                          'Data' = combined_wot$Data)))
combined_wot$Data <- factor(combined_wot$Data, levels = c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE103334', 'GSE87375_A', 'GSE87375_B', 'GSE107858', 'GSE107863'))
man_colors <- setNames(RColorBrewer::brewer.pal(n=8, 'Dark2'), levels(combined_wot$Data))
all.methods <- c('raw_deconv', 'raw_sctransform', 'raw_dca_nb', 'raw_dca_zinb', 'scimpute_deconv', 'scimpute_sctransform', 'drimpute_deconv', 'drimpute_sctransform')
perm.methods <- combn(all.methods, 2)
perm.methods <- data.table('Method.A' = perm.methods[1,],
                           'Method.B' = perm.methods[2,])
combined_wot <- combined_wot[combined_wot$Method.A %in% perm.methods$Method.A & combined_wot$Method.B %in% perm.methods$Method.B,]

ht_list <- sapply(unique(combined_wot$Data), simplify = FALSE, function(s){
  local_combined_wot <- combined_wot %>%
    filter(Data == s)
  
  cluster_stats <- local_combined_wot %>%
    distinct() %>%
    group_by(Method.A, Method.B) %>%
    summarise(med=median(Cor),
              sd=var(Cor))
  cluster_stats <- rbindlist(list(cluster_stats,
                                  data.table('Method.A'=cluster_stats$Method.B,
                                             'Method.B'=cluster_stats$Method.A,
                                             'med'=cluster_stats$med,
                                             'sd'=cluster_stats$sd))) %>%
    distinct()
  cluster_stats$Method.A <- factor(cluster_stats$Method.A,
                                   levels = sort(union(cluster_stats$Method.A, cluster_stats$Method.B)))
  cluster_stats$Method.B <- factor(cluster_stats$Method.B,
                                   levels = sort(union(cluster_stats$Method.A, cluster_stats$Method.B)))
  med_ft <- reshape2::dcast(data=cluster_stats, Method.A~Method.B, value.var = 'med')
  med_mat <- data.matrix(med_ft[,-1])
  rownames(med_mat) <- as.character(med_ft[,1])
  sd_ft <- reshape2::dcast(data=cluster_stats, Method.A~Method.B, value.var = 'sd')
  sd_mat <- data.matrix(sd_ft[,-1])
  rownames(sd_mat) <- sd_ft[,1]
  
  combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
  combined_stats[upper.tri(combined_stats)] <- med_mat[upper.tri(med_mat)]
  combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
  return(combined_stats)
})
ht_top <- do.call('cbind', ht_list[1:4])
ht_bottom <- do.call('cbind', ht_list[5:8])

sd_max <- max(sapply(ht_list, function(x){max(x[lower.tri(x)])}))
sd_min <- 0
ari_col_fun <- colorRamp2(breaks=seq(0,1,0.001), colors=viridis(length(seq(0,1,0.001))))
sd_col_fun <- colorRamp2(breaks=seq(0,0.03,0.001), colors=magma(length(seq(0,0.03,0.001))))

top_cell_fun <- function(j, i, x, y, width, height, fill){
  grid.rect(x = x, y = y, width = width, height = height, 
            gp = gpar(col = "grey", fill = NA))
  temp_idx <- j %% 8
  if(temp_idx == 0)
    temp_idx <- 8
  
  if(i == temp_idx) {
    grid.text('', x = x, y = y)
  } else if(i > temp_idx) {
    grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
              gp = gpar(fill = sd_col_fun(ht_top[i, j]), col = NA))
  } else {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = ari_col_fun(ht_top[i, j])))
  }
}

bottom_cell_fun <- function(j, i, x, y, width, height, fill){
  grid.rect(x = x, y = y, width = width, height = height, 
            gp = gpar(col = "grey", fill = NA))
  temp_idx <- j %% 8
  if(temp_idx == 0)
    temp_idx <- 8
  
  if(i == temp_idx) {
    grid.text('', x = x, y = y)
  } else if(i > temp_idx) {
    grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
              gp = gpar(fill = sd_col_fun(ht_bottom[i, j]), col = NA))
  } else {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = ari_col_fun(ht_bottom[i, j])))
  }
}


#### Annotations ####
RowData <- sapply(stringi::stri_split(rownames(ht_top), fixed = '_'), function(s){s[1]})
RowNormalization <- sapply(stringi::stri_split(rownames(ht_top), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
ColData <- sapply(stringi::stri_split(colnames(ht_top), fixed = '_'), function(s){s[1]})
ColNormalization <- sapply(stringi::stri_split(colnames(ht_top), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})

row_data_image_paths <- ifelse(RowData == 'raw', '',
                           ifelse(RowData == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
row_norm_image_paths <- ifelse(RowNormalization == 'deconv', 'plots/deconv_icon.png',
                           ifelse(RowNormalization == 'sctransform', 'plots/sctransform_icon.png',
                                  ifelse(RowNormalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
col_data_image_paths <- ifelse(ColData == 'raw', '',
                               ifelse(ColData == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
col_norm_image_paths <- ifelse(ColNormalization == 'deconv', 'plots/deconv_icon.png',
                               ifelse(ColNormalization == 'sctransform', 'plots/sctransform_icon.png',
                                      ifelse(ColNormalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))

ha = HeatmapAnnotation(Data=anno_image(col_data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                       Normalization=anno_image(col_norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                       height = unit(1.2, 'cm'),
                       show_annotation_name = FALSE,
                       show_legend = FALSE)

hr = HeatmapAnnotation(Data=anno_image(row_data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                       Normalization=anno_image(row_norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                       width = unit(1.2, 'cm'),
                       show_annotation_name = FALSE,
                       show_legend = FALSE, 
                       which = 'row')


h_top <- Heatmap(ht_top,
             cell_fun = top_cell_fun,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_column_names = FALSE,
             show_row_names = FALSE,
             show_heatmap_legend = FALSE,
             top_annotation = ha,
             left_annotation = hr,
             column_split = do.call('c', sapply(unique(combined_wot$Data)[1:4], simplify = FALSE, function(s){rep(s,8)})),
             rect_gp = gpar(type = "none"),
             column_title = c('NSCLC Alectinib Treatment', 
                              'NSCLC Lorlatinib Treatment', 
                              'NSCLC Crizotinib Treatment',
                              'Neurodegeneration',
                              'Pancreatic Maturation alpha',
                              'Pancreatic Maturation beta',
                              'MCF7 E2 Treatment',
                              'T47D E2 Treatment')[1:4],
             column_title_gp = gpar(fontface='bold'),
             column_gap = unit(0.75, 'cm'),
             na_col = 'white')

h_bott <- Heatmap(ht_bottom,
                 cell_fun = bottom_cell_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 show_heatmap_legend = FALSE,
                 top_annotation = ha,
                 left_annotation = hr,
                 column_split = do.call('c', sapply(unique(combined_wot$Data)[5:8], simplify = FALSE, function(s){rep(s,8)})),
                 rect_gp = gpar(type = "none"),
                 column_title = c('NSCLC Alectinib Treatment', 
                                  'NSCLC Lorlatinib Treatment', 
                                  'NSCLC Crizotinib Treatment',
                                  'Neurodegeneration',
                                  'Pancreatic Maturation alpha',
                                  'Pancreatic Maturation beta',
                                  'MCF7 E2 Treatment',
                                  'T47D E2 Treatment')[5:8],
                 column_title_gp = gpar(fontface='bold'),
                 column_gap = unit(0.75, 'cm'),
                 na_col = 'white')

require(cowplot)
require(grid)
ari_leg <- Legend(col_fun = ari_col_fun, 
                  title='Spearman Correlation', 
                  title_position='leftcenter-rot', 
                  legend_height=unit(6,'cm'),
                  at = c(0, 0.25, 0.5, 0.75, 1))

sd_leg <- Legend(col_fun = sd_col_fun,
                 title = 'Variation',
                 title_position = 'leftcenter-rot',
                 legend_height=unit(6,'cm'))
pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))

p1 <- grid.grabExpr(draw(h_top))
p2 <- grid.grabExpr(draw(h_bott))
p3 <- grid.grabExpr(draw(pd))

p_combined <- plot_grid(p1, p2, ncol = 1, align = "hv")
p_combined_v <- plot_grid(p_combined, p3, ncol=2, rel_widths = c(0.9,0.1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
save_plot(plot=p_combined_v, filename = sprintf('plots/wot/Heatmap_WOT_V2.pdf'), base_width = 10, base_height = 6)

h_combined <- h_bott %v% h_top


pdf(sprintf('plots/wot/Heatmap_WOT_V2.pdf'), width = 12, height = 8)
draw(h_combined, annotation_legend_list=pd, ht_gap = unit(2, "cm"))
dev.off()

#### With Icons (keep just in case) ####
require(cowplot)
all.methods <- c('raw_deconv', 'raw_sctransform', 'raw_dca_nb', 'raw_dca_zinb', 'scimpute_deconv', 'scimpute_sctransform', 'drimpute_deconv', 'drimpute_sctransform')
for(method.a.idx in 1:(length(all.methods)-1)){
  method.a <- all.methods[method.a.idx]
  combined_plots <- list()
  for(method.b.idx in (method.a.idx+1):length(all.methods)){
    method.b <- all.methods[method.b.idx]
    plot_ft <- combined_wot %>%
      filter(Method.A == method.a & Method.B == method.b) %>%
      droplevels()
    
    p <- ggplot(plot_ft, aes(x=MAD, y=Cor, fill=Data)) +
      geom_boxplot(position = position_dodge2()) +
      theme_classic() +
      theme(axis.text = element_text(face='bold', size=10),
            strip.text.x = element_text(face='bold', size=10),
            strip.background = element_blank(),
            legend.position = 'none') +
      ylim(0,1) +
      ylab('Spearman Correlation') +
      xlab('') +
      facet_wrap(~Lev) +
      scale_fill_manual(values = man_colors,
                        breaks = names(man_colors),
                        labels = c('NSCLC Alectinib Treatment', 
                                   'NSCLC Lorlatinib Treatment', 
                                   'NSCLC Crizotinib Treatment',
                                   'Neurodegeneration',
                                   'Pancreatic Maturation alpha',
                                   'Pancreatic Maturation beta',
                                   'MCF7 E2 Treatment',
                                   'T47D E2 Treatment')) +
      guides(fill=guide_legend(title=''))
    
    ## Add icon
    l_data <- sapply(stringi::stri_split(method.b, fixed = '_'), function(s){s[1]})
    l_norm <- sapply(stringi::stri_split(method.b, fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
    
    data_image_paths <- ifelse(l_data == 'raw', '',
                               ifelse(l_data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
    norm_image_paths <- ifelse(l_norm == 'deconv', 'plots/deconv_icon.png',
                               ifelse(l_norm == 'sctransform', 'plots/sctransform_icon.png',
                                      ifelse(l_norm == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
    
    if(l_data != 'raw'){
      p.t <- ggdraw() + 
        draw_plot(p) +
        draw_image(data_image_paths, x = 0.94, y = 0.95, hjust = 1, vjust = 1, halign = 1, valign = 1,
                   width = 0.025) +
        draw_image(norm_image_paths, x = 0.99, y = 0.95, hjust = 1, vjust = 1, halign = 1, valign = 1,
                   width = 0.025)
    }else{
      p.t <- ggdraw() + 
        draw_plot(p) +
        draw_image(norm_image_paths, x = 0.99, y = 0.95, hjust = 1, vjust = 1, halign = 1, valign = 1,
                   width = 0.025)
    }
    combined_plots[[length(combined_plots)+1]] <- p.t
  }
  
  # Get Legend #
  plot_ft <- combined_wot %>%
    filter(Method.A == 'raw_deconv' & Method.B == 'raw_sctransform') %>%
    droplevels()
  
  p <- ggplot(plot_ft, aes(x=MAD, y=Cor, fill=Data)) +
    geom_boxplot(position = position_dodge2()) +
    theme_classic() +
    theme(axis.text = element_text(face='bold', size=10),
          strip.text.x = element_text(face='bold', size=10),
          strip.background = element_blank()) +
    ylim(0,1) +
    ylab('Spearman Correlation') +
    xlab('') +
    facet_wrap(~Lev) +
    scale_fill_manual(values = man_colors,
                      breaks = names(man_colors),
                      labels = c('NSCLC Alectinib Treatment', 
                                 'NSCLC Lorlatinib Treatment', 
                                 'NSCLC Crizotinib Treatment',
                                 'Neurodegeneration',
                                 'Pancreatic Maturation alpha',
                                 'Pancreatic Maturation beta',
                                 'MCF7 E2 Treatment',
                                 'T47D E2 Treatment')) +
    guides(fill=guide_legend(title=''))
  
  legend <- get_legend(
    # create some space to the left of the legend
    p + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  prow <- plot_grid(plotlist = combined_plots, ncol=1, rel_widths = 18, rel_heights = 16)
  lprow <- plot_grid(prow, legend, rel_widths = c(5, 1))
  save_plot(lprow, file=sprintf('~/sc_eval/plots/CombinedPlots.pdf'), base_height = 18, base_width = 12)
}

## Slingshot Calculate Stats ##
local_ent <- function(SpearmanCor=NULL, Jacc=NULL, Traj.a=NULL, Traj.b=NULL){
  df <- data.table('SpearmanCor' = SpearmanCor,
                   'Jacc' = Jacc,
                   'Traj.a' = Traj.a,
                   'Traj.b'= Traj.b)
  df_mat_cor <- data.matrix(reshape2::dcast(df, Traj.a~Traj.b, value.var = 'SpearmanCor')[,-1])
  df_mat_jacc <- data.matrix(reshape2::dcast(df, Traj.a~Traj.b, value.var = 'Jacc')[,-1])
  df_mat <- 0.5*(df_mat_cor + df_mat_jacc)
  x <- as.vector(df_mat)
  if(length(x) == 1){
    return(x)
  }else{
    x <- x/sum(x)
    return(-1*sum(x*log2(ifelse(x==0, 1e-12, x)))/log2(length(x)))  
  }
}

local_ent_sep <- function(Variable=NULL, Traj.a=NULL, Traj.b=NULL){
  x <- as.vector(Variable)
  if(length(x) == 1){
    x <- c(x, 1-x)
  }else{
    x <- x/sum(x)
  }
  return(-1*sum(x*log2(x), na.rm = TRUE)/log2(length(x)))
  
  df <- data.table('Val' = Variable,
                   'Traj.a' = Traj.a,
                   'Traj.b'= Traj.b)
  df_mat <- data.matrix(reshape2::dcast(df, Traj.a~Traj.b, value.var = 'Val')[,-1])
  x <- as.vector(df_mat)
  if(length(x) == 1){
    return(1-x)
  }else{
    x <- x/sum(x)
    return(-1*sum(x*log2(ifelse(x==0, 1e-12, x)))/log2(length(x)))  
  }
}

for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE107863', 'GSE107858', 'GSE87375_A', 'GSE87375_B', 'GSE103334')){
  comp_res <- sapply(1:4, simplify = FALSE, function(i){
    local_comp_res <- sapply(1:3, simplify = FALSE, function(j){
      traject_res <- list()
      for(d_type in c('raw', 'scimpute', 'drimpute')){
        for(n_type in c('deconv', 'sctransform', 'dca_nb', 'dca_zinb')){
          for(r_type in c('umap', 'tsne', 'vae', 'dm', 'paga_umap')){
            if((d_type == 'scimpute' || d_type == 'drimpute') && (n_type == 'dca_nb' || n_type == 'dca_zinb'))
              next
            
            tag <- sprintf('~/sc_eval/data/v2/%s_MADS%d_LevS%d/%s_MADS%d_LevS%d_%s_%s_%s_Slingshot.rds', title_sub, i, j, title_sub, i, j, d_type, n_type, r_type)
            tryCatch({
              traj_res <- readRDS(tag)
              traject_res[[length(traject_res)+1]] <- traj_res[[1]]
              names(traject_res)[length(traject_res)] <- sprintf('%s_%s_%s', d_type, n_type, r_type)
            },
            error=function(e){
              message(sprintf('Error when reading..%s', tag))
            })
          }
        }
      }
      combined_comp_res <- list()
      for(m in 1:(length(traject_res)-1)){
        res.1 <- traject_res[[m]]
        res.1 <- subset(res.1, res.1$Id != 'Root')
        ncurv.1 <- length(grep(pattern='curve', colnames(res.1)))
        for(k in m:length(traject_res)){
          res.2 <- traject_res[[k]]
          res.2 <- subset(res.2, res.2$Id != 'Root')
          ncurv.2 <- length(grep(pattern='curve', colnames(res.2)))
          
          rank_cor <- sapply(1:ncurv.1, simplify = FALSE, function(t){
            temp.1 <- na.omit(setNames(res.1[,which(colnames(res.1) == paste0('curve', t))], res.1$Id))
            local_rank_cor <- sapply(1:ncurv.2, simplify = FALSE, function(z){
              temp.2 <- na.omit(setNames(res.2[,which(colnames(res.2) == paste0('curve', z))], res.2$Id))
              common.ids <- intersect(names(temp.1), names(temp.2))
              if(length(common.ids) < 10)
                return(NA)
              
              temp.1.filt <- temp.1[match(common.ids, table=names(temp.1))]
              temp.2.filt <- temp.2[match(common.ids, table=names(temp.2))]
              local_cor_res <- (cor(temp.1.filt, temp.2.filt, method='spearman')+1)/2
              local_jcc_res <- length(common.ids)/length(union(names(temp.1), names(temp.2)))
              return(c(local_cor_res, local_jcc_res, t, z))
            })
            local_rank_cor_ft <- as.data.frame(do.call('rbind', local_rank_cor))
            colnames(local_rank_cor_ft) <- c('SpearmanCor', 'Jacc', 'Traj.a', 'Traj.b')
            local_rank_cor_ft$Data.a <- rep(sprintf('S%d',i), nrow(local_rank_cor_ft))
            local_rank_cor_ft$Data.b <- rep(sprintf('S%d',j), nrow(local_rank_cor_ft))
            local_rank_cor_ft$Traj.a <- paste0('traj.', local_rank_cor_ft$Traj.a)
            local_rank_cor_ft$Traj.b <- paste0('traj.', local_rank_cor_ft$Traj.b)
            return(local_rank_cor_ft)
          })
          rank_cor <- do.call('rbind', rank_cor)
          rank_cor$Method.a <- rep(names(traject_res)[m], nrow(rank_cor))
          rank_cor$Method.b <- rep(names(traject_res)[k], nrow(rank_cor))
          combined_comp_res[[length(combined_comp_res)+1]] <- as.data.table(rank_cor)
        }
      }
      combined_comp_res <- rbindlist(combined_comp_res)
      combined_comp_res <- combined_comp_res %>%
        filter(Method.a != Method.b) %>%
        distinct()
      
      require(dplyr)
      combined_comp_res_filt <- combined_comp_res %>%
        group_by(Method.a, Method.b) %>%
        summarise(SpearmanEnt=local_ent_sep(SpearmanCor, Traj.a, Traj.b),
                  JaccEnt=local_ent_sep(Jacc, Traj.a, Traj.b))
      return(combined_comp_res_filt)
    })
    names(local_comp_res) <- paste0('LevS', 1:3)
    return(local_comp_res)
  })
  names(comp_res) <- paste0('MADS', 1:4)
  saveRDS(comp_res, sprintf('~/sc_eval/data/v2/%s_Slingshot_Comparison_Ent.rds', title_sub))
}



## Slinghost Plot ##
for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE107863', 'GSE107858', 'GSE87375_A', 'GSE87375_B', 'GSE103334')){
  comp_res <- readRDS(sprintf('~/sc_eval/data/v2/%s_Slingshot_Comparison.rds', title_sub))
  comp_res_combined <- sapply(names(comp_res), simplify = FALSE, function(s){
    local_comp_res_combined <- sapply(names(comp_res[[s]]), simplify = FALSE, function(k){
      temp <- comp_res[[s]][[k]]
      temp$Levels <- rep(k, nrow(temp))
      return(temp)
    })
    local_comp_res_combined <- rbindlist(local_comp_res_combined)
    local_comp_res_combined$MADS <- rep(s, nrow(local_comp_res_combined))
    return(local_comp_res_combined)
  })
  comp_res_combined <- rbindlist(comp_res_combined)
  comp_res_combined <- rbindlist(list(comp_res_combined,
                                      data.table('Method.a'=comp_res_combined$Method.b,
                                                 'Method.b'=comp_res_combined$Method.a,
                                                 'Jacc'=comp_res_combined$Jacc,
                                                 'SpearmanCor'=comp_res_combined$SpearmanCor,
                                                 'Levels'=comp_res_combined$Levels,
                                                 'MADS'=comp_res_combined$MADS))) %>%
    distinct()
  p <- ggplot(comp_res_combined, aes(x=Method.a, y=Method.b, fill=SpearmanCor)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6),
          axis.text.y = element_text(size=6)) +
    facet_grid(Levels~MADS) +
    scale_fill_viridis_c(begin = 0, end = 1, values = seq(0,1,0.01))
  
  local_stats <- comp_res_combined %>%
    group_by(Method.a, Method.b) %>%
    summarise(SpearmanCorMed = median(SpearmanCor, na.rm = TRUE),
              SpearmanCorVar = var(SpearmanCor, na.rm = TRUE),
              JaccMed = median(Jacc, na.rm = TRUE),
              JaccVar = var(Jacc, na.rm = TRUE))
  
  local_stats <- rbindlist(list(local_stats,
                                data.table('Method.a'=local_stats$Method.b,
                                           'Method.b'=local_stats$Method.a,
                                           'SpearmanCorMed'=local_stats$SpearmanCorMed,
                                           'SpearmanCorVar'=local_stats$SpearmanCorVar,
                                           'JaccMed'=local_stats$JaccMed,
                                           'JaccVar'=local_stats$JaccVar)))
  local_stats$Method.a <- factor(local_stats$Method.a,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  local_stats$Method.b <- factor(local_stats$Method.b,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  
  for(v in c('SpearmanCorMed', 'JaccMed')){
    if(v == 'SpearmanCorMed'){
      med_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'SpearmanCorMed')
      med_mat <- data.matrix(med_ft[,-1])
      rownames(med_mat) <- med_ft$Method.a
      med_mat <- forceSymmetric(med_mat)
      
      sd_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'SpearmanCorVar')
      sd_mat <- data.matrix(sd_ft[,-1])
      rownames(sd_mat) <- sd_ft$Method.a
    }else{
      med_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'JaccMed')
      med_mat <- data.matrix(med_ft[,-1])
      rownames(med_mat) <- med_ft$Method.a
      med_mat <- forceSymmetric(med_mat)
      
      sd_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'JaccVar')
      sd_mat <- data.matrix(sd_ft[,-1])
      rownames(sd_mat) <- sd_ft$Method.a
    }
    try({
      combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
      combined_stats[upper.tri(combined_stats)] <- med_mat[upper.tri(med_mat)]
      combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
      sd_max <- max(combined_stats[lower.tri(combined_stats)])
      sd_min <- 0
      ari_col_fun <- colorRamp2(breaks=seq(0,1,0.001), colors=viridis(length(seq(0,1,0.001))))
      sd_col_fun <- colorRamp2(breaks=seq(0,sd_max,0.001), colors=magma(length(seq(0,sd_max,0.001))))
      
      cell_fun <- function(j, i, x, y, width, height, fill){
        grid.rect(x = x, y = y, width = width, height = height, 
                  gp = gpar(col = "grey", fill = NA))
        if(i == j) {
          grid.text('', x = x, y = y)
        } else if(i > j) {
          grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
                    gp = gpar(col = 'grey', fill = sd_col_fun(combined_stats[i, j])))
        } else {
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "grey", fill = ari_col_fun(combined_stats[i, j])))
        }
      }
      
      ## Annotations ##
      Data <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){s[1]})
      Normalization <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
      DimensionRed <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[4:length(s)], collapse = '_'), paste(s[3:length(s)], collapse = '_'))})
      
      data_image_paths <- ifelse(Data == 'raw', '',
                                 ifelse(Data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
      norm_image_paths <- ifelse(Normalization == 'deconv', 'plots/deconv_icon.png',
                                 ifelse(Normalization == 'sctransform', 'plots/sctransform_icon.png',
                                        ifelse(Normalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
      red_image_paths <- ifelse(DimensionRed == 'umap', 'plots/umap_icon.png',
                                ifelse(DimensionRed == 'paga_umap', 'plots/paga_umap_icon.png',
                                       ifelse(DimensionRed == 'vae', 'plots/vae_icon.png',
                                              ifelse(DimensionRed == 'tsne', 'plots/tsne_icon.png', 'plots/dm_icon.png'))))
      
      ha = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                             Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                             DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                             height = unit(1.2, 'cm'),
                             show_annotation_name = FALSE,
                             show_legend = FALSE)
      
      hr = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                             Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                             DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                             width = unit(1.2, 'cm'),
                             show_annotation_name = FALSE,
                             show_legend = FALSE, 
                             which = 'row')
      
      
      h <- Heatmap(combined_stats,
                   cell_fun = cell_fun,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   column_names_gp = gpar(fontsize=10),
                   row_names_gp = gpar(fontsize=10),
                   column_names_rot = 60,
                   show_heatmap_legend = FALSE,
                   top_annotation = ha,
                   left_annotation = hr,
                   rect_gp = gpar(type = "none"),
                   column_names_side = 'top',
                   na_col = 'white',
                   use_raster=FALSE)
      
      if(v=='SpearmanCorMed'){
        ari_leg <- Legend(col_fun = ari_col_fun, 
                          title='Spearman Correlation', 
                          title_position='leftcenter-rot', 
                          legend_height=unit(6,'cm'),
                          at = c(0, 0.25, 0.5, 0.75, 1))
      }else{
        ari_leg <- Legend(col_fun = ari_col_fun, 
                          title='Jaccard Similarity', 
                          title_position='leftcenter-rot', 
                          legend_height=unit(6,'cm'),
                          at = c(0, 0.25, 0.5, 0.75, 1))
      }
      
      sd_leg <- Legend(col_fun=sd_col_fun,
                       title = 'Variation',
                       title_position = 'leftcenter-rot',
                       legend_height = unit(6,'cm'))
      
      pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))
      cairo_pdf(sprintf('plots/slingshot/%s_SlingshotStats_%s_Heatmap.pdf', title_sub, v), width = 10, height = 8)
      draw(h, annotation_legend_list=pd)
      dev.off()
    })
  }
}

# Slingshot Plot Entropy #
for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE107863', 'GSE107858', 'GSE87375_A', 'GSE87375_B', 'GSE103334')){
  comp_res <- readRDS(sprintf('~/sc_eval/data/v2/%s_Slingshot_Comparison_Ent.rds', title_sub))
  comp_res_combined <- sapply(names(comp_res), simplify = FALSE, function(s){
    local_comp_res_combined <- sapply(names(comp_res[[s]]), simplify = FALSE, function(k){
      temp <- comp_res[[s]][[k]]
      temp$Levels <- rep(k, nrow(temp))
      return(temp)
    })
    local_comp_res_combined <- rbindlist(local_comp_res_combined)
    local_comp_res_combined$MADS <- rep(s, nrow(local_comp_res_combined))
    return(local_comp_res_combined)
  })
  comp_res_combined <- rbindlist(comp_res_combined)
  comp_res_combined <- rbindlist(list(comp_res_combined,
                                      data.table('Method.a'=comp_res_combined$Method.b,
                                                 'Method.b'=comp_res_combined$Method.a,
                                                 'SpearmanEnt'=comp_res_combined$SpearmanEnt,
                                                 'JaccEnt' = comp_res_combined$JaccEnt,
                                                 'Levels'=comp_res_combined$Levels,
                                                 'MADS'=comp_res_combined$MADS))) %>%
    distinct()
  local_stats <- comp_res_combined %>%
    group_by(Method.a, Method.b) %>%
    summarise(SpearmanEntM = mean(SpearmanEnt, na.rm = TRUE),
              JaccEntM = mean(JaccEnt, na.rm = TRUE),
              SpearmanEntVar = var(SpearmanEnt, na.rm = TRUE),
              JaccEntVar = var(JaccEnt, na.rm=TRUE)) %>%
    distinct() %>%
    filter(Method.a != Method.b)
  
  p <- ggplot(comp_res_combined, aes(log10(1-SpearmanEnt), log10(1-JaccEnt))) +
    geom_point(size=1, alpha=0.3) +
    theme_classic() +
    geom_vline(xintercept = log10(0.7), linetype='dashed') +
    geom_hline(yintercept = log10(0.7), linetype='dashed') +
    xlab('log10(1-Spearman Entropy)') +
    ylab('log10(1- Jaccard Entropy)')
  ggsave(p, filename = sprintf('~/sc_eval/plots/trajectory/slingshot/Scatter_%s.pdf', title_sub), width = 5, height = 5)
  
  comp_res_combined_sig <- comp_res_combined %>%
    filter(1-SpearmanEnt > 0.7 & 1-JaccEnt > 0.7)
  
  local_stats <- rbindlist(list(local_stats,
                                data.table('Method.a'=local_stats$Method.b,
                                           'Method.b'=local_stats$Method.a,
                                           'Med'=local_stats$Med,
                                           'EntVar'=local_stats$EntVar))) %>%
    distinct()
  local_stats$Method.a <- factor(local_stats$Method.a,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  local_stats$Method.b <- factor(local_stats$Method.b,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  
  med_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'Med')
  med_mat <- data.matrix(med_ft[,-1])
  rownames(med_mat) <- med_ft$Method.a
  med_mat <- forceSymmetric(med_mat)
      
  sd_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'EntVar')
  sd_mat <- data.matrix(sd_ft[,-1])
  rownames(sd_mat) <- sd_ft$Method.a
  try({
    combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
    combined_stats[upper.tri(combined_stats)] <- 1-med_mat[upper.tri(med_mat)]
    combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
    sd_max <- max(combined_stats[lower.tri(combined_stats)])
    ari_max <- round(max(1-med_mat, na.rm=TRUE) + 0.001, digits = 2)
    ari_breaks <- pretty(seq(0,ari_max,0.001))
    sd_min <- 0
    ari_col_fun <- colorRamp2(breaks=ari_breaks, colors=viridis(length(ari_breaks)))
    sd_col_fun <- colorRamp2(breaks=seq(0,sd_max,0.001), colors=magma(length(seq(0,sd_max,0.001))))
      
    cell_fun <- function(j, i, x, y, width, height, fill){
      grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = "grey", fill = NA))
      if(i == j) {
        grid.text('', x = x, y = y)
      } else if(i > j) {
        grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
                  gp = gpar(col = 'grey', fill = sd_col_fun(combined_stats[i, j])))
      } else {
        grid.rect(x = x, y = y, width = width, height = height, 
                  gp = gpar(col = "grey", fill = ari_col_fun(combined_stats[i, j])))
      }
    }
      
    ## Annotations ##
    Data <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){s[1]})
    Normalization <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
    DimensionRed <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[4:length(s)], collapse = '_'), paste(s[3:length(s)], collapse = '_'))})
      
    data_image_paths <- ifelse(Data == 'raw', '',
                               ifelse(Data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
    norm_image_paths <- ifelse(Normalization == 'deconv', 'plots/deconv_icon.png',
                               ifelse(Normalization == 'sctransform', 'plots/sctransform_icon.png',
                                      ifelse(Normalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
    red_image_paths <- ifelse(DimensionRed == 'umap', 'plots/umap_icon.png',
                              ifelse(DimensionRed == 'paga_umap', 'plots/paga_umap_icon.png',
                                     ifelse(DimensionRed == 'vae', 'plots/vae_icon.png',
                                            ifelse(DimensionRed == 'tsne', 'plots/tsne_icon.png', 'plots/dm_icon.png'))))
    
    ha = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           height = unit(1.2, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE)
    
    hr = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           width = unit(1.2, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE, 
                           which = 'row')
      
      
    h <- Heatmap(combined_stats,
                 cell_fun = cell_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 column_names_gp = gpar(fontsize=10),
                 row_names_gp = gpar(fontsize=10),
                 column_names_rot = 60,
                 show_heatmap_legend = FALSE,
                 top_annotation = ha,
                 left_annotation = hr,
                 rect_gp = gpar(type = "none"),
                 column_names_side = 'top',
                 na_col = 'white',
                 use_raster=FALSE)
      
    ari_leg <- Legend(col_fun = ari_col_fun, 
                      title='(1-Entropy)', 
                      title_position='leftcenter-rot', 
                      legend_height=unit(6,'cm'))

    sd_leg <- Legend(col_fun=sd_col_fun,
                     title = 'Variation',
                     title_position = 'leftcenter-rot',
                     legend_height = unit(6,'cm'))
      
    pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))
    cairo_pdf(sprintf('plots/slingshot/%s_SlingshotStats_Heatmap_Ent.pdf', title_sub), width = 10, height = 8)
    draw(h, annotation_legend_list=pd)
    dev.off()
  })
}


for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE107863', 'GSE107858', 'GSE87375_A', 'GSE87375_B', 'GSE103334')){
  comp_res <- readRDS(sprintf('~/sc_eval/data/v2/%s_Slingshot_Comparison_Ent.rds', title_sub))
  comp_res_combined <- sapply(names(comp_res), simplify = FALSE, function(s){
    local_comp_res_combined <- sapply(names(comp_res[[s]]), simplify = FALSE, function(k){
      temp <- comp_res[[s]][[k]]
      temp$Levels <- rep(k, nrow(temp))
      return(temp)
    })
    local_comp_res_combined <- rbindlist(local_comp_res_combined)
    local_comp_res_combined$MADS <- rep(s, nrow(local_comp_res_combined))
    return(local_comp_res_combined)
  })
  comp_res_combined <- rbindlist(comp_res_combined)
  comp_res_combined <- rbindlist(list(comp_res_combined,
                                      data.table('Method.a'=comp_res_combined$Method.b,
                                                 'Method.b'=comp_res_combined$Method.a,
                                                 'SpearmanEnt'=comp_res_combined$SpearmanEnt,
                                                 'JaccEnt'=comp_res_combined$JaccEnt,
                                                 'Levels'=comp_res_combined$Levels,
                                                 'MADS'=comp_res_combined$MADS))) %>%
    distinct() %>%
    filter(Method.a != Method.b)
  
  local_stats <- comp_res_combined %>%
    group_by(Method.a, Method.b) %>%
    summarise(Med = median(SpearmanEnt, na.rm = TRUE),
              EntVar = var(SpearmanEnt, na.rm = TRUE))
  
  local_stats <- rbindlist(list(local_stats,
                                data.table('Method.a'=local_stats$Method.b,
                                           'Method.b'=local_stats$Method.a,
                                           'Med'=local_stats$Med,
                                           'EntVar'=local_stats$EntVar))) %>%
    distinct()
  local_stats$Method.a <- factor(local_stats$Method.a,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  local_stats$Method.b <- factor(local_stats$Method.b,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  
  med_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'Med')
  med_mat <- data.matrix(med_ft[,-1])
  rownames(med_mat) <- med_ft$Method.a
  med_mat <- forceSymmetric(med_mat)
  
  sd_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'EntVar')
  sd_mat <- data.matrix(sd_ft[,-1])
  rownames(sd_mat) <- sd_ft$Method.a
  
  try({
    combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
    combined_stats[upper.tri(combined_stats)] <- 1-med_mat[upper.tri(med_mat)]
    combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
    sd_max <- max(combined_stats[lower.tri(combined_stats)])
    ari_max <- round(max(1-med_mat, na.rm=TRUE) + 0.001, digits = 2)
    ari_breaks <- pretty(seq(0,ari_max,0.001))
    sd_min <- 0
    ari_col_fun <- colorRamp2(breaks=ari_breaks, colors=viridis(length(ari_breaks)))
    sd_col_fun <- colorRamp2(breaks=seq(0,sd_max,0.001), colors=magma(length(seq(0,sd_max,0.001))))
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = "grey", fill = NA))
      if(i == j) {
        grid.text('', x = x, y = y)
      } else if(i > j) {
        grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
                  gp = gpar(col = 'grey', fill = sd_col_fun(combined_stats[i, j])))
      } else {
        grid.rect(x = x, y = y, width = width, height = height, 
                  gp = gpar(col = "grey", fill = ari_col_fun(combined_stats[i, j])))
      }
    }
    
    ## Annotations ##
    Data <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){s[1]})
    Normalization <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
    DimensionRed <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[4:length(s)], collapse = '_'), paste(s[3:length(s)], collapse = '_'))})
    
    data_image_paths <- ifelse(Data == 'raw', '',
                               ifelse(Data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
    norm_image_paths <- ifelse(Normalization == 'deconv', 'plots/deconv_icon.png',
                               ifelse(Normalization == 'sctransform', 'plots/sctransform_icon.png',
                                      ifelse(Normalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
    red_image_paths <- ifelse(DimensionRed == 'umap', 'plots/umap_icon.png',
                              ifelse(DimensionRed == 'paga_umap', 'plots/paga_umap_icon.png',
                                     ifelse(DimensionRed == 'vae', 'plots/vae_icon.png',
                                            ifelse(DimensionRed == 'tsne', 'plots/tsne_icon.png', 'plots/dm_icon.png'))))
    
    ha = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           height = unit(1.0, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE)
    
    hr = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           width = unit(1.0, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE, 
                           which = 'row')
    
    
    h <- Heatmap(combined_stats,
                 cell_fun = cell_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 column_names_gp = gpar(fontsize=10),
                 row_names_gp = gpar(fontsize=10),
                 column_names_rot = 60,
                 show_heatmap_legend = FALSE,
                 top_annotation = ha,
                 left_annotation = hr,
                 rect_gp = gpar(type = "none"),
                 column_names_side = 'top',
                 na_col = 'white',
                 use_raster=FALSE)
    
    ari_leg <- Legend(col_fun = ari_col_fun, 
                      title='(1-Entropy)', 
                      title_position='leftcenter-rot', 
                      legend_height=unit(4,'cm'), 
                      legend_gp = gpar(fontsize=12))
    
    sd_leg <- Legend(col_fun=sd_col_fun,
                     title = 'Variation',
                     title_position = 'leftcenter-rot',
                     legend_height = unit(4,'cm'),
                     legend_gp = gpar(fontsize=12))
    
    pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))
    cairo_pdf(sprintf('plots/slingshot/%s_SlingshotStats_Heatmap_Ent.pdf', title_sub), width = 6, height = 5)
    draw(h, annotation_legend_list=pd)
    dev.off()
  })
}


#### Variation-Ent ####
for(title_sub in c('Alectinib', 'Lorlatinib', 'Crizotinib', 'GSE107863', 'GSE107858', 'GSE87375_A', 'GSE87375_B', 'GSE103334')){
  comp_res <- readRDS(sprintf('~/sc_eval/data/v2/%s_Slingshot_Comparison_Ent.rds', title_sub))
  comp_res_combined <- sapply(names(comp_res), simplify = FALSE, function(s){
    local_comp_res_combined <- sapply(names(comp_res[[s]]), simplify = FALSE, function(k){
      temp <- comp_res[[s]][[k]]
      temp$Levels <- rep(k, nrow(temp))
      return(temp)
    })
    local_comp_res_combined <- rbindlist(local_comp_res_combined)
    local_comp_res_combined$MADS <- rep(s, nrow(local_comp_res_combined))
    return(local_comp_res_combined)
  })
  comp_res_combined <- rbindlist(comp_res_combined)
  comp_res_combined <- rbindlist(list(comp_res_combined,
                                      data.table('Method.a'=comp_res_combined$Method.b,
                                                 'Method.b'=comp_res_combined$Method.a,
                                                 'SpearmanEnt'=comp_res_combined$SpearmanEnt,
                                                 'JaccEnt'=comp_res_combined$JaccEnt,
                                                 'Levels'=comp_res_combined$Levels,
                                                 'MADS'=comp_res_combined$MADS))) %>%
    distinct() %>%
    filter(Method.a != Method.b)
  
  local_stats <- comp_res_combined %>%
    group_by(Method.a, Method.b) %>%
    summarise(Med = median(SpearmanEnt, na.rm = TRUE),
              EntVar = var(SpearmanEnt, na.rm = TRUE))
  
  p <- ggplot(local_stats, aes(x=Med, y=EntVar)) +
    geom_point()
  
  local_stats <- rbindlist(list(local_stats,
                                data.table('Method.a'=local_stats$Method.b,
                                           'Method.b'=local_stats$Method.a,
                                           'Med'=local_stats$Med,
                                           'EntVar'=local_stats$EntVar))) %>%
    distinct()
  local_stats$Method.a <- factor(local_stats$Method.a,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  local_stats$Method.b <- factor(local_stats$Method.b,
                                 levels = sort(union(local_stats$Method.a, local_stats$Method.b)))
  
  med_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'Med')
  med_mat <- data.matrix(med_ft[,-1])
  rownames(med_mat) <- med_ft$Method.a
  med_mat <- forceSymmetric(med_mat)
  
  sd_ft <- reshape2::dcast(data=local_stats, Method.a~Method.b, value.var = 'EntVar')
  sd_mat <- data.matrix(sd_ft[,-1])
  rownames(sd_mat) <- sd_ft$Method.a
  
  try({
    combined_stats <- matrix(NA, ncol=ncol(med_mat), nrow=nrow(med_mat), dimnames = list(rownames(med_mat), colnames(med_mat)))
    combined_stats[upper.tri(combined_stats)] <- 1-med_mat[upper.tri(med_mat)]
    combined_stats[lower.tri(combined_stats)] <- sd_mat[lower.tri(sd_mat)]
    sd_max <- max(combined_stats[lower.tri(combined_stats)])
    ari_max <- round(max(1-med_mat, na.rm=TRUE) + 0.001, digits = 2)
    ari_breaks <- pretty(seq(0,ari_max,0.001))
    sd_min <- 0
    ari_col_fun <- colorRamp2(breaks=ari_breaks, colors=viridis(length(ari_breaks)))
    sd_col_fun <- colorRamp2(breaks=seq(0,sd_max,0.001), colors=magma(length(seq(0,sd_max,0.001))))
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = "grey", fill = NA))
      if(i == j) {
        grid.text('', x = x, y = y)
      } else if(i > j) {
        grid.rect(x = x, y = y, width=width, height=height,#abs(combined_stats[i, j])/2 * min(unit.c(width, height)), 
                  gp = gpar(col = 'grey', fill = sd_col_fun(combined_stats[i, j])))
      } else {
        grid.rect(x = x, y = y, width = width, height = height, 
                  gp = gpar(col = "grey", fill = ari_col_fun(combined_stats[i, j])))
      }
    }
    
    ## Annotations ##
    Data <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){s[1]})
    Normalization <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[2:3],collapse='_'), s[2])})
    DimensionRed <- sapply(stringi::stri_split(colnames(combined_stats), fixed = '_'), function(s){ifelse(s[2] == 'dca', paste(s[4:length(s)], collapse = '_'), paste(s[3:length(s)], collapse = '_'))})
    
    data_image_paths <- ifelse(Data == 'raw', '',
                               ifelse(Data == 'scimpute', 'plots/scimpute_icon.png', 'plots/drimpute_icon.png'))
    norm_image_paths <- ifelse(Normalization == 'deconv', 'plots/deconv_icon.png',
                               ifelse(Normalization == 'sctransform', 'plots/sctransform_icon.png',
                                      ifelse(Normalization == 'dca_zinb', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')))
    red_image_paths <- ifelse(DimensionRed == 'umap', 'plots/umap_icon.png',
                              ifelse(DimensionRed == 'paga_umap', 'plots/paga_umap_icon.png',
                                     ifelse(DimensionRed == 'vae', 'plots/vae_icon.png',
                                            ifelse(DimensionRed == 'tsne', 'plots/tsne_icon.png', 'plots/dm_icon.png'))))
    
    ha = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           height = unit(1.0, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE)
    
    hr = HeatmapAnnotation(Data=anno_image(data_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           Normalization=anno_image(norm_image_paths, space=unit(0.1,'mm'), border=FALSE),
                           DimensionReduction=anno_image(red_image_paths, space=unit(0.1, 'mm'), border=FALSE),
                           width = unit(1.0, 'cm'),
                           show_annotation_name = FALSE,
                           show_legend = FALSE, 
                           which = 'row')
    
    
    h <- Heatmap(combined_stats,
                 cell_fun = cell_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 column_names_gp = gpar(fontsize=10),
                 row_names_gp = gpar(fontsize=10),
                 column_names_rot = 60,
                 show_heatmap_legend = FALSE,
                 top_annotation = ha,
                 left_annotation = hr,
                 rect_gp = gpar(type = "none"),
                 column_names_side = 'top',
                 na_col = 'white',
                 use_raster=FALSE)
    
    ari_leg <- Legend(col_fun = ari_col_fun, 
                      title='(1-Entropy)', 
                      title_position='leftcenter-rot', 
                      legend_height=unit(4,'cm'), 
                      legend_gp = gpar(fontsize=12))
    
    sd_leg <- Legend(col_fun=sd_col_fun,
                     title = 'Variation',
                     title_position = 'leftcenter-rot',
                     legend_height = unit(4,'cm'),
                     legend_gp = gpar(fontsize=12))
    
    pd = packLegend(ari_leg, sd_leg, direction = 'vertical', row_gap =  unit(1, "cm"))
    cairo_pdf(sprintf('plots/slingshot/%s_SlingshotStats_Heatmap_Ent.pdf', title_sub), width = 6, height = 5)
    draw(h, annotation_legend_list=pd)
    dev.off()
  })
}


impute_icons <- c('plots/scimpute_icon.png', 'plots/drimpute_icon.png')
norm_icons <- c('plots/deconv_icon.png', 'plots/sctransform_icon.png', 'plots/dca_zinb_icon.png', 'plots/dca_nb_icon.png')
dim_icons <- c('plots/umap_icon.png','plots/paga_umap_icon.png','plots/vae_icon.png','plots/tsne_icon.png', 'plots/dm_icon.png')

scimpute <- readPNG(impute_icons[1])
drimpute <- readPNG(impute_icons[2])

grid.newpage()
grid.raster(scimpute, x = unit(0.12, "npc"), y = unit(0.25, "npc"), width=0.2)
grid.raster(drimpute, x = unit(0.12, "npc"), y = unit(0.75, "npc"), width=0.2)
grid.text('DrImpute')