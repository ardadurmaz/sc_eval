library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd('~/Research/sc_eval/')

title_sub <- 'EGEOD103334'
sling <- readRDS(sprintf('data/%s_Slingshot_Comparison.rds', title_sub))
monoc <- readRDS(sprintf('data/%s_Monocle_Comparison.rds', title_sub))

max_val <- round(max(sapply(monoc, max), sapply(sling, max)), digits = 2)
min_val <- 0

## Slingshot ##
for(dtype in c('raw', 'scimpute')){
  local.mat <- list()
  local.annot <- list()
  for(ntype in c('deconv', 'sctransform', 'dca')){
    if(dtype=='scimpute' && ntype=='dca')
      next
    for(rtype in c('umap', 'tsne', 'dhaka', 'paga_umap')){
      tag <- sprintf('%s_%s_%s_Slingshot.rds', dtype, ntype, rtype)
      comp_mat <- sling[[tag]]
      annot_mat <- data.frame('Data' = rep(dtype, ncol(comp_mat)),
                              'Norm' = rep(ntype, ncol(comp_mat)),
                              'Dim' = rep(rtype, ncol(comp_mat)))
      local.mat[[length(local.mat) + 1]] <- comp_mat
      local.annot[[length(local.annot) + 1]] <- annot_mat
    }
  }
  combined.mat <- do.call('cbind', local.mat)
  combined.annot <- do.call('rbind', local.annot)
  local.mapp.norm <- setNames(c('Deconvolution', 'ScTransform', 'DCA'), c('deconv', 'sctransform', 'dca'))
  local.mapp.dim <- setNames(c('UMAP', 't-SNE', 'Dhaka', 'Paga+UMAP'), c('umap', 'tsne', 'dhaka', 'paga_umap'))
  split.annot <- data.frame('Normalization' = local.mapp.norm[match(combined.annot$Norm, table=names(local.mapp.norm))],
                            'Dimension Reduction' = local.mapp.dim[match(combined.annot$Dim, table=names(local.mapp.dim))])
 
  ## Plot Heatmap ##
  col_fun <- colorRamp2(seq(0,max_val,0.001), viridis(n=length(seq(0,max_val,0.001))))
  combined.mat[combined.mat == 0] <- NA
  ht <- Heatmap(combined.mat,
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(title = 'GW Distance', legend_height = unit(4, 'cm'), title_position = 'leftcenter-rot'),
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_split = split.annot,
                column_title = "%s\n%s",
                column_title_gp = gpar(fontsize=8, fontface='bold'),
                na_col = 'white')
  require(Cairo)
  CairoPDF(file=sprintf('~/Research/sc_eval/%s_Slingshot_%s_Comparison.pdf', title_sub, dtype), width = 12, height = 2)
  draw(ht)
  dev.off()
  
}

## Monocle ##
for(dtype in c('raw', 'scimpute')){
  local.mat <- list()
  local.annot <- list()
  for(ntype in c('deconv', 'sctransform', 'dca')){
    if(dtype=='scimpute' && ntype=='dca')
      next
  
    tag <- sprintf('%s_%s_Monocle', dtype, ntype)
    comp_mat <- monoc[[tag]]
    annot_mat <- data.frame('Data' = rep(dtype, ncol(comp_mat)),
                            'Norm' = rep(ntype, ncol(comp_mat)))
    local.mat[[length(local.mat) + 1]] <- comp_mat
    local.annot[[length(local.annot) + 1]] <- annot_mat
  
  }
  combined.mat <- do.call('cbind', local.mat)
  combined.annot <- do.call('rbind', local.annot)
  local.mapp.norm <- setNames(c('Deconvolution', 'ScTransform', 'DCA'), c('deconv', 'sctransform', 'dca'))
  split.annot <- data.frame('Normalization' = local.mapp.norm[match(combined.annot$Norm, table=names(local.mapp.norm))])
  
  ## Plot Heatmap ##
  col_fun <- colorRamp2(seq(0,max_val,0.001), viridis(n=length(seq(0,max_val,0.001))))
  combined.mat[combined.mat == 0] <- NA
  ht <- Heatmap(combined.mat,
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(title = 'GW Distance', legend_height = unit(4, 'cm'), title_position = 'leftcenter-rot'),
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_split = split.annot,
                column_title = "%s",
                column_title_gp = gpar(fontsize=10, fontface='bold'),
                na_col = 'white')
  require(Cairo)
  CairoPDF(file=sprintf('~/Research/sc_eval/%s_Monocle_%s_Comparison.pdf', title_sub, dtype), width = 8, height = 2)
  draw(ht)
  dev.off()
}
