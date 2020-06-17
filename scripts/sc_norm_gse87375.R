library(Matrix)

## Deconvolution ##
for(d in c('a', 'b')){
  
  # Select Drug
  if(d == 'a'){
    in_file <- '~/sc_eval/data/gse87375/scimpute_res/scimpute_a_out/scimpute_a_outscimpute_count.rds'
    out.file <- '~/sc_eval/data/gse87375/scimpute_res/deconv_res/gse87375_a_scimpute_deconv.rds'
    
  }else if(d == 'b'){
    in_file <- '~/sc_eval/data/gse87375/scimpute_res/scimpute_b_out/scimpute_b_outscimpute_count.rds'
    out.file <- '~/sc_eval/data/gse87375/scimpute_res/deconv_res/gse87375_b_scimpute_deconv.rds'
  }
  
  # Load  Mat 
  count.mat <- round(readRDS(in_file))
  labels <- sc.data$labels[match(colnames(count.mat), table = sc.data$labels$id),]
  
  count.mat <- count.mat[rowSums(count.mat > 0) > 5,]
  
  # Deconvolution
  require(scran)
  sce <- SingleCellExperiment(list(counts = count.mat))
  local.clusters <- quickCluster(sce, 
                                 min.size = 100, 
                                 method = 'hclust', 
                                 d = NA, 
                                 use.ranks = TRUE, 
                                 min.mean = 1)
  sce <- computeSumFactors(sce, clusters = local.clusters)
  
  require(scater)
  sce <- normalize(sce)
  norm.mat <- logcounts(sce)
  saveRDS(list('counts' = count.mat,
               'labels' = labels,
               'norm_mat' = norm.mat), 
          file = out.file)
}

## SCTransform ##
for(d in c('a', 'b')){
  
  # Select Drug
  if(d == 'a'){
    in_file <- '~/sc_eval/data/gse87375/scimpute_res/scimpute_a_out/scimpute_a_outscimpute_count.rds'
    out.file <- '~/sc_eval/data/gse87375/scimpute_res/sctransform_res/gse87375_a_scimpute_sctransform.rds'
    
  }else if(d == 'b'){
    in_file <- '~/sc_eval/data/gse87375/scimpute_res/scimpute_b_out/scimpute_b_outscimpute_count.rds'
    out.file <- '~/sc_eval/data/gse87375/scimpute_res/sctransform_res/gse87375_b_scimpute_sctransform.rds'
  }
  
  # Load  Mat 
  count.mat <- round(readRDS(in_file))
  labels <- sc.data$labels[match(colnames(count.mat), table = sc.data$labels$id),]
  
  count.mat <- count.mat[rowSums(count.mat > 0) > 5,]
  
  # SCtransform
  require(sctransform)
  vst_out <- sctransform::vst(as.matrix(count.mat), 
                              latent_var = c('log_umi'), 
                              return_gene_attr = TRUE, 
                              return_cell_attr = TRUE)
  norm.mat <- vst_out$y
  norm.mat <- norm.mat[!vst_out$model_pars_outliers,]
  
  
  # Save Results
  saveRDS(list('counts' = count.mat,
               'labels' = labels,
               'norm_mat' = norm.mat), 
          file = out.file)
}
