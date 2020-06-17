library(Matrix)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))

setwd('~/sc_eval/data/')
set.seed(42)
CombinedData <- readRDS('MCF7_Counts.rds')

## Run ##
for(i in 6:length(CombinedData)){
  
  local_sceval <- CombinedData[[i]]
  
  ## Create Sample Data ##
  #sample_idx <- sample(1:ncol(counts(local_sceval, 'raw')), size = 0.4*ncol(counts(local_sceval, 'raw')), replace = FALSE)

  #count_mat <- round(as.matrix(counts(local_sceval, 'raw')[,sample_idx]))
  #count_mat <- count_mat[rowSums(count_mat > 0) >= 50,]
  #annot_data <- data.frame('Id' = colnames(count_mat),
  #                         'Label' = local_sceval@annot$Label[sample_idx],
  #                         'Days' = local_sceval@annot$Days[sample_idx],
  #                         'Color' = local_sceval@annot$Color[sample_idx])
  #sample_data <- new('SCEVAL', 
  #                   title = sprintf('SampleData_%d', i), 
  #                   counts = list('raw' = Matrix(count_mat)), 
  #                   annot = annot_data)
  #local_sceval <- sample_data
  ########################
  
  log_file <- paste0(local_sceval@title, '.log')
  wd <- sprintf('%s', local_sceval@title)
  
  ## Setup Dir ##
  if(!dir.exists(wd))
    tryCatch({
      dir.create(wd)
    }, error = function(err){
      message('Failed to create directory: ', wd)
      next
    })
  setwd(wd)
  
  ## Impute ##
  cat(sprintf('SC_EVAL [%s]: Running imputation...\n', date()), file = log_file, append = TRUE)
  
  counts(local_sceval) <- list('scimpute' = sc_impute(counts(local_sceval, 'raw'), method = 'scimpute'))
  
  ## Normalize ##
  cat(sprintf('SC_EVAL [%s]: Running normalization...\n', date()), file = log_file, append = TRUE)
  
  norm(local_sceval) <- list('raw_deconv' = sc_norm(counts(local_sceval, 'raw'), method = 'deconv'))
  norm(local_sceval) <- list('raw_sctransform' = sc_norm(counts(local_sceval, 'raw'), method = 'sctransform'))
  norm(local_sceval) <- list('raw_dca' = sc_norm(counts(local_sceval, 'raw'), method = 'dca'))
  norm(local_sceval) <- list('scimpute_deconv' = sc_norm(counts(local_sceval, 'scimpute'), method = 'deconv'))
  norm(local_sceval) <- list('scimpute_sctransform' = sc_norm(counts(local_sceval, 'scimpute'), method = 'sctransform'))
  
  
  ## Dimension Reduction ##
  cat(sprintf('SC_EVAL [%s]: Running dimension reduction...\n', date()), file = log_file, append = TRUE)
  
  redDim(local_sceval) <- list('raw_deconv_umap' = sc_dimR(norm(local_sceval, 'raw_deconv'), method = 'umap'))
  redDim(local_sceval) <- list('raw_deconv_tsne' = sc_dimR(norm(local_sceval, 'raw_deconv'), method = 'tsne'))
  redDim(local_sceval) <- list('raw_deconv_dhaka' = sc_dimR(norm(local_sceval, 'raw_deconv'), method = 'dhaka'))
  redDim(local_sceval) <- list('raw_deconv_dm' = sc_dimR(norm(local_sceval, 'raw_deconv'), method = 'dm'))
  
  redDim(local_sceval) <- list('raw_sctransform_umap' = sc_dimR(norm(local_sceval, 'raw_sctransform'), method = 'umap'))
  redDim(local_sceval) <- list('raw_sctransform_tsne' = sc_dimR(norm(local_sceval, 'raw_sctransform'), method = 'tsne'))
  redDim(local_sceval) <- list('raw_sctransform_dhaka' = sc_dimR(norm(local_sceval, 'raw_sctransform'), method = 'dhaka'))
  redDim(local_sceval) <- list('raw_sctransform_dm' = sc_dimR(norm(local_sceval, 'raw_sctransform'), method = 'dm'))

  redDim(local_sceval) <- list('raw_dca_umap' = sc_dimR(norm(local_sceval, 'raw_dca'), method = 'umap'))
  redDim(local_sceval) <- list('raw_dca_tsne' = sc_dimR(norm(local_sceval, 'raw_dca'), method = 'tsne'))
  redDim(local_sceval) <- list('raw_dca_dhaka' = sc_dimR(norm(local_sceval, 'raw_dca'), method = 'dhaka'))
  redDim(local_sceval) <- list('raw_dca_dm' = sc_dimR(norm(local_sceval, 'raw_dca'), method = 'dm'))
  
  redDim(local_sceval) <- list('scimpute_deconv_umap' = sc_dimR(norm(local_sceval, 'scimpute_deconv'), method = 'umap'))
  redDim(local_sceval) <- list('scimpute_deconv_tsne' = sc_dimR(norm(local_sceval, 'scimpute_deconv'), method = 'tsne'))
  redDim(local_sceval) <- list('scimpute_deconv_dhaka' = sc_dimR(norm(local_sceval, 'scimpute_deconv'), method = 'dhaka'))
  redDim(local_sceval) <- list('scimpute_deconv_dm' = sc_dimR(norm(local_sceval, 'scimpute_deconv'), method = 'dm'))
  
  redDim(local_sceval) <- list('scimpute_sctransform_umap' = sc_dimR(norm(local_sceval, 'scimpute_sctransform'), method = 'umap'))
  redDim(local_sceval) <- list('scimpute_sctransform_tsne' = sc_dimR(norm(local_sceval, 'scimpute_sctransform'), method = 'tsne'))
  redDim(local_sceval) <- list('scimpute_sctransform_dhaka' = sc_dimR(norm(local_sceval, 'scimpute_sctransform'), method = 'dhaka'))
  redDim(local_sceval) <- list('scimpute_sctransform_dm' = sc_dimR(norm(local_sceval, 'scimpute_sctransform'), method = 'dm'))
  
  ## Trajectory Inference ##
  cat(sprintf('SC_EVAL [%s]: Running trajectory inference...\n', date()), file = log_file, append = TRUE)
  local_sceval <- sc_traj(local_sceval, method = 'slingshot')
  local_sceval <- sc_traj(local_sceval, method = 'monocle')
  
  cat(sprintf('SC_EVAL [%s]: Saving results...\n', date()), file = log_file, append = TRUE)
  saveRDS(local_sceval, file = sprintf('%s_ProcessedData.rds', local_sceval@title))

  setwd('~/sc_eval/data')
  gc()
}
