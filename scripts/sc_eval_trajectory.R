library(Matrix)
library(parallel)


setwd('~/sc_eval/data/')
set.seed(42)
CombinedData <- readRDS('GSE107863_Data.rds')

# Select Lorlatinib #
#CombinedData <- CombinedData[8:14]
cl <- makeCluster(7)

run_res <- parSapply(cl, CombinedData, simplify = FALSE, function(local_sceval){
  suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
  suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
  
  log_file <- paste0(local_sceval@title, '.log')
  wd <- sprintf('~/sc_eval/data/%s', local_sceval@title)
  setwd(wd)
  
  ## Impute ##
  counts(local_sceval) <- list('scimpute' = readRDS(sprintf('%s_ScImpute.rds', local_sceval@title)))
  
  ## Normalize ##
  norm(local_sceval) <- list('raw_deconv' = readRDS(sprintf('%s_RawDeconv.rds', local_sceval@title)))
  norm(local_sceval) <- list('raw_sctransform' = readRDS(sprintf('%s_RawScTransform.rds', local_sceval@title)))
  norm(local_sceval) <- list('raw_dca' = readRDS(sprintf('%s_RawDCA.rds', local_sceval@title)))
  norm(local_sceval) <- list('scimpute_deconv' = readRDS(sprintf('%s_ScImputeDeconv.rds', local_sceval@title)))
  norm(local_sceval) <- list('scimpute_sctransform' = readRDS(sprintf('%s_ScImputeScTransform.rds', local_sceval@title)))
  
  ## Dim Reduce ##
  for(d_type in c('raw', 'scimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca')){
      for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
        if(d_type == 'scimpute' && n_type == 'dca')
          next
        tag <- sprintf('%s_%s_%s', d_type, n_type, r_type)
        redDim(local_sceval) <- setNames(list(readRDS(sprintf('%s_%s.rds', local_sceval@title, tag))), tag)
      }
    }
  }
  
  ## Trajectory Inference ##
  cat(sprintf('SC_EVAL [%s]: Running trajectory inference...\n', date()), file = log_file, append = TRUE)
  local_sceval <- sc_traj(local_sceval, method = 'slingshot')
  #local_sceval <- sc_traj(local_sceval, method = 'monocle')
  #local_sceval <- sc_traj(local_sceval, method = 'OT')
  gc()
})
stopCluster(cl)
cat(sprintf('trajectory inference done..'))
