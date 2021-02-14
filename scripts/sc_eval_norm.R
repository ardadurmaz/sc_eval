library(Matrix)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))

setwd('~/sc_eval/data/')
set.seed(42)
CombinedData <- readRDS('E2_Treatment_T47D.rds')

## Run ##
for(i in 1:length(CombinedData)){
  
  local_sceval <- CombinedData[[i]]
  log_file <- paste0(local_sceval@title, '.log')
  wd <- sprintf('~/sc_eval/data/%s', local_sceval@title)
  setwd(wd)

  ## Impute ##
  counts(local_sceval) <- list('scimpute' = readRDS(sprintf('%s_ScImpute.rds', local_sceval@title)))
  counts(local_sceval) <- list('drimpute' = readRDS(sprintf('%s_DrImpute.rds', local_sceval@title)))
  
  
  
  ## Normalize ##
  cat(sprintf('SC_EVAL [%s]: Running normalization...\n', date()), file = log_file, append = TRUE)
  norm(local_sceval) <- list('raw_deconv' = sc_norm(counts(local_sceval, 'raw'), method = 'deconv'))
  saveRDS(norm(local_sceval, 'raw_deconv'), file = sprintf('%s_RawDeconv.rds', local_sceval@title))
  
  norm(local_sceval) <- list('raw_sctransform' = sc_norm(counts(local_sceval, 'raw'), method = 'sctransform'))
  saveRDS(norm(local_sceval, 'raw_sctransform'), file = sprintf('%s_RawScTransform.rds', local_sceval@title))
  
  norm(local_sceval) <- list('raw_dca' = sc_norm(counts(local_sceval, 'raw'), method = 'dca_nb'))
  saveRDS(norm(local_sceval, 'raw_dca'), file = sprintf('%s_RawDCA_NB.rds', local_sceval@title))

  norm(local_sceval) <- list('raw_dca_zinb' = sc_norm(counts(local_sceval, 'raw'), method = 'dca_zinb'))
  saveRDS(norm(local_sceval, 'raw_dca_zinb'), file = sprintf('%s_RawDCA_ZINB.rds', local_sceval@title))
  
  norm(local_sceval) <- list('scimpute_deconv' = sc_norm(counts(local_sceval, 'scimpute'), method = 'deconv'))
  saveRDS(norm(local_sceval, 'scimpute_deconv'), file = sprintf('%s_ScImputeDeconv.rds', local_sceval@title))
  
  norm(local_sceval) <- list('scimpute_sctransform' = sc_norm(counts(local_sceval, 'scimpute'), method = 'sctransform'))
  saveRDS(norm(local_sceval, 'scimpute_sctransform'), file = sprintf('%s_ScImputeScTransform.rds', local_sceval@title))
  
  norm(local_sceval) <- list('drimpute_deconv' = sc_norm(counts(local_sceval, 'drimpute'), method = 'deconv'))
  saveRDS(norm(local_sceval, 'drimpute_deconv'), file = sprintf('%s_DrImputeDeconv.rds', local_sceval@title))
  
  norm(local_sceval) <- list('drimpute_sctransform' = sc_norm(counts(local_sceval, 'drimpute'), method = 'sctransform'))
  saveRDS(norm(local_sceval, 'drimpute_sctransform'), file = sprintf('%s_DrImputeScTransform.rds', local_sceval@title))
  
}
cat(sprintf('normalization done..'))
