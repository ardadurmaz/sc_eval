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
  saveRDS(counts(local_sceval, 'scimpute'), file = sprintf('%s_ScImpute.rds', local_sceval@title))
  counts(local_sceval) <- list('drimpute' = sc_impute(counts(local_sceval, 'raw'), method = 'drimpute'))
  saveRDS(counts(local_sceval, 'drimpute'), file = sprintf('%s_DrImpute.rds', local_sceval@title))
  
  gc()
  setwd('~/sc_eval/data/')
}
cat(sprintf('imputation done..'))
