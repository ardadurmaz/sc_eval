library(Matrix)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))

setwd('~/sc_eval/data/')
set.seed(42)
CombinedData <- readRDS('PancreaticMaturation_B.rds')

## Run ##
for(i in 1:length(CombinedData)){
  
  local_sceval <- CombinedData[[i]]
  log_file <- paste0(local_sceval@title, '.log')
  wd <- sprintf('~/sc_eval/data/%s', local_sceval@title)
  setwd(wd)
  
  ## Impute ##
  counts(local_sceval) <- list('scimpute' = readRDS(sprintf('%s_ScImpute.rds', local_sceval@title)))
  
  ## Normalize ##
  norm(local_sceval) <- list('raw_deconv' = readRDS(sprintf('%s_RawDeconv.rds', local_sceval@title)))
  norm(local_sceval) <- list('raw_sctransform' = readRDS(sprintf('%s_RawScTransform.rds', local_sceval@title)))
  norm(local_sceval) <- list('raw_dca_nb' = readRDS(sprintf('%s_RawDCA_NB.rds', local_sceval@title)))
  norm(local_sceval) <- list('raw_dca_zinb' = readRDS(sprintf('%s_RawDCA_ZINB.rds', local_sceval@title)))
  norm(local_sceval) <- list('scimpute_deconv' = readRDS(sprintf('%s_ScImputeDeconv.rds', local_sceval@title)))
  norm(local_sceval) <- list('scimpute_sctransform' = readRDS(sprintf('%s_ScImputeScTransform.rds', local_sceval@title)))
  norm(local_sceval) <- list('drimpute_deconv' = readRDS(sprintf('%s_DrImputeDeconv.rds', local_sceval@title)))
  norm(local_sceval) <- list('drimpute_sctransform' = readRDS(sprintf('%s_DrImputeScTransform.rds', local_sceval@title)))
  
  ## Dim Reduce ##
  cat(sprintf('SC_EVAL [%s]: Running dimension reduction...\n', date()), file = log_file, append = TRUE)
  for(d_type in c('raw', 'scimpute', 'drimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca_nb', 'dca_zinb')){
      for(r_type in c('umap', 'tsne', 'vae', 'dm', 'paga_umap')){
        if((d_type == 'scimpute' || d_type == 'drimpute') && (n_type == 'dca_nb' || n_type == 'dca_zinb'))
          next
        cat(sprintf('[%s: %s %s %s]\n', local_sceval@title, d_type, n_type, r_type))
        tag <- sprintf('%s_%s_%s', d_type, n_type, r_type)
        m_tag <- sprintf('%s_%s', d_type, n_type)
        temp_dim <- sc_dimR(norm(local_sceval, m_tag), method = r_type)
        if(nrow(temp_dim) != ncol(norm(local_sceval, m_tag)))
          warning('Reduced dimensions do not match input matrix')
        temp_list <- list(temp_dim)
        names(temp_list) <- tag
        redDim(local_sceval) <- temp_list
        saveRDS(redDim(local_sceval, tag), file = sprintf('%s_%s.rds', local_sceval@title, tag))
      }
    }
  }
}
cat(sprintf('dimension reduction done..'))
