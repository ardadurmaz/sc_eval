library(Matrix)
library(xgboost)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

title_sub <- 'Lorlatinib'

## Cluster Results ##
clust_res <- readRDS(sprintf('~/sc_eval/data/%s_Clusters.rds', title_sub))

## Within workflow ##
all_res <- list()
for(d_type in c('raw', 'scimpute')){
  for(n_type in c('deconv', 'sctransform', 'dca')){
    for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      cat(sprintf('Processing: %s %s %s\n', d_type, n_type, r_type))
      imp_res <- sapply(1:7, simplify = FALSE, function(i){
        title <- paste0(title_sub, sprintf('_S%d', i))
        local_clust_res <- clust_res[[i]][[which(names(clust_res[[i]]) == sprintf('%s_%s_%s', d_type, n_type, r_type))]]
        
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
        
        xgb.train.data <- xgb.DMatrix(data = as.matrix(t(expr_res)),
                                      label = local_clust_res - 1)
        xgb.params <- list('subsample' = 0.8,
                           'colsample_bytree' = 1,
                           'objective' = 'multi:softmax',
                           'num_class' = length(unique(local_clust_res)),
                           'eval_metric' = 'mlogloss',
                           'max_depth' = 6,
                           'eta' = 0.1,
                           'min_child_weight' = 10)
        xgb.fit <- xgb.train(params = xgb.params,
                             data = xgb.train.data,
                             nrounds = 50,
                             verbose = 1)
        local_imp_res <- xgb.importance(model=xgb.fit)
        return(setNames(local_imp_res$Gain, local_imp_res$Feature))
      })
      all_res[[length(all_res)+1]] <- imp_res
      names(all_res)[length(all_res)] <- sprintf("%s_%s_%s", d_type, n_type, r_type)
    }
  }
}
saveRDS(all_res, file = sprintf('~/sc_eval/data/%s_Importance.rds', title_sub))
quit(save='no')
