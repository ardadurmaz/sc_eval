library(Matrix)
library(ggplot2)
library(dbscan)
library(xgboost)
library(leiden)
library(parallel)

get_nn <- function(x=NULL, k=15){
  local.d <- as.matrix(dist(x))
  diag(local.d) <- NA
  local.e <- do.call('rbind', sapply(1:nrow(local.d), simplify = FALSE, function(i){
    local.c <- colnames(local.d)[order(local.d[i,], decreasing = FALSE)[1:k]]
    return(data.frame('cell.a' = rep(rownames(local.d)[i], length(local.c)),
                      'cell.b' = local.c))
  }))
  return(local.e)
}

get_clust <- function(edge.ft=NULL){
  require(igraph)
  require(leiden)
  local.net <- graph_from_data_frame(edge.ft, directed = FALSE)
  clust.res <- leiden(local.net, n_iterations = -1, resolution_parameter = 0.1)
  names(clust.res) <- V(local.net)$name
  return(clust.res)
}

set.seed(42)

title_sub <- 'Alectinib'
cl <- makeCluster(7)
clusterExport(cl, varlist = c('get_nn', 'get_clust', 'title_sub'))
clust_res <- parSapply(cl, 1:7, simplify = FALSE, function(i){
  library(Matrix)
  library(ggplot2)
  library(dbscan)
  library(xgboost)
  library(leiden)
  library(parallel)
  
  suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
  suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
  
  local_clust_res <- list()
  for(d_type in c('raw', 'scimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca')){
      for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
        if(d_type == 'scimpute' && n_type == 'dca')
          next
        
        title <- paste0(title_sub, sprintf('_S%d', i))
        tag <- sprintf('~/sc_eval/data/%s/%s_%s_%s_%s.rds', title, title, d_type, n_type, r_type)
        dim_res <- readRDS(tag)
        rownames(dim_res) <- paste0('cell.', 1:nrow(dim_res))
        
        ## Cluster ##
        local.clust <- get_clust(edge.ft = get_nn(dim_res, k = 15))
        local_clust_res[[length(local_clust_res) + 1]] <- local.clust
        names(local_clust_res)[length(local_clust_res)] <- sprintf('%s_%s_%s', d_type, n_type, r_type)
      }
    }
  }
  return(local_clust_res)
})
stopCluster(cl)

saveRDS(clust_res, file = sprintf('~/sc_eval/data/%s_Clusters.rds', title_sub))
quit(save='no')
