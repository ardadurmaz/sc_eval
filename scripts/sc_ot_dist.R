library(Matrix)
library(ggplot2)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

title_sub <- 'GSE107863'

## Within workflow ##
all_res <- list()
for(d_type in c('raw', 'scimpute')){
  for(n_type in c('deconv', 'sctransform', 'dca')){
    for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      
      traject_res <- sapply(1:7, simplify = FALSE, function(i){
        title <- paste0(title_sub, sprintf('_S%d', i))
        tag <- sprintf('~/sc_eval/data/%s/%s_%s_%s_%s_Slingshot.rds', title, title, d_type, n_type, r_type)
        traj_res <- readRDS(tag)[[1]]
        traj_res <- subset(traj_res, is.na(traj_res$Label))
        return(traj_res)
      })
      comp_res <- sapply(1:(length(traject_res)-1), simplify = FALSE, function(i){
        local_comp_res <- sapply((i+1):length(traject_res), simplify = TRUE, function(j){
          res.1 <- traject_res[[i]]
          res.2 <- traject_res[[j]]
          
          res.1 <- res.1[sample(1:nrow(res.1), size = 15, replace = FALSE),]
          res.2 <- res.2[sample(1:nrow(res.2), size = 15, replace = FALSE),]
          dist.1 <- as.matrix(dist(as.matrix(res.1[,1:2])))
          dist.2 <- as.matrix(dist(as.matrix(res.2[,1:2])))
          dist.1 <- dist.1/max(dist.1)
          dist.2 <- dist.2/max(dist.2)
          
          write.table(dist.1, file = '~/sc_eval/data/DistMat_A.tsv', col.names = FALSE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)
          write.table(dist.2, file = '~/sc_eval/data/DistMat_B.tsv', col.names = FALSE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)
          
          gw_dist <- as.numeric(system("python3 /home/arda/sc_eval/scripts/calculate_gw.py", intern=TRUE))
          return(gw_dist)
        })
        return(local_comp_res)
      })
      dist_res <- matrix(0, ncol=7, nrow=7)
      for(i in 1:6){
        dist_res[i,(i+1):7] <- comp_res[[i]]
      }
      all_res[[length(all_res)+1]] <- dist_res
      names(all_res)[length(all_res)] <- sprintf('%s_%s_%s_Slingshot.rds', d_type, n_type, r_type)
    }
  }
}
saveRDS(all_res, file = sprintf('~/sc_eval/data/%s_Slingshot_Comparison.rds', title_sub))



all_res <- list()
for(d_type in c('raw', 'scimpute')){
  for(n_type in c('deconv', 'sctransform', 'dca')){
    if(d_type == 'scimpute' && n_type == 'dca')
      next
      
    traject_res <- sapply(1:7, simplify = FALSE, function(i){
      title <- paste0(title_sub, sprintf('_S%d', i))
      tag <- sprintf('~/sc_eval/data/%s/%s_%s_%s_Monocle2.rds', title, title, d_type, n_type)
      traj_res <- readRDS(tag)[[1]]
      traj_res <- subset(traj_res, is.na(traj_res$Label))
      return(traj_res)
    })
    comp_res <- sapply(1:(length(traject_res)-1), simplify = FALSE, function(i){
      local_comp_res <- sapply((i+1):length(traject_res), simplify = TRUE, function(j){
        res.1 <- traject_res[[i]]
        res.2 <- traject_res[[j]]
          
        res.1 <- res.1[sample(1:nrow(res.1), size = 15, replace = FALSE),]
        res.2 <- res.2[sample(1:nrow(res.2), size = 15, replace = FALSE),]
        dist.1 <- as.matrix(dist(as.matrix(res.1[,1:2])))
        dist.2 <- as.matrix(dist(as.matrix(res.2[,1:2])))
        dist.1 <- dist.1/max(dist.1)
        dist.2 <- dist.2/max(dist.2)
          
        write.table(dist.1, file = '~/sc_eval/data/DistMat_A.tsv', col.names = FALSE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)
        write.table(dist.2, file = '~/sc_eval/data/DistMat_B.tsv', col.names = FALSE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE)
          
        gw_dist <- as.numeric(system("python3 /home/arda/sc_eval/scripts/calculate_gw.py", intern=TRUE))
        return(gw_dist)
      })
      return(local_comp_res)
    })
    dist_res <- matrix(0, ncol=7, nrow=7)
    for(i in 1:6){
      dist_res[i,(i+1):7] <- comp_res[[i]]
    }
    all_res[[length(all_res)+1]] <- dist_res
    names(all_res)[length(all_res)] <- sprintf('%s_%s_Monocle', d_type, n_type)
  }
}
saveRDS(all_res, file = sprintf('~/sc_eval/data/%s_Monocle_Comparison.rds', title_sub))
