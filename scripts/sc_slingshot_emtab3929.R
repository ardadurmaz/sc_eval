library(Matrix)
library(Rfast)
library(uwot)
library(dbscan)
library(slingshot)

all_labels <- c('embryonic day 3', 'embryonic day 4', 'embryonic day 5', 'embryonic day 6', 'embryonic day 7')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 5, 'Set2')[as.numeric(as.factor(all_labels))], all_labels)

setwd('~/sc_eval/data/emtab3929_data/')

for(base_name in c('raw', 'scimpute')){
  for(norm_name in c('deconv', 'sctransform')){
    for(dim_name in c('umap', 'dhaka')){
      setwd(sprintf('~/sc_eval/data/emtab3929_data/%s_res/%s_res', base_name, norm_name))
      sc_data <- readRDS(sprintf('emtab3929_%s_%s.rds', base_name, norm_name))
      if(dim_name == 'umap'){
        dim_res <- read.table(sprintf('%s_res/emtab3929_%s_%s_%s.txt', dim_name, base_name, norm_name, dim_name), sep = '\t', comment.char = '', header = TRUE)  
        
        ## Cluster
        clust_res <- hdbscan(dim_res[,1:3], minPts = 50)
        clust_labs <- setNames(paste0('Cluster-', clust_res$cluster), dim_res$Id)
        clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
        
        ## Remove outlier
        dim_res <- dim_res[match(names(clust_labs), table = dim_res$Id),]
        
        ## Slingshot
        lin_res <- getLineages(data = dim_res[,1:3], clusterLabels = clust_labs)
        curv_res <- getCurves(lin_res)
        
        ## 3D Plot
        curv_dims <- do.call('rbind', sapply(curv_res@curves, simplify = FALSE, function(x){
          local_pos <- x$s
          local_w <- x$w[match(rownames(local_pos), table = names(x$w))]
          local_ord <- x$ord
          
          local_pos <- local_pos[local_ord,]
          local_w <- local_w[local_ord]
          
          idx <- which(local_w > 0)
          local_pos <- local_pos[idx,]
          return(local_pos)
        }))
        
        orig_dims <- as.data.frame(dim_res[,1:3], make.names = FALSE)
        orig_dims$color <- common_colors[match(dim_res$Label, table = names(common_colors))]
        orig_dims$type <- rep('CB', nrow(orig_dims))
        orig_dims$id <- dim_res$Id
        orig_dims$label <- dim_res$Label
        colnames(orig_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
        
        curv_dims <- as.data.frame(curv_dims, row.names = NULL)
        curv_dims$color <- rep('#000000', nrow(curv_dims))
        curv_dims$type <- rep('Tree', nrow(curv_dims))
        curv_dims$id <- rep('Root', nrow(curv_dims))
        curv_dims$label <- rep(NA, nrow(curv_dims))
        colnames(curv_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
        
        
        combined_data <- rbind(orig_dims, curv_dims)
        write.table(combined_data, 
                    file = sprintf('umap_res/slingshot_res/emtab3929_%s_%s_umap_slingshot.txt', base_name, norm_name),
                    sep = '\t', 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    append = FALSE,
                    quote = FALSE)
        
        # Plot
        system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                       sprintf('umap_res/slingshot_res/emtab3929_%s_%s_umap_slingshot.txt', base_name, norm_name),
                       sprintf('umap_res/slingshot_res/emtab3929_%s_%s_umap_slingshot.pdf', base_name, norm_name)))
        
      }else if(dim_name == 'dhaka'){
        dim_res <- read.table(sprintf('%s_res/emtab3929_%s_%s_%s.txt.txt', dim_name, base_name, norm_name, dim_name), sep = '', comment.char = '', header = FALSE)  
        dim_res$Id <- sc_data$labels$Id
        dim_res$Label <- sc_data$labels$Label
        
        ## Cluster
        clust_res <- hdbscan(dim_res[,1:3], minPts = 50)
        clust_labs <- setNames(paste0('Cluster-', clust_res$cluster), dim_res$Id)
        clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
        
        ## Remove outlier
        dim_res <- dim_res[match(names(clust_labs), table = dim_res$Id),]
        
        ## Slingshot
        lin_res <- getLineages(data = dim_res[,1:3], clusterLabels = clust_labs)
        curv_res <- getCurves(lin_res)
        
        ## 3D Plot
        curv_dims <- do.call('rbind', sapply(curv_res@curves, simplify = FALSE, function(x){
          local_pos <- x$s
          local_w <- x$w[match(rownames(local_pos), table = names(x$w))]
          local_ord <- x$ord
          
          local_pos <- local_pos[local_ord,]
          local_w <- local_w[local_ord]
          
          idx <- which(local_w > 0)
          local_pos <- local_pos[idx,]
          return(local_pos)
        }))
        
        orig_dims <- as.data.frame(dim_res[,1:3], make.names = FALSE)
        orig_dims$color <- common_colors[match(dim_res$Label, table = names(common_colors))]
        orig_dims$type <- rep('CB', nrow(orig_dims))
        orig_dims$id <- dim_res$Id
        orig_dims$label <- dim_res$Label
        colnames(orig_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
        
        curv_dims <- as.data.frame(curv_dims, row.names = NULL)
        curv_dims$color <- rep('#000000', nrow(curv_dims))
        curv_dims$type <- rep('Tree', nrow(curv_dims))
        curv_dims$id <- rep('Root', nrow(curv_dims))
        curv_dims$label <- rep(NA, nrow(curv_dims))
        colnames(curv_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
        
        
        combined_data <- rbind(orig_dims, curv_dims)
        write.table(combined_data, 
                    file = sprintf('dhaka_res/slingshot_res/emtab3929_%s_%s_dhaka_slingshot.txt', base_name, norm_name),
                    sep = '\t', 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    append = FALSE,
                    quote = FALSE)
        
        # Plot
        system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                       sprintf('dhaka_res/slingshot_res/emtab3929_%s_%s_dhaka_slingshot.txt', base_name, norm_name),
                       sprintf('dhaka_res/slingshot_res/emtab3929_%s_%s_dhaka_slingshot.pdf', base_name, norm_name)))
        
      }
    }
  }
}
