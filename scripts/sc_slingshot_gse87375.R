library(Matrix)
library(Rfast)
library(uwot)
library(dbscan)
library(slingshot)

all_labels <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 7, 'Dark2'), all_labels)

## UMAP ##
for(data_name in c('a', 'b')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('deconv', 'sctransform')){
      cat(sprintf('\n*** Processing: %s, %s, %s ***\n', data_name, base_name, norm_name))
      setwd(sprintf('~/sc_eval/data/gse87375/%s_res/%s_res/umap_res', base_name, norm_name))
      dim_res <- read.table(file = sprintf('gse87375_%s_%s_%s_umap.txt', data_name, base_name, norm_name),
                            sep = '\t',
                            header = TRUE,
                            comment = '',
                            stringsAsFactors = FALSE)
      ## Cluster
      clust_res <- hdbscan(dim_res[,1:3], minPts = 10)
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
                  file = sprintf('slingshot_res/gse87375_%s_%s_%s_umap_slingshot.txt', data_name, base_name, norm_name),
                  sep = '\t', 
                  col.names = TRUE, 
                  row.names = FALSE, 
                  append = FALSE,
                  quote = FALSE)
      
      # Plot
      system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                     sprintf('slingshot_res/gse87375_%s_%s_%s_umap_slingshot.txt', data_name, base_name, norm_name),
                     sprintf('slingshot_res/gse87375_%s_%s_%s_umap_slingshot.pdf', data_name, base_name, norm_name)))
      
    }
  }
}

## Dhaka ##
for(data_name in c('a', 'b')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('deconv', 'sctransform')){
      cat(sprintf('\n*** Processing: %s, %s, %s ***\n', data_name, base_name, norm_name))
      setwd(sprintf('~/sc_eval/data/gse87375/%s_res/%s_res/dhaka_res', base_name, norm_name))
      sc_data <- readRDS(sprintf('../gse87375_%s_%s_%s.rds', data_name, base_name, norm_name))
      dim_res <- read.table(file = sprintf('gse87375_%s_%s_%s_dhaka.txt.txt', data_name, base_name, norm_name),
                            sep = '',
                            header = FALSE,
                            comment = '',
                            stringsAsFactors = FALSE)
      colnames(dim_res) <- c('Dim.1', 'Dim.2', 'Dim.3')
      dim_res$Id <- colnames(sc_data$norm_mat)  
      dim_res$Label <- sc_data$labels$time
      
      ## Cluster
      clust_res <- hdbscan(dim_res[,1:3], minPts = 5)
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
                  file = sprintf('slingshot_res/gse87375_%s_%s_%s_dhaka_slingshot.txt', data_name, base_name, norm_name),
                  sep = '\t', 
                  col.names = TRUE, 
                  row.names = FALSE, 
                  append = FALSE,
                  quote = FALSE)
      
      # Plot
      system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                     sprintf('slingshot_res/gse87375_%s_%s_%s_dhaka_slingshot.txt', data_name, base_name, norm_name),
                     sprintf('slingshot_res/gse87375_%s_%s_%s_dhaka_slingshot.pdf', data_name, base_name, norm_name)))
      
    }
  }
}
