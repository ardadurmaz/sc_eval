library(Matrix)
library(Rfast)
library(uwot)
library(dbscan)
library(slingshot)

## UMAP ##
all.labels <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1',
                'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3')
common.colors <- setNames(ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all.labels))],
                          all.labels)

for(data_name in c('alevin', 'cellranger')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('sctransform', 'deconv')){
      
      for(drug in c('alectinib', 'lorlatinib', 'crizotinib')){
        message(sprintf('Processing: [%s, %s, %s, umap, %s]', data_name, base_name, norm_name, drug))
        
        if(drug == 'alectinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/alec/umap_res', data_name, base_name, norm_name))
          dim_res <- read.table(sprintf('%s_%s_%s_%s_umap.txt', data_name, drug, base_name, norm_name), header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '')
        }else if(drug == 'lorlatinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/lor/umap_res', data_name, base_name, norm_name))
          dim_res <- read.table(sprintf('%s_%s_%s_%s_umap.txt', data_name, drug, base_name, norm_name), header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '')
        }else{
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/criz/umap_res', data_name, base_name, norm_name))
          dim_res <- read.table(sprintf('%s_%s_%s_%s_umap.txt', data_name, drug, base_name, norm_name), header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '')
        }
        
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
        orig_dims$color <- common.colors[match(dim_res$Label, table = names(common.colors))]
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
                    file = sprintf('slingshot_res/%s_%s_%s_%s_umap_slingshot.txt', data_name, drug, base_name, norm_name),
                    sep = '\t', 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    append = FALSE,
                    quote = FALSE)
        
        # Plot
        system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                       sprintf('slingshot_res/%s_%s_%s_%s_umap_slingshot.txt', data_name, drug, base_name, norm_name),
                       sprintf('slingshot_res/%s_%s_%s_%s_umap_slingshot.pdf', data_name, drug, base_name, norm_name)))
      }
    }
  }
}

## DHAKA ##
for(data_name in c('alevin', 'cellranger')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('sctransform', 'deconv')){
      
      for(drug in c('alectinib', 'lorlatinib', 'crizotinib')){
        message(sprintf('Processing: [%s, %s, %s, %s, %s]', data_name, base_name, norm_name, 'dhaka', drug))
        
        if(drug == 'alectinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/alec/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          dhaka_res <- read.table(sprintf('dhaka_res/%s_%s_%s_%s_dhaka.txt', data_name, drug, base_name, norm_name),
                                  header = FALSE,
                                  sep = '')
          embedd <- data.frame('Dim.1' = dhaka_res[,1],
                               'Dim.2' = dhaka_res[,2],
                               'Dim.3' = dhaka_res[,3],
                               'Color' = common.colors[match(sc_data$labels, table = names(common.colors))],
                               'Label' = sc_data$labels,
                               'Id' = colnames(sc_data$norm_mat))
          
          clust_res <- hdbscan(embedd[,1:3], minPts = 50)
          clust_labs <- setNames(paste0('Cluster-', clust_res$cluster), embedd$Id)
          clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
          
          ## Remove outlier
          embedd <- embedd[match(names(clust_labs), table = embedd$Id),]
          
          ## Slingshot
          lin_res <- getLineages(data = embedd[,1:3], clusterLabels = clust_labs)
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
          
          orig_dims <- as.data.frame(embedd[,1:3], make.names = FALSE)
          orig_dims$color <- common.colors[match(embedd$Label, table = names(common.colors))]
          orig_dims$type <- rep('CB', nrow(orig_dims))
          orig_dims$id <- embedd$Id
          orig_dims$label <- embedd$Label
          colnames(orig_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          curv_dims <- as.data.frame(curv_dims, row.names = NULL)
          curv_dims$color <- rep('#000000', nrow(curv_dims))
          curv_dims$type <- rep('Tree', nrow(curv_dims))
          curv_dims$id <- rep('Root', nrow(curv_dims))
          curv_dims$label <- rep(NA, nrow(curv_dims))
          colnames(curv_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          
          combined_data <- rbind(orig_dims, curv_dims)
          write.table(combined_data, 
                      sprintf('dhaka_res/slingshot_res/%s_%s_%s_%s_dhaka_slingshot.txt', data_name, drug, base_name, norm_name),
                      sep = '\t', 
                      col.names = TRUE, 
                      row.names = FALSE, 
                      append = FALSE,
                      quote = FALSE)
        }else if(drug == 'lorlatinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/lor/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          dhaka_res <- read.table(sprintf('dhaka_res/%s_%s_%s_%s_dhaka.txt', data_name, drug, base_name, norm_name),
                                  header = FALSE,
                                  sep = '')
          embedd <- data.frame('Dim.1' = dhaka_res[,1],
                               'Dim.2' = dhaka_res[,2],
                               'Dim.3' = dhaka_res[,3],
                               'Color' = common.colors[match(sc_data$labels, table = names(common.colors))],
                               'Label' = sc_data$labels,
                               'Id' = colnames(sc_data$norm_mat))
          clust_res <- hdbscan(embedd[,1:3], minPts = 50)
          clust_labs <- setNames(paste0('Cluster-', clust_res$cluster), embedd$Id)
          clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
          
          ## Remove outlier
          embedd <- embedd[match(names(clust_labs), table = embedd$Id),]
          
          ## Slingshot
          lin_res <- getLineages(data = embedd[,1:3], clusterLabels = clust_labs)
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
          
          orig_dims <- as.data.frame(embedd[,1:3], make.names = FALSE)
          orig_dims$color <- common.colors[match(embedd$Label, table = names(common.colors))]
          orig_dims$type <- rep('CB', nrow(orig_dims))
          orig_dims$id <- embedd$Id
          orig_dims$label <- embedd$Label
          colnames(orig_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          curv_dims <- as.data.frame(curv_dims, row.names = NULL)
          curv_dims$color <- rep('#000000', nrow(curv_dims))
          curv_dims$type <- rep('Tree', nrow(curv_dims))
          curv_dims$id <- rep('Root', nrow(curv_dims))
          curv_dims$label <- rep(NA, nrow(curv_dims))
          colnames(curv_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          
          combined_data <- rbind(orig_dims, curv_dims)
          write.table(combined_data, 
                      sprintf('dhaka_res/slingshot_res/%s_%s_%s_%s_dhaka_slingshot.txt', data_name, drug, base_name, norm_name),
                      sep = '\t', 
                      col.names = TRUE, 
                      row.names = FALSE, 
                      append = FALSE,
                      quote = FALSE)
        }else{
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/criz/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          dhaka_res <- read.table(sprintf('dhaka_res/%s_%s_%s_%s_dhaka.txt', data_name, drug, base_name, norm_name),
                                  header = FALSE,
                                  sep = '')
          embedd <- data.frame('Dim.1' = dhaka_res[,1],
                               'Dim.2' = dhaka_res[,2],
                               'Dim.3' = dhaka_res[,3],
                               'Color' = common.colors[match(sc_data$labels, table = names(common.colors))],
                               'Label' = sc_data$labels,
                               'Id' = colnames(sc_data$norm_mat))
          clust_res <- hdbscan(embedd[,1:3], minPts = 50)
          clust_labs <- setNames(paste0('Cluster-', clust_res$cluster), embedd$Id)
          clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
          
          ## Remove outlier
          embedd <- embedd[match(names(clust_labs), table = embedd$Id),]
          if(nrow(embedd) < 10){
            next
          }
          ## Slingshot
          lin_res <- getLineages(data = embedd[,1:3], clusterLabels = clust_labs)
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
          
          orig_dims <- as.data.frame(embedd[,1:3], make.names = FALSE)
          orig_dims$color <- common.colors[match(embedd$Label, table = names(common.colors))]
          orig_dims$type <- rep('CB', nrow(orig_dims))
          orig_dims$id <- embedd$Id
          orig_dims$label <- embedd$Label
          colnames(orig_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          curv_dims <- as.data.frame(curv_dims, row.names = NULL)
          curv_dims$color <- rep('#000000', nrow(curv_dims))
          curv_dims$type <- rep('Tree', nrow(curv_dims))
          curv_dims$id <- rep('Root', nrow(curv_dims))
          curv_dims$label <- rep(NA, nrow(curv_dims))
          colnames(curv_dims) <- c('Dim.1', 'Dim.2', 'Dim.3', 'Color', 'Type', 'Id', 'Label')
          
          
          combined_data <- rbind(orig_dims, curv_dims)
          write.table(combined_data, 
                      sprintf('dhaka_res/slingshot_res/%s_%s_%s_%s_dhaka_slingshot.txt', data_name, drug, base_name, norm_name),
                      sep = '\t', 
                      col.names = TRUE, 
                      row.names = FALSE, 
                      append = FALSE,
                      quote = FALSE)
        }
        # Plot
        system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                       sprintf('dhaka_res/slingshot_res/%s_%s_%s_%s_dhaka_slingshot.txt', data_name, drug, base_name, norm_name),
                       sprintf('dhaka_res/slingshot_res/%s_%s_%s_%s_dhaka_slingshot.pdf', data_name, drug, base_name, norm_name)))
      }
    }
  }
}
