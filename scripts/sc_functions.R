sc_norm <- function(counts = NULL, method = 'deconv'){
  if(is.null(counts))
    stop('Provide count matrix')
  
  cat(sprintf("%s\t%.4f\n", method, as.vector(system.time({
    if(method == 'deconv'){
      suppressMessages(require(scran))
      suppressMessages(require(scater))
      message('\n*** Deconvolution ***\n')
      sce <- SingleCellExperiment(list(counts = counts))
      local_clusters <- quickCluster(sce, 
                                     min.size = round(0.1*ncol(counts)), 
                                     method = 'igraph', 
                                     d = NA, 
                                     use.ranks = TRUE, 
                                     min.mean = 1)
      sce <- computeSumFactors(sce, clusters = local_clusters)
      sce <- normalize(sce)
      norm_mat <- as.matrix(logcounts(sce))
    }else if(method == 'sctransform'){
      suppressMessages(require(sctransform))
      message('\n*** ScTransform ***\n')
      vst_out <- sctransform::vst(as.matrix(counts), 
                                  latent_var = c('log_umi'), 
                                  return_gene_attr = FALSE, 
                                  return_cell_attr = FALSE, 
                                  show_progress = FALSE)
      norm_mat <- as.matrix(vst_out$y)
    }else if(method == 'dca'){
      message('\n*** DCA ***\n')
      if(dir.exists('DCA_OUT') && file.exists('DCA_OUT/mean.tsv')){
        dca_out <- read.table('DCA_OUT/mean.tsv', header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)
      }else{
        write.table(as.matrix(counts), file='CountMat.tsv', sep = '\t', col.names = TRUE, row.names = TRUE, append = FALSE, quote = FALSE)
        system(command = sprintf('dca -d 0.1 CountMat.tsv DCA_OUT'))
        dca_out <- read.table('DCA_OUT/mean.tsv', header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)
      }
      row_names <- dca_out[,1]
      dca_out <- Matrix(data.matrix(dca_out[,-1]))
      rownames(dca_out) <- row_names
      
      # Normalize
      lib_size <- colSums(dca_out, na.rm = TRUE)
      norm_mat <- log(1 + (dca_out %*% Matrix::Diagonal(x = (lib_size/median(lib_size))^-1)))
      norm_mat <- as.matrix(norm_mat)
    }
  }))[1]), file = "RunTime.txt", append = TRUE)
  return(norm_mat)
}

sc_dimR <- function(expr = NULL, method = 'umap'){
  if(is.null(expr))
    stop('Provide normalized expressions')
  
  ## Check for 0 Variation ##
  idx <- apply(expr, 1, function(x){var(x, na.rm = TRUE)}) > 0
  expr <- expr[idx,]
  
  cat(sprintf("%s\t%.4f\n", method, as.vector(system.time({
    if(method=='umap'){
      suppressMessages(require(uwot))
      message('\n*** UMAP ***\n')
      dim_res <- umap(as.matrix(t(expr)),
                      n_neighbors = 10,
                      n_components = 2,
                      n_epochs = 800,
                      scale = TRUE,
                      init = 'spectral',
                      min_dist = 0.01,
                      spread = 1,
                      pca = NULL,
                      n_threads = 20)
    }else if(method=='tsne'){
      suppressMessages(require(Rtsne))
      message('\n*** t-SNE ***\n')
      tsne_res <- Rtsne(scale(as.matrix(t(expr))),
                        dims = 2,
                        perplexity = 10,
                        theta = 0.1,
                        partial_pca = FALSE,
                        pca_center = FALSE,
                        pca_scale = FALSE,
                        num_threads = 20)
      dim_res <- tsne_res$Y
    }else if(method=='paga_umap'){
      message('\n*** PAGA ***\n')
      out_file <- sprintf('%s/NormExpr.txt', getwd())
      write.table(t(as.matrix(expr)),
                  file = out_file, 
                  col.names = FALSE, 
                  row.names = FALSE,
                  sep = '\t',
                  append = FALSE,
                  quote = FALSE)
      system(command=sprintf('python3 ~/sc_eval/scripts/sc_paga.py %s', getwd()), intern = FALSE)
      dim_res <- read.table(sprintf('%s/UMAP_Paga.txt', getwd()), header = FALSE, sep = '\t')
    }else if(method=='dhaka'){
      message('\n*** Dhaka ***\n')
      out_file <- sprintf('%s/TempExpr.txt', getwd())
      dhaka_out <- sprintf('%s/TempExprDhaka', getwd())
      write.table(t(as.matrix(expr)),
                  file = out_file, 
                  col.names = FALSE, 
                  row.names = FALSE,
                  sep = '\t',
                  append = FALSE,
                  quote = FALSE)
      local_wd <- getwd()
      local_tag <- TRUE
      count <- 0
      while(local_tag && count < 5){
        setwd('~/Dhaka/autoencoder/')
        system(command = sprintf('python3 -c \'import sys; from Dhaka import Dhaka; Dhaka(input_datafile=sys.argv[1], output_datafile=sys.argv[2], epochs=10, batch_size=50, N_starts=1, to_plot=0, verbose=False)\' %s %s', out_file, dhaka_out),
               intern = FALSE)
        setwd(local_wd)
        try({
          dim_res <- read.table(paste0(dhaka_out, '.txt'), header = FALSE, sep ='')
          if(sum(apply(dim_res, 1, function(x){any(is.na(x))})) == 0){
            local_tag <- FALSE
          }else{
            count <- count + 1
            cat(sprintf('Dhaka Retry\n'))
          }
        })
      }
    }else if(method=='dm'){
      require(destiny)
      message('\n*** Diffusion Map ***\n')
      dm_res <- DiffusionMap(t(as.matrix(expr)), density_norm = FALSE, n_eigs=2)
      dim_res <- dm_res@eigenvectors[,1:2]
    }
  }))[1]), file = 'RunTime.txt', append = TRUE)
  return(dim_res)
}

sc_impute <- function(count_mat = NULL, method = 'scimpute'){
  if(is.null(count_mat))
    stop('Provide count matrix for imputation')
  
  if(method=='scimpute'){
    suppressMessages(require(scImpute))
    suppressMessages(require(scran))
    message('\n*** ScImpute ***\n')
    if(!file.exists('scimputeOut/scimpute_count.rds')){
      saveRDS(as.matrix(count_mat), file = 'TempCounts.rds')
      message('Clustering...')
      local_clusters <- quickCluster(count_mat, 
                                     min.size = round(0.1*ncol(count_mat)), 
                                     method = 'igraph', 
                                     d = NA, 
                                     use.ranks = TRUE, 
                                     min.mean = 1)
      cat(sprintf('ScImpute\t%.4f\n', as.vector(system.time({
        scimpute(count_path = 'TempCounts.rds',
                 type = 'count',
                 infile = 'rds',
                 outfile = 'rds',
                 out_dir = 'scimputeOut/',
                 labeled = TRUE,
                 labels = paste0('Cluster-', local_clusters),
                 drop_thre = 0.5,
                 ncores = 20)
      }))[1]), file = 'RunTime.txt', append = TRUE)
    }
    imputed_counts <- readRDS('scimputeOut/scimpute_count.rds')
    return(Matrix(round(as.matrix(imputed_counts))))
  }
}

sc_traj <- function(obj = NULL, method = 'slingshot'){
  if(is.null(obj))
    stop('Provide data')
  if(method=='slingshot'){
    suppressMessages(require(dbscan))
    suppressMessages(require(scran))
    suppressMessages(require(scater))
    suppressMessages(require(slingshot))
    
    message('\n*** Slingshot ***\n')
    for(dim_meth in names(redDim(obj))){
      message(sprintf('[Slingshot: %s]', dim_meth))
      n_dim <- ncol(redDim(obj, dim_meth))
      local_dim_res <- cbind(redDim(obj, dim_meth), annot(obj))
      
      if(grepl(pattern='raw', dim_meth)){
        count_mat <- counts(obj, 'raw')
      }else{
        count_mat <- counts(obj, 'scimpute')
      }
      local_clusters <- quickCluster(count_mat, 
                                     min.size = round(0.1*ncol(count_mat)), 
                                     method = 'igraph',
                                     use.ranks = TRUE, 
                                     min.mean = 1)
      clust_labs <- setNames(paste0('Cluster-', local_clusters), colnames(count_mat))
      
      ## Remove outlier
      clust_labs <- clust_labs[!(clust_labs %in% 'Cluster-0')]
      local_dim_res <- local_dim_res[match(names(clust_labs), table = local_dim_res$Id),]
      
      ## Slingshot
      cat(sprintf('Slingshot\t%.4f\n', as.vector(system.time({
        lin_res <- getLineages(data = local_dim_res[,1:2], clusterLabels = clust_labs)
        curv_res <- getCurves(lin_res)
      }))[1]), file = 'RunTime.txt', append = TRUE)
      
      
      ## 3D/2D Plot
      curv_w <- sapply(curv_res@curves, function(x){
        x$w
      })

      all_curv <- list()
      for(i in 1:length(curv_res@curves)){
        x <- curv_res@curves[[i]]
        local_pos <- x$s
        local_ord <- x$ord
        local_pos <- local_pos[local_ord,]
        
        curv_dims <- as.data.frame(local_pos, row.names = NULL)
        curv_dims$color <- rep('#000000', nrow(curv_dims))
        curv_dims$type <- rep(paste0('Tree.', i), nrow(curv_dims))
        curv_dims$id <- rep('Root', nrow(curv_dims))
        curv_dims$label <- rep(NA, nrow(curv_dims))
        colnames(curv_dims) <- c(paste0('Dim.', 1:n_dim), 'Color', 'Type', 'Id', 'Label')
        all_curv[[i]] <- curv_dims
      }
      curv_dims <- do.call('rbind', all_curv)
      curv_dims <- cbind(curv_dims, as.data.frame(matrix(NA, ncol=ncol(slingPseudotime(curv_res)), nrow=nrow(curv_dims))))
      
      orig_dims <- as.data.frame(local_dim_res[,1:n_dim], make.names = FALSE)
      orig_dims$color <- local_dim_res$Color
      orig_dims$type <- rep('CB', nrow(orig_dims))
      orig_dims$id <- local_dim_res$Id
      orig_dims$label <- local_dim_res$Label
      colnames(orig_dims) <- c(paste0('Dim.', 1:n_dim), 'Color', 'Type', 'Id', 'Label')
      orig_dims <- cbind(orig_dims, slingPseudotime(curv_res))
      colnames(curv_dims) <- colnames(orig_dims)
      combined_data <- rbind(orig_dims, curv_dims)
      
      temp <- list(combined_data)
      names(temp) <- sprintf('%s_slingshot', dim_meth)
      traject(obj) <- temp
      saveRDS(temp, file = sprintf('%s_%s_Slingshot.rds', obj@title, dim_meth))
    }
  }else if(method=='monocle'){
    #suppressMessages(require(densityClust))
    #suppressMessages(require(dbscan))
    #suppressMessages(require(monocle))
    suppressMessages(require(gmodels))
    suppressMessages(require(DDRTree))
    
    message('\n*** Monocle 2 ***\n')
    for(norm_name in names(norm(obj))){
      norm_mat <- as.matrix(norm(obj, norm_name))
      annot_ft <- annot(obj)
      pca_res <- fast.prcomp(t(norm_mat))
      red_dim <- DDRTree(t(pca_res$x[,1:50]), k = 2, maxIter = 20)
      tree_coord <- t(red_dim$Y)
      data_coord <- t(red_dim$Z)
      
      data_coord <- data.frame('Dim.1' = data_coord[,1],
                               'Dim.2' = data_coord[,2],
                               #'Dim.3' = data_coord[,3],
                               'Id' = colnames(norm_mat),
                               'Label' = annot_ft$Label[match(colnames(norm_mat), table=annot_ft$Id)],
                               'Color' = annot_ft$Color[match(colnames(norm_mat), table = annot_ft$Id)],
                               stringsAsFactors = FALSE)
      
      tree_coord <- data.frame('Dim.1' = tree_coord[,1],
                               'Dim.2' = tree_coord[,2],
                               #'Dim.3' = tree_coord[,3],
                               'Id' = paste0('Tree.', 1:nrow(tree_coord)),
                               'Label' = rep(NA, nrow(tree_coord)),
                               'Color' = rep('#000000', nrow(tree_coord)),
                               stringsAsFactors = FALSE)
      
      combined_data <- rbind(data_coord, tree_coord)
      
      temp <- list(combined_data)
      names(temp) <- sprintf('%s_monocle', norm_name)
      traject(obj) <- temp
      saveRDS(temp, file = sprintf('%s_%s_Monocle2.rds', obj@title, norm_name))
    }
  }else if(method == 'OT'){
    message('\n*** Waddington-OT ***\n')
    for(norm_name in names(norm(obj))){
      expr_mat <- as.matrix(norm(obj, norm_name))
      write.table(t(expr_mat),
                  file = 'NormMat.txt',
                  sep = '\t',
                  col.names = TRUE,
                  row.names = TRUE,
                  append = FALSE,
                  quote = FALSE)
      
      annot_data <- annot(obj)
      annot_data <- data.frame('id' = annot_data$Id,
                               'day' = annot_data$Days)
      write.table(annot_data, 
                  file = 'cell_days.txt', 
                  sep = '\t', 
                  col.names = TRUE, 
                  row.names = FALSE,
                  append = FALSE,
                  quote = FALSE)
      system(command = sprintf('wot optimal_transport --matrix NormMat.txt --cell_days cell_days.txt --out %s_%s_WOT', obj@title, norm_name), intern = FALSE)
    }
  }
  return(obj)
}

