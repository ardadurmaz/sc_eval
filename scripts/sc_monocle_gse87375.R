library(Matrix)
library(scran)
library(scater)
library(fastcluster)
library(Rtsne)
library(monocle)

all_labels <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 7, 'Dark2'), all_labels)

for(data_name in c('a', 'b')){
  for(base_name in c('raw', 'scimpute')){
    setwd(sprintf('~/sc_eval/data/gse87375/%s_res/deconv_res', base_name))
    sc_data <- readRDS(sprintf('gse87375_%s_%s_deconv.rds', data_name, base_name))
    
    count_mat <- round(sc_data$counts)
    labels <- sc_data$labels
    
    message('Running t-SNE...')
    tsne_res <- Rtsne(as.matrix(t(sc_data$norm_mat)),
                      perplexity = 10,
                      theta = 0.1,
                      partial_pca = FALSE,
                      num_threads = 20)
    tsne_data <- tsne_res$Y
    rownames(tsne_data) <- colnames(count_mat)
    
    feat.data <- data.frame('gene_short_name' = rownames(count_mat),
                            row.names = rownames(count_mat))
    pheno.data <- data.frame('id' = colnames(count_mat),
                             'labels' = labels,
                             row.names = colnames(count_mat))
    
    pd <- new("AnnotatedDataFrame", 
              data = pheno.data)
    fd <- new("AnnotatedDataFrame",
              data = feat.data)
    
    cds <- newCellDataSet(count_mat, 
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    reducedDimA(cds) <- t(tsne_data)
    
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    cds <- clusterCells(cds, 
                        method = 'densityPeak', 
                        verbose = TRUE,
                        inspect_rho_sigma = FALSE)
    
    diff_test_res <- differentialGeneTest(cds, 
                                          fullModelFormulaStr = "~Cluster", 
                                          cores = 20)
    ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
    cds <- setOrderingFilter(cds, ordering_genes)
    
    cds <- reduceDimension(cds, 
                           max_components = 3, 
                           norm_method = 'log',
                           scaling = TRUE,
                           reduction_method = 'DDRTree')
    cds <- orderCells(cds)
    
    tree_coord <- t(reducedDimK(cds))
    data_coord <- t(reducedDimS(cds))
    
    data_coord <- data.frame('Dim.1' = data_coord[,1],
                             'Dim.2' = data_coord[,2],
                             'Dim.3' = data_coord[,3],
                             'Id' = colnames(cds),
                             'Label' = labels$time[match(rownames(data_coord), 
                                                         table = colnames(count_mat))],
                             'Color' = common_colors[match(labels$time[match(rownames(data_coord), 
                                                                             table = colnames(count_mat))], 
                                                           table = names(common_colors))],
                             stringsAsFactors = FALSE)
    
    tree_coord <- data.frame('Dim.1' = tree_coord[,1],
                             'Dim.2' = tree_coord[,2],
                             'Dim.3' = tree_coord[,3],
                             'Id' = paste0('Tree.', 1:nrow(tree_coord)),
                             'Label' = rep(NA, nrow(tree_coord)),
                             'Color' = rep('#000000', nrow(tree_coord)),
                             stringsAsFactors = FALSE)
    
    combined_data <- rbind(data_coord, tree_coord)
    write.table(combined_data, 
                file = sprintf('~/sc_eval/data/gse87375/%s_res/monocle_res/gse87375_%s_%s_monocle.txt', base_name, base_name, data_name),
                sep = '\t', 
                col.names = TRUE, 
                row.names = FALSE, 
                append = FALSE,
                quote = FALSE)
    
    # Plot
    system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                   sprintf('~/sc_eval/data/gse87375/%s_res/monocle_res/gse87375_%s_%s_monocle.txt', base_name, base_name, data_name),
                   sprintf('~/sc_eval/data/gse87375/%s_res/monocle_res/gse87375_%s_%s_monocle.pdf', base_name, base_name, data_name)))
    gc()
    
  }
}