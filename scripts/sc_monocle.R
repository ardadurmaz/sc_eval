library(Matrix)
library(scran)
library(scater)
library(fastcluster)
library(Rtsne)
library(monocle)

all.labels <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1',
                'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3')
common.colors <- setNames(ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all.labels))],
                          all.labels)

for(data_name in c('alevin', 'cellranger')){
  for(base_name in c('raw', 'scimpute')){
    for(drug in c('alectinib', 'lorlatinib', 'crizotinib')){
      message(sprintf('Processing: [%s, %s, %s]', data_name, base_name, drug))
      if(drug == 'alectinib'){
        setwd(sprintf('~/sc_eval/data/%s_data/%s_res/deconv_res/alec/', data_name, base_name))
        sc.data <- readRDS(sprintf('%s_%s_%s_deconv.rds', data_name, drug, base_name))
      }else if(drug == 'lorlatinib'){
        setwd(sprintf('~/sc_eval/data/%s_data/%s_res/deconv_res/lor/', data_name, base_name))
        sc.data <- readRDS(sprintf('%s_%s_%s_deconv.rds', data_name, drug, base_name))
      }else{
        setwd(sprintf('~/sc_eval/data/%s_data/%s_res/deconv_res/criz/', data_name, base_name))
        sc.data <- readRDS(sprintf('%s_%s_%s_deconv.rds', data_name, drug, base_name))
      }
      count_mat <- round(sc.data$counts)
      labels <- sc.data$labels
      
      message('Running t-SNE...')
      tsne_res <- Rtsne(as.matrix(t(sc.data$norm_mat)),
                        perplexity = 30,
                        theta = 0.1,
                        partial_pca = TRUE,
                        pca_center = TRUE,
                        pca_scale = TRUE,
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
                               'Label' = labels[match(rownames(data_coord), table = colnames(count_mat))],
                               'Color' = common.colors[match(labels[match(rownames(data_coord), table = colnames(count_mat))],
                                                             table = names(common.colors))],
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
                  file = sprintf('~/sc_eval/data/%s_data/%s_res/monocle_res/%s/%s_%s_monocle.txt', data_name, base_name, drug, data_name, drug),
                  sep = '\t', 
                  col.names = TRUE, 
                  row.names = FALSE, 
                  append = FALSE,
                  quote = FALSE)
      
      # Plot
      system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                     sprintf('~/sc_eval/data/%s_data/%s_res/monocle_res/%s/%s_%s_monocle.txt', data_name, base_name, drug, data_name, drug),
                     sprintf('~/sc_eval/data/%s_data/%s_res/monocle_res/%s/%s_%s_monocle.pdf', data_name, base_name, drug, data_name, drug)))
      gc()
    }
  }
}
