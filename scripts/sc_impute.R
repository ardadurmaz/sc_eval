library(Matrix)
library(scran)
library(scImpute)

for(drug in c('alectinib', 'lorlatinib', 'crizotinib')){
  message('Processing drug: ', drug)
  if(drug == 'alectinib'){
    sc.data <- readRDS('~/sc_eval/data/cellranger_data/raw_res/deconv_res/alec/cellranger_alectinib_deconv.rds')
    count.mat <- as.matrix(sc.data$counts)
    saveRDS(count.mat, file = '~/sc_eval/data/temp_data/temp_count.rds')
    
    message('Clustering...')
    clusters <- quickCluster(count.mat, 
                             min.size = 30, 
                             method = 'igraph', 
                             use.ranks = FALSE)
    
    message('Imputation...')
    out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                        type = 'count',
                        infile = 'rds',
                        outfile = 'rds',
                        out_dir = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_alectinib_out',
                        labeled = TRUE,
                        labels = paste0('Cluster-', clusters),
                        drop_thre = 0.3,
                        ncores = 20)
    saveRDS(out.idx, file = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_alectinib_outidx.rds')
  }else if(drug == 'lorlatinib'){
    sc.data <- readRDS('~/sc_eval/data/cellranger_data/raw_res/deconv_res/lor/cellranger_lorlatinib_deconv.rds')
    count.mat <- as.matrix(sc.data$counts)
    saveRDS(count.mat, file = '~/sc_eval/data/temp_data/temp_count.rds')
    
    message('Clustering...')
    clusters <- quickCluster(count.mat, 
                             min.size = 30, 
                             method = 'igraph', 
                             use.ranks = FALSE)
    
    message('Imputation...')
    out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                        type = 'count',
                        infile = 'rds',
                        outfile = 'rds',
                        out_dir = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_lorlatinib_out',
                        labeled = TRUE,
                        labels = paste0('Cluster-', clusters),
                        drop_thre = 0.3,
                        ncores = 20)
    saveRDS(out.idx, file = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_lorlatinib_outidx.rds')
    
  }else{
    sc.data <- readRDS('~/sc_eval/data/cellranger_data/raw_res/deconv_res/criz/cellranger_crizotinib_deconv.rds')
    count.mat <- as.matrix(sc.data$counts)
    saveRDS(count.mat, file = '~/sc_eval/data/temp_data/temp_count.rds')
    
    message('Clustering...')
    clusters <- quickCluster(count.mat, 
                             min.size = 30, 
                             method = 'igraph', 
                             use.ranks = FALSE)
    
    message('Imputation...')
    out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                        type = 'count',
                        infile = 'rds',
                        outfile = 'rds',
                        out_dir = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_crizotinib_out',
                        labeled = TRUE,
                        labels = paste0('Cluster-', clusters),
                        drop_thre = 0.3,
                        ncores = 20)
    saveRDS(out.idx, file = '~/sc_eval/data/cellranger_data/scimpute_res/scimpute_crizotinib_outidx.rds')
  }
}
