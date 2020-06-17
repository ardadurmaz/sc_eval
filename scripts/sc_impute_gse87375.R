library(Matrix)
library(scran)
library(scImpute)

for(data_name in c('a', 'b')){

  if(data_name == 'a'){
    sc_data <- readRDS('~/sc_eval/data/gse87375/raw_res/deconv_res/gse87375_raw_deconv_data.rds')
    label_data <- sc_data$labels
    label_data <- subset(label_data, label_data$type == data_name)
    count_mat <- sc_data$counts[,match(label_data$id, table = colnames(sc_data$counts))]
    
    saveRDS(count_mat, file = '~/sc_eval/data/temp_data/temp_count.rds')
    
    message('Clustering...')
    clusters <- quickCluster(count_mat, 
                             min.size = 10, 
                             method = 'igraph', 
                             use.ranks = FALSE)
    
    message('Imputation...')
    out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                        type = 'count',
                        infile = 'rds',
                        outfile = 'rds',
                        out_dir = '~/sc_eval/data/gse87375/scimpute_res/scimpute_a_out',
                        labeled = TRUE,
                        labels = paste0('Cluster-', clusters),
                        drop_thre = 0.3,
                        ncores = 20)
    saveRDS(out.idx, file = '~/sc_eval/data/gse87375/scimpute_res/scimpute_a_outidx.rds')
  }else if(data_name == 'b'){
    sc_data <- readRDS('~/sc_eval/data/gse87375/raw_res/deconv_res/gse87375_raw_deconv_data.rds')
    label_data <- sc_data$labels
    label_data <- subset(label_data, label_data$type == data_name)
    count_mat <- sc_data$counts[,match(label_data$id, table = colnames(sc_data$counts))]
    
    saveRDS(count_mat, file = '~/sc_eval/data/temp_data/temp_count.rds')
    
    message('Clustering...')
    clusters <- quickCluster(count_mat, 
                             min.size = 10, 
                             method = 'igraph', 
                             use.ranks = FALSE)
    
    message('Imputation...')
    out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                        type = 'count',
                        infile = 'rds',
                        outfile = 'rds',
                        out_dir = '~/sc_eval/data/gse87375/scimpute_res/scimpute_b_out',
                        labeled = TRUE,
                        labels = paste0('Cluster-', clusters),
                        drop_thre = 0.3,
                        ncores = 20)
    saveRDS(out.idx, file = '~/sc_eval/data/gse87375/scimpute_res/scimpute_b_outidx.rds')
  }
}