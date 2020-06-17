library(Matrix)

for(data_name in c('alevin', 'cellranger')){
  for(base_name in c('raw', 'scimpute')){
    for(norm_name in c('sctransform', 'deconv')){
      
      for(drug in c('alectinib', 'lorlatinib', 'crizotinib')){
        message(sprintf('Processing: [%s, %s, %s, %s, %s]', data_name, base_name, norm_name, 'dhaka', drug))
        
        if(drug == 'alectinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/alec/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          write.table(round(as.matrix(t(sc_data$norm_mat)), digits=3), 
                      file = sprintf(sprintf('dhaka_res/%s_%s_%s_%s_mat.txt', data_name, drug, base_name, norm_name)),
                      sep = '\t',
                      col.names = FALSE,
                      row.names = FALSE,
                      append = FALSE,
                      quote = FALSE)
        }else if(drug == 'lorlatinib'){
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/lor/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          write.table(round(as.matrix(t(sc_data$norm_mat)), digits=3), 
                      file = sprintf(sprintf('dhaka_res/%s_%s_%s_%s_mat.txt', data_name, drug, base_name, norm_name)),
                      sep = '\t',
                      col.names = FALSE,
                      row.names = FALSE,
                      append = FALSE,
                      quote = FALSE)
        }else{
          setwd(sprintf('~/sc_eval/data/%s_data/%s_res/%s_res/criz/', data_name, base_name, norm_name))
          sc_data <- readRDS(sprintf('%s_%s_%s_%s.rds', data_name, drug, base_name, norm_name))
          write.table(round(as.matrix(t(sc_data$norm_mat)), digits=3), 
                      file = sprintf(sprintf('dhaka_res/%s_%s_%s_%s_mat.txt', data_name, drug, base_name, norm_name)),
                      sep = '\t',
                      col.names = FALSE,
                      row.names = FALSE,
                      append = FALSE,
                      quote = FALSE)
        }
      }
    }
  }
}
