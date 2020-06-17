library(Matrix)

all.labels <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1',
                'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3')
common.colors <- setNames(ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all.labels))],
                          all.labels)


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
          write.table(embedd, 
                      file = sprintf('dhaka_res/%s_%s_%s_%s_dhaka_meta.txt', data_name, drug, base_name, norm_name),
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
          write.table(embedd, 
                      file = sprintf('dhaka_res/%s_%s_%s_%s_dhaka_meta.txt', data_name, drug, base_name, norm_name),
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
          write.table(embedd, 
                      file = sprintf('dhaka_res/%s_%s_%s_%s_dhaka_meta.txt', data_name, drug, base_name, norm_name),
                      sep = '\t',
                      col.names = TRUE,
                      row.names = FALSE,
                      append = FALSE,
                      quote = FALSE)
        }
        # Plot
        system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                       sprintf('dhaka_res/%s_%s_%s_%s_dhaka_meta.txt', data_name, drug, base_name, norm_name),
                       sprintf('dhaka_res/%s_%s_%s_%s_dhaka_meta.pdf', data_name, drug, base_name, norm_name)))
      }
    }
  }
}
