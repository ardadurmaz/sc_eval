setwd('~/sc_eval/')
runtime_files <- list.files(pattern='runtime.+txt$',ignore.case = TRUE,recursive = TRUE)
runtime_data <- sapply(runtime_files,simplify = FALSE,function(f){
  read.table(f)
})
runtime_data <- do.call('rbind', runtime_data)
colnames(runtime_data) <- c('Tool', 'Time')
runtime_data$Data <- as.character(sapply(rownames(runtime_data), function(s){unlist(strsplit(s, split = '/', fixed = TRUE))[2]}))

runtime_data_filt <- droplevels(subset(runtime_data, runtime_data$Tool == 'Monocle' | runtime_data$Tool == 'Slingshot'))

library(ggplot2)
library(ggsci)

ggplot(runtime_data_filt, aes(x = Tool, y = log10(Time), fill = Data)) +
  geom_boxplot() +
  xlab('') + 
  ylab('log10(Time)') +
  theme_classic() +
  theme(axis.text = element_text(size = 12, face = 'bold')) +
  scale_fill_rickandmorty(labels = c('Alevin_Alectinib' = 'Alectinib (Alevin)',
                                     'Alevin_Crizotinib' = 'Crizotinib (Alevin)',
                                     'Alevin_Lorlatinib' = 'Lorlatinib (Alevin)',
                                     'Cellranger_Alectinib' = 'Alectinib (Cellranger)',
                                     'Cellranger_Crizotinib' = 'Crizotinib (Cellranger)',
                                     'Cellranger_Lorlatinib' = 'Lorlatinib (Cellranger)',
                                     'EGEOD103334' = 'GSE103334',
                                     'EMTAB3929' = 'E-MTAB-3929',
                                     'GSE87375_A' = 'GSE87375 A',
                                     'GSE87375_B' = 'GSE87375 B'))
