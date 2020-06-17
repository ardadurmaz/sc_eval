library(Matrix)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

title_sub <- 'Alectinib'
for(i in 1:5){
  title <- paste0(title_sub, sprintf('_S%d', i))
  setwd(sprintf('~/sc_eval/data/%s', title))
  proc_data <- readRDS(sprintf('%s_ProcessedData.rds', title))
  
  for(norm_name in names(norm(proc_data))){
    expr_mat <- as.matrix(norm(proc_data, norm_name))
    write.table(t(expr_mat),
                file = 'WOT_NormMat.txt',
                sep = '\t',
                col.names = TRUE,
                row.names = TRUE,
                append = FALSE,
                quote = FALSE)
    
    annot_data <- annot(proc_data)
    annot_data <- data.frame('id' = annot_data$Id,
                             'day' = annot_data$Days)
    write.table(annot_data, 
                file = 'WOT_CellDays.txt', 
                sep = '\t', 
                col.names = TRUE, 
                row.names = FALSE,
                append = FALSE,
                quote = FALSE)
    system(command = 'wot optimal_transport --matrix NormMat.txt --cell_days CellDays.txt --out wot_res', intern = FALSE)
  }
  
}




