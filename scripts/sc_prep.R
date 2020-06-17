require(Matrix)
source('~/sc_eval/scripts/sc_utils.R')

setwd('~/sc_eval/data/')

## Data 1 ##
cellranger_data <- readRDS('cellranger_data/cellranger_data.rds')
all.labels <- c('H3122-1', 
                'H3122-2', 
                'H3122-Alec-4h', 
                'H3122-Alec-48h', 
                'H3122-Alec-3W', 
                'H3122-erAlec-1', 
                'H3122-GFP-erAlec',
                'H3122-Lor-48h', 
                'H3122-Lor-3W', 
                'H3122-erLor', 
                'H3122-erLor-1',
                'H3122-Cz-48h', 
                'H3122-Cz-3W', 
                'H3122-erCz', 
                'H3122-erCz-3')

common_colors <- setNames(ggsci::pal_simpsons('springfield')(15)[as.numeric(as.factor(all.labels))], all.labels)
days_map <- c('H3122-1' = 0, 
              'H3122-2' = 0, 
              'H3122-Alec-3W' = 21, 
              'H3122-Alec-48h'= 2, 
              'H3122-Alec-4h' = 0.15, 
              'H3122-erAlec-1' = 112, 
              'H3122-GFP-erAlec' = 112,
              'H3122-Lor-3W' = 21, 
              'H3122-Lor-48h'= 2, 
              'H3122-erLor-1' = 112, 
              'H3122-erLor' = 112,
              'H3122-Cz-3W' = 21, 
              'H3122-Cz-48h'= 2, 
              'H3122-erCz-3' = 112, 
              'H3122-erCz' = 112)

data_1 <- list()
for(d in c('Alectinib', 'Lorlatinib', 'Crizotinib')){
  if(d == 'Alectinib'){
    local.labels <- c('H3122-1', 
                      'H3122-2', 
                      'H3122-Alec-4h', 
                      'H3122-Alec-48h', 
                      'H3122-Alec-3W', 
                      'H3122-erAlec-1', 
                      'H3122-GFP-erAlec')
  }else if(d == 'Lorlatinib'){
    local.labels <- c('H3122-1', 
                      'H3122-2',  
                      'H3122-Lor-48h', 
                      'H3122-Lor-3W', 
                      'H3122-erLor-1', 
                      'H3122-erLor')
  }else{
    local.labels <- c('H3122-1', 
                      'H3122-2',  
                      'H3122-Cz-48h', 
                      'H3122-Cz-3W', 
                      'H3122-erCz', 
                      'H3122-erCz-3')
  }
  local.idx <- cellranger_data$labels %in% local.labels
  local_counts <- cellranger_data$counts[,local.idx]
  row_idx <- rowSums(local_counts > 0) > 0
  local_counts <- local_counts[row_idx,]
  col_idx <- colSums(local_counts > 0) > 0
  local_counts <- local_counts[,col_idx]

  ## Levels ##
  test_thr_col <- seq(-3, 0, 0.5)
  test_thr_row <- seq(4, 7, 0.5)
  for(i in 1:length(test_thr_col)){
    local_local_counts <- local_counts
    col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
    local_local_counts <- local_local_counts[,as.vector(col_idx)]
    row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
    local_local_counts <- local_local_counts[as.vector(row_idx),]
    
    local_labels <- cellranger_data$labels[match(colnames(local_local_counts), 
                                                 table = colnames(cellranger_data$counts))]
    annot_data <- data.frame('Id' = colnames(local_local_counts),
                             'Label' = local_labels,
                             'Days' = days_map[match(local_labels, table=names(days_map))],
                             'Color' = common_colors[match(local_labels, 
                                                           table = names(common_colors))])
    data_1[[length(data_1) + 1]] <- new('SCEVAL', 
                                        title = sprintf("%s_S%d", d, i), 
                                        counts = list('raw' = Matrix(local_local_counts)), 
                                        annot = annot_data)
  }
  rm(local_counts, local.idx, row_idx)
  gc()
}
saveRDS(data_1, file = 'CombinedData.rds')

#################################
## Data 2 ##
all_labels <- c('embryonic day 3', 'embryonic day 4', 'embryonic day 5', 'embryonic day 6', 'embryonic day 7')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 5, 'Set2')[as.numeric(as.factor(all_labels))], all_labels)

days_map <- c('embryonic day 3' = 3,
              'embryonic day 4' = 4,
              'embryonic day 5' = 5,
              'embryonic day 6' = 6,
              'embryonic day 7' = 7)

setwd('~/sc_eval/data/E-MTAB-3929/')

count_mat <- data.matrix(read.table('counts.txt', sep = '\t', check.names = FALSE))
annot <- read.table('E-MTAB-3929.sdrf.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '', check.names = FALSE)
annot_data <- data.frame('Id' = annot$`Source Name`,
                         'Label' = annot$`Characteristics[developmental stage]`)

annot_data <- data.frame('Id' = annot_data$Id,
                         'Label' = annot_data$Label,
                         'Days' = days_map[match(annot_data$Label, table = names(days_map))],
                         'Color' = common_colors[match(annot_data$Label, table = names(common_colors))])


data_2 <- list()
test_thr_col <- seq(-3, 0, 0.5)
test_thr_row <- seq(3, 6, 0.5)

for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_2[[length(data_2) + 1]] <- new('SCEVAL', 
                                      title = sprintf("%s_S%d", 'EMTAB3929', i), 
                                      counts = list('raw' = Matrix(local_local_counts)), 
                                      annot = local_annot_data)
}
saveRDS(data_2, file='EMTAB3929_Data.rds')

## Data 3 ##
all_labels <- c('E17.5', 'P0', 'P3', 'P9', 'P15', 'P18', 'P60')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 7, 'Dark2'), all_labels)

days_map <- c('E17.5' = 0,
              'P0' = 4,
              'P3' = 7,
              'P9' = 13,
              'P15' = 19,
              'P18' = 22,
              'P60' = 64)

### Type A ###
sc_data <- readRDS('gse87375/raw_res/deconv_res/gse87375_a_raw_deconv.rds')
count_mat <- round(as.matrix(sc_data$counts))
annot_data <- data.frame('Id' = colnames(count_mat),
                         'Label' = sc_data$labels$time,
                         'Days' = days_map[match(sc_data$labels$time, table = names(days_map))],
                         'Color' = common_colors[match(sc_data$labels$time, table = names(common_colors))])

data_3_a <- list()
test_thr_col <- seq(-4, -1, 0.5)
test_thr_row <- seq(2, 5, 0.5)
for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_3_a[[length(data_3_a) + 1]] <- new('SCEVAL', 
                                      title = sprintf("%s_S%d", 'GSE87375_A', i), 
                                      counts = list('raw' = Matrix(local_local_counts)), 
                                      annot = local_annot_data)
}
saveRDS(data_3_a, file = 'GSE87375_A_Data.rds')


### Type B ###
sc_data <- readRDS('gse87375/raw_res/deconv_res/gse87375_b_raw_deconv.rds')
count_mat <- round(as.matrix(sc_data$counts))
annot_data <- data.frame('Id' = colnames(count_mat),
                         'Label' = sc_data$labels$time,
                         'Days' = days_map[match(sc_data$labels$time, table = names(days_map))],
                         'Color' = common_colors[match(sc_data$labels$time, table = names(common_colors))])

data_3_b <- list()
test_thr_col <- seq(-4, -1, 0.5)
test_thr_row <- seq(2, 5, 0.5)
for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_3_b[[length(data_3_b) + 1]] <- new('SCEVAL', 
                                          title = sprintf("%s_S%d", 'GSE87375_B', i), 
                                          counts = list('raw' = Matrix(local_local_counts)), 
                                          annot = local_annot_data)
}
saveRDS(data_3_b, file = 'GSE87375_B_Data.rds')


## Data 4 ##
all_labels <- c('Week-0', 'Week-1', 'Week-2', 'Week-6')
common_colors <- setNames(RColorBrewer::brewer.pal(n = 4, 'Dark2'), all_labels)
days_map <- c('Week-0' = 0,
              'Week-1' = 7,
              'Week-2' = 14,
              'Week-6' = 42)

sc_data <- readRDS('E-GEOD-103334/egeod103334_raw.rds')
count_mat <- round(as.matrix(sc_data$counts))
annot_data <- data.frame('Id' = colnames(count_mat),
                         'Label' = sc_data$labels,
                         'Days' = days_map[match(sc_data$labels, table = names(days_map))],
                         'Color' = common_colors[match(sc_data$labels, table = names(common_colors))])
data_4 <- list()
test_thr_col <- seq(-4, -1, 0.5)
test_thr_row <- seq(2, 5, 0.5)

for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_4[[length(data_4) + 1]] <- new('SCEVAL', 
                                          title = sprintf("%s_S%d", 'EGEOD103334', i), 
                                          counts = list('raw' = Matrix(local_local_counts)), 
                                          annot = local_annot_data)
}
saveRDS(data_4, file = 'EGEOD103334_Data.rds')

## Data 5 ##
count_mat <- readRDS('MCF7_Counts.rds')
sample_data <- read.csv('MCF7_SraRunTable.txt', stringsAsFactors = FALSE)
idx <- match(colnames(count_mat), table = sample_data$Run)

annot_data <- data.frame('Id' = colnames(count_mat),
                         'Label' = gsub(pattern = '^.+?_', replacement = '', perl = TRUE, x = sample_data$Library.Name)[idx],
                         stringsAsFactors = FALSE)
bulk_idx <- grepl(pattern = "0k", annot_data$Label)
annot_data <- annot_data[!bulk_idx,]
count_mat <- count_mat[,match(annot_data$Id, table = colnames(count_mat))]

count_mat <- count_mat[rowSums(count_mat > 0) > 0,]
annot_data$Days <- as.numeric(gsub(pattern = 'h', replacement = '', annot_data$Label))
annot_data$Color <- RColorBrewer::brewer.pal(4, 'Dark2')[as.numeric(as.factor(annot_data$Label))]

data_5 <- list()
test_thr_col <- seq(-3, -0, 0.5)
test_thr_row <- seq(1, 4, 0.5)

for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_5[[length(data_5) + 1]] <- new('SCEVAL', 
                                      title = sprintf("%s_S%d", 'GSE107858', i), 
                                      counts = list('raw' = Matrix(local_local_counts)), 
                                      annot = local_annot_data)
}
saveRDS(data_5, file = 'GSE107858_Data.rds')

## Data 6 ##
count_mat <- readRDS('T47D_Counts.rds')
sample_data <- read.csv('T47D_SraRunTable.txt', stringsAsFactors = FALSE)
idx <- match(colnames(count_mat), table = sample_data$Run)

annot_data <- data.frame('Id' = colnames(count_mat),
                         'Label' = gsub(pattern = '^.+?_', replacement = '', perl = TRUE, x = sample_data$Library.Name)[idx])
bulk_idx <- grepl(pattern = "population", annot_data$Label)
annot_data <- annot_data[!bulk_idx,]
count_mat <- count_mat[,match(annot_data$Id, table = colnames(count_mat))]
annot_data$Label <- gsub(pattern = '_.+$', replacement = '', perl = TRUE, x = annot_data$Label)
count_mat <- count_mat[rowSums(count_mat > 0) > 0,]

annot_data$Days <- as.numeric(gsub(pattern = 'h', replacement = '', annot_data$Label))
annot_data$Color <- RColorBrewer::brewer.pal(4, 'Dark2')[as.numeric(as.factor(annot_data$Label))]

data_6 <- list()
test_thr_col <- seq(-3, -0, 0.5)
test_thr_row <- seq(1, 4, 0.5)

for(i in 1:length(test_thr_col)){
  local_local_counts <- count_mat
  col_idx <- scale(log(colSums(local_local_counts > 0))) > test_thr_col[i]
  local_local_counts <- local_local_counts[,as.vector(col_idx)]
  row_idx <- log(rowSums(local_local_counts > 0)) > test_thr_row[i]
  local_local_counts <- local_local_counts[as.vector(row_idx),]
  
  local_annot_data <- annot_data[match(colnames(local_local_counts), table=annot_data$Id),]
  data_6[[length(data_6) + 1]] <- new('SCEVAL', 
                                      title = sprintf("%s_S%d", 'GSE107863', i), 
                                      counts = list('raw' = Matrix(local_local_counts)), 
                                      annot = local_annot_data)
}
saveRDS(data_6, file = 'GSE107863_Data.rds')
