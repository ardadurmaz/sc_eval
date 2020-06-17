library(stringi)
setwd('~/sc_eval/data/gse63818/')

expr_data <- read.table('gse63818_exp.tsv', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
expr_mat <- data.matrix(expr_data[,-1])
rownames(expr_mat) <- expr_data[,1]

expr_mat_pgc <- expr_mat[,grep(pattern = '_PGC_', colnames(expr_mat))]
annot_data <- data.frame('id' = colnames(expr_mat_pgc),
                         'sex' = sapply(stri_split(colnames(expr_mat_pgc), fixed = '_'), function(x){x[1]}),
                         'time' = sapply(stri_split(colnames(expr_mat_pgc), fixed = '_'), function(x){x[3]}))

idx <- apply(expr_mat_pgc, 1, function(x){sum(x > 0.1) > 3})
expr_mat_pgc <- expr_mat_pgc[idx,]

all.data <- list('labels' = annot_data,
                 'norm_mat' = log2(expr_mat_pgc + 0.25))

saveRDS(all.data, file = '~/sc_eval/data/gse63818/gse63818_norm.rds')
