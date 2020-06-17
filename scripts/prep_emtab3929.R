library(Matrix)

setwd('~/sc_eval/data/E-MTAB-3929/')

count_mat <- data.matrix(read.table('counts.txt', sep = '\t', check.names = FALSE))
annot <- read.table('E-MTAB-3929.sdrf.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '', check.names = FALSE)
annot_data <- data.frame('Id' = annot$`Source Name`,
                         'Label' = annot$`Characteristics[developmental stage]`)

## Filter ##
idx <- rowSums(count_mat > 0) > 30
count_mat <- count_mat[idx,]

## Deconvolution ##
require(scran)
sce <- SingleCellExperiment(list(counts = as.matrix(count_mat)))
local.clusters <- quickCluster(sce, 
                               min.size = 20, 
                               method = 'igraph', 
                               d = NA, 
                               use.ranks = TRUE, 
                               min.mean = 1)
sce <- computeSumFactors(sce, clusters = local.clusters)

require(scater)
sce <- normalize(sce)
norm_mat <- logcounts(sce)
saveRDS(list('counts' = count_mat,
             'labels' = annot_data,
             'norm_mat' = norm_mat), 
        file = '~/sc_eval/data/emtab3929/embtab3929_raw_deconv.rds')

## Sctransform ##
require(sctransform)
vst_out <- sctransform::vst(count_mat, 
                            latent_var = c('log_umi'), 
                            return_gene_attr = FALSE, 
                            return_cell_attr = TRUE)
norm.mat <- vst_out$y
norm.mat <- norm.mat[!vst_out$model_pars_outliers,]


# Save Results
saveRDS(list('counts' = count_mat,
             'labels' = annot_data,
             'norm_mat' = norm.mat), 
        file = '~/sc_eval/data/emtab3929/embtab3929_raw_sctransform.rds')


## Scimpute ##
require(scImpute)
saveRDS(count_mat, file = '~/sc_eval/data/temp_data/temp_count.rds')

message('Clustering...')
local.clusters <- quickCluster(count_mat, 
                               min.size = 20, 
                               method = 'hclust', 
                               d = NA, 
                               use.ranks = TRUE, 
                               min.mean = 1)
message('Imputation...')
out.idx <- scimpute(count_path = '~/sc_eval/data/temp_data/temp_count.rds',
                    type = 'count',
                    infile = 'rds',
                    outfile = 'rds',
                    out_dir = '~/sc_eval/data/emtab3929/embtab3929_scimpute__out',
                    labeled = TRUE,
                    labels = paste0('Cluster-', local.clusters),
                    drop_thre = 0.3,
                    ncores = 20)

## Deconvolution Scimpute ##
count_mat <- readRDS('~/sc_eval/data/emtab3929/embtab3929_scimpute__outscimpute_count.rds')
annot <- read.table('E-MTAB-3929.sdrf.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '', check.names = FALSE)
annot_data <- data.frame('Id' = annot$`Source Name`,
                         'Label' = annot$`Characteristics[developmental stage]`)

sce <- SingleCellExperiment(list(counts = round(as.matrix(count_mat))))
local.clusters <- quickCluster(sce, 
                               min.size = 20, 
                               method = 'hclust', 
                               d = NA, 
                               use.ranks = TRUE, 
                               min.mean = 1)
sce <- computeSumFactors(sce, clusters = local.clusters)

require(scater)
sce <- normalize(sce)
norm_mat <- logcounts(sce)
saveRDS(list('counts' = count_mat,
             'labels' = annot_data,
             'norm_mat' = norm_mat), 
        file = '~/sc_eval/data/emtab3929/embtab3929_scimpute_deconv.rds')

## Sctransform Scimpute ##
count_mat <- readRDS('~/sc_eval/data/emtab3929/embtab3929_scimpute__outscimpute_count.rds')
annot <- read.table('E-MTAB-3929.sdrf.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '', check.names = FALSE)
annot_data <- data.frame('Id' = annot$`Source Name`,
                         'Label' = annot$`Characteristics[developmental stage]`)
vst_out <- sctransform::vst(round(count_mat), 
                            latent_var = c('log_umi'), 
                            return_gene_attr = FALSE, 
                            return_cell_attr = TRUE)
norm.mat <- vst_out$y
norm.mat <- norm.mat[!vst_out$model_pars_outliers,]


# Save Results
saveRDS(list('counts' = count_mat,
             'labels' = annot_data,
             'norm_mat' = norm.mat), 
        file = '~/sc_eval/data/emtab3929/embtab3929_scimpute_sctransform.rds')
