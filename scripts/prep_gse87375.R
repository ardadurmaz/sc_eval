library(scater)
library(scran)
library(Matrix)

setwd('data/gse87375/')
count_mat <- read.table('GSE87375_Single_Cell_RNA-seq_Gene_Read_Count.txt', 
                        header = TRUE, 
                        sep = '\t', 
                        stringsAsFactors = FALSE)
temp_genes <- count_mat$ID
count_mat <- round(data.matrix(count_mat[,-c(1,2,3)]))
rownames(count_mat) <- temp_genes
rm(temp_genes)

# Filter Mat
count_mat <- count_mat[!apply(count_mat,1,function(x){all(x==0)}),]

# Spike In #
is.spike <- grep("ENSMUSG", invert = TRUE, rownames(count_mat)) 
sce <- SingleCellExperiment(list(counts = count_mat))
isSpike(sce, "spike") <- is.spike
sce <- computeSpikeFactors(sce, general.use = TRUE)
sce <- normalize(sce)
norm_mat <- logcounts(sce)
is.spike <- grep("ENSMUSG", invert = TRUE, rownames(norm_mat))
norm_mat <- norm_mat[-is.spike,]

# Deconvolution #
count_mat <- count_mat[grep("ENSMUSG", rownames(count_mat)),]
sce <- SingleCellExperiment(list(counts = count_mat))
clusters <- quickCluster(sce, min.size=10)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- normalize(sce)
norm_mat <- logcounts(sce)

# Sctransform #
require(sctransform)
count_mat <- count_mat[grep("ENSMUSG", rownames(count_mat)),]
vst_out <- sctransform::vst(as.matrix(count_mat), 
                            latent_var = c('log_umi'), 
                            return_gene_attr = TRUE, 
                            return_cell_attr = TRUE)
norm_mat <- vst_out$y
norm_mat <- norm_mat[!vst_out$model_pars_outliers,]

annot <- data.frame('id' = colnames(norm_mat),
                    'type' = substr(colnames(norm_mat), start = 1, stop = 1),
                    'time' = as.character(sapply(substr(colnames(norm_mat), start = 2, stop = 100000), function(s){unlist(strsplit(s, split = '_'))[1]})))


saveRDS(list('counts' = count_mat,
             'labels' = annot,
             'norm_mat' = norm_mat), 
        file = 'gse87375_sctransform_data.rds')
