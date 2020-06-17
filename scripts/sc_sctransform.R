library(Matrix)
library(sctransform)

# Imputed Data
sc.mat <- readRDS('scimpute_crizotinib_out/scimpute_count.rds')
sc.out.idx <- readRDS('data/imputation.out.crizotinib.idx.rds')

# Raw Data
#sc.mat <- readRDS('data/sc.counts.filtered.crizotinib.rds')

vst_out <- sctransform::vst(round(as.matrix(sc.mat)), 
                            latent_var = c('log_umi'), 
                            return_gene_attr = TRUE, 
                            return_cell_attr = TRUE, 
                            show_progress = TRUE, 
                            n_genes = NULL)
norm.mat <- vst_out$y
norm.mat <- norm.mat[!vst_out$model_pars_outliers,]
saveRDS(norm.mat, file = 'data/sc.mat.scimpute.sctransform.crizotinib.rds')

#rs.var <- vst_out$gene_attr$residual_variance[!vst_out$model_pars_outliers]

# Plot residual variance
# require(ggrepel)
# plot.temp <- data.frame('id' = rownames(norm.mat),
#                         'var' = rs.var,
#                         'symbol'= sc.data$features$V2[match(rownames(norm.mat), table = sc.data$features$V1)],
#                         'alpha' = rs.var / max(rs.var),
#                         stringsAsFactors = FALSE)
# 
# 
# 
# p <- ggplot(plot.temp, aes(x = id, y = var)) +
#   geom_point(aes(alpha = alpha), size = 1) +
#   theme_minimal() +
#   ggtitle('All Drugs') +
#   theme_classic() +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')) +
#   geom_label_repel(data = subset(plot.temp, plot.temp$var > 5),
#                    aes(label = symbol), size = 4) +
#   guides(alpha = FALSE, 
#          color = FALSE)
# ggsave(p, filename = 'plots/residual_variance.pdf', width = 12, height = 10)

# Save
#sc.data$sc_transform <- norm.mat
#saveRDS(sc.data, file = 'data/sc.data.sctransform.rds')

