library(Matrix)
library(biomaRt)

ppi <- read.table('data/stringdb/9606.protein.links.v11.0.txt', sep = '', header = TRUE, stringsAsFactors = FALSE)
ppi$protein1 <- gsub(pattern = '^9606.', replacement = '', ppi$protein1)
ppi$protein2 <- gsub(pattern = '^9606.', replacement = '', ppi$protein2)

ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_peptide_id', 'ensembl_gene_id'),
                 filters = 'ensembl_peptide_id',
                 values = union(ppi$protein1, ppi$protein2))

ppi.mapped <- na.omit(data.frame('gene.a' = mapping$ensembl_gene_id[match(ppi$protein1, table = mapping$ensembl_peptide_id)],
                                 'gene.b' = mapping$ensembl_gene_id[match(ppi$protein2, table = mapping$ensembl_peptide_id)],
                                 'score' = ppi$combined_score,
                                 stringsAsFactors = FALSE))
ppi.mapped <- rbind(ppi.mapped,
                    data.frame('gene.a' = ppi.mapped$gene.b,
                               'gene.b' = ppi.mapped$gene.a,
                               'score' = ppi.mapped$score,
                               stringsAsFactors = FALSE))
ppi.mapped <- ppi.mapped[order(ppi.mapped$score, decreasing = TRUE),]
idx <- duplicated(ppi.mapped[,1:2])
ppi.mapped <- ppi.mapped[!idx,]
rownames(ppi.mapped) <- 1:nrow(ppi.mapped)

all.ids <- sort(union(ppi.mapped$gene.a, ppi.mapped$gene.b))
row.idx <- match(ppi.mapped$gene.a, table = all.ids)
col.idx <- match(ppi.mapped$gene.b, table = all.ids)
vals <- ppi.mapped$score

net <- sparseMatrix(i = row.idx, j = col.idx, x = vals,
                    dims = c(length(all.ids), length(all.ids)),
                    dimnames = list(all.ids, all.ids),
                    use.last.ij = TRUE)
net <- Matrix::forceSymmetric(net)
saveRDS(net, file = 'data/stringdb/stringdb_net.rds')
