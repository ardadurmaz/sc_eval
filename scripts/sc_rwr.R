library(Matrix)

net <- readRDS('data/stringdb/stringdb_net.rds')
sc_data <- readRDS('data/alevin_data/raw_res/deconv_res/alec/alevin_alectinib_deconv.rds')
expr_mat <- sc_data$norm_mat
rownames(expr_mat) <- gsub(pattern = '\\.\\d+$', replacement = '', rownames(expr_mat))

## Format ##
expr_mat <- as.matrix(expr_mat)
expr_mat <- expr_mat[match(rownames(net), table = rownames(expr_mat)),]
rownames(expr_mat) <- rownames(net)
expr_mat[is.na(expr_mat)] <- 0
expr_mat <- Matrix(expr_mat)


## Scale Expression ##
expr_mat <- t(scale(t(expr_mat)))
expr_mat <- 1/(1+exp(-1*expr_mat))
expr_mat[is.na(expr_mat)] <- 0

expr_mat <- expr_mat %*% Matrix::Diagonal(x = Matrix::colSums(expr_mat)^-1)
ppi.n <- net %*% Matrix::Diagonal(x = Matrix::colSums(net)^-1)

## RWR ##
r <- 0.5
delta <- 1
count <- 1

p0 <- expr_mat
pt <- Matrix(1/ncol(ppi.n), ncol = ncol(p0), nrow = nrow(p0))

while(delta > 1e-12 && count < 100){
  message('Processing: ', count, ' Delta: ', delta)
  px <- (1-r) * ppi.n %*% pt + r * p0
  delta <- max(Matrix::colSums(abs(px - pt)))
  pt <- px
  count <- count + 1
}
saveRDS(pt, file = 'data/alevin_alectinib_deconv_stringdb.rds')
