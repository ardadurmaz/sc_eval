library(Matrix)


suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))

run_title <- 'Alevin_Alectinib'
setwd(sprintf('~/sc_eval/data/%s', run_title))

proc_data <- readRDS(sprintf('%s_ProcessedData.rds', run_title))

## Trajectory Plots ##
for(local_n in names(traject(proc_data))){
  local_res <- traject(proc_data, local_n)
  write.table(local_res, 
              file = 'TrajectoryPlotData.txt', 
              sep = '\t', 
              col.names = TRUE, 
              row.names = FALSE, 
              append = FALSE, 
              quote = FALSE)
  if(grepl(pattern = 'slingshot', local_n)){
    system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                   'TrajectoryPlotData.txt',
                   sprintf('%s_%s.pdf slingshot', run_title, local_n)))
  }else{
    system(sprintf('python3 ~/sc_eval/scripts/plot_3d.py %s %s',
                   'TrajectoryPlotData.txt',
                   sprintf('%s_%s.pdf monocle', run_title, local_n)))
  }
}
quit(save='no')


## Temp ##
library(scran)
sce <- SingleCellExperiment(list(counts = counts(proc_data, 'raw')))
local_clusters <- quickCluster(sce, 
                               min.size = 20, 
                               method = 'igraph', 
                               d = NA, 
                               use.ranks = TRUE, 
                               min.mean = 1)
sce <- computeSumFactors(sce, clusters = local_clusters)

## Variance Trend
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.3), use.spikes=FALSE)
decomp <- decomposeVar(sce, var.fit)
top.hvgs <- order(decomp$bio, decreasing=TRUE)
head(decomp[top.hvgs,])

sce <- normalize(sce, return_log = FALSE)
norm_counts <- normcounts(sce)
means <- rowMeans(norm_counts)
cv2 <- apply(norm_counts, 1, var)/(means^2)
dm.stat <- DM(means, cv2)


norm_mat <- as.matrix(logcounts(sce))

norm_expr <- norm(proc_data, 'raw_deconv')
cv_res <- improvedCV2(norm_expr,)
