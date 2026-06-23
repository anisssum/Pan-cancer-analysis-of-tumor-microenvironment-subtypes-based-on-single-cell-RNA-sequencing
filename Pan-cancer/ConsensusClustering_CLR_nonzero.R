#!/usr/bin/Rscript

library(ConsensusClustering)
library(compositions)
library(zCompositions)

data <- read.csv('/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/analysis_results_harmony2/cell_composition_matrix.csv', stringsAsFactors=FALSE, head = TRUE, row.names = 1, check.names = FALSE)
prop <- sweep(data, 2, colSums(data), "/")
prop <- prop + min(prop[prop>0])/2
clr_mat <- clr(t(prop))

Adj <- adj_mat(clr_mat,
               method="euclidean")

CM = consensus_matrix(Adj,
                      max.cluster = 6, 
                      resample.ratio = 0.7, 
                      max.itter = 500, 
                      clustering.method = 'pam')

Result <- cc_cluster_count(CM, plot.cdf = TRUE)
LogitScore = Result[["LogitScore"]]
PAC = Result[["PAC"]]
deltaA = Result[["deltaA"]]
CMavg = Result[["CMavg"]]

Kopt = Result[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))

save.image(file='consensus_matrix_CRC_ImmuneFilter_compact_clr_temp2.RData')
