library(matrixStats)

rank_genes <- function(data, labels, rank='s2n'){
  d_mat <- as.matrix(data)
  control_means <- rowMeans(d_mat[, labels$Group=="c"])
  disease_means <- rowMeans(d_mat[, labels$Group=="d"])
  if (rank == 's2n'){
    control_sds <- rowSds(d_mat[, labels$Group=="c"])
    disease_sds <- rowSds(d_mat[, labels$Group=="d"])
    ranks <- (disease_means - control_means)/(control_sds+disease_sds)
  }else if(rank=='lfc'){
    ranks <- log(disease_means/control_means)
  }
  ranks
}