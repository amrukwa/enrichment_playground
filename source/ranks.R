library(matrixStats)

rank_genes <- function(data, labels, rank='s2n'){
  d_mat <- as.matrix(data)
  control <- d_mat[, labels$Group=="c"]
  disease <- d_mat[, labels$Group=="d"]
  control_means <- rowMeans(control)
  disease_means <- rowMeans(disease)
  if (rank == 's2n'){
    control_sds <- rowSds(control)
    disease_sds <- rowSds(disease)
    ranks <- (disease_means - control_means)/(control_sds+disease_sds)
  }else if(rank=='lfc'){
    ranks <- log(disease_means/control_means)
  }
  ranks
}