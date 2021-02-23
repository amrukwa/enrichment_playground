source("source/ranks.R")
library(foreach)
library(doParallel)

ES_i <- function(hit, norm_rank, miss_increment, prev_es){
  if(hit){
    es <- prev_es + norm_rank
  }else{
    es <- prev_es - miss_increment
  }
  es
}


get_ES <- function(data, geneset, labels, miss_increment, is_hit, rank = "s2n", absolute=TRUE){
  data$ranks <- rank_genes(data, labels, rank)
  if (absolute){
    data$ranks <- abs(data$ranks)
  }
  data = data[order(-data$ranks), ]
  N_R <- sum(data[is_hit, "ranks"])
  es <- ES_i(is_hit[1], (data[1, "ranks"]/N_R), miss_increment, 0)
  ES <- es
  for (i in 2:nrow(data)){
    es <- ES_i(is_hit[i], (data[i, "ranks"]/N_R), miss_increment, es)
    if(abs(es) > abs(ES)){
      ES <- es
    }
  }
  ES
}


gsea <- function(data, geneset, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  clusterExport(cl, c("get_ES", "rank_genes", "ES_i"))
  registerDoParallel(cl)
  
  N <- nrow(data)
  N_H <- length(geneset$GENES$ID)
  miss_increment <- 1/(N-N_H)
  is_hit <- row.names(data)  %in% geneset$GENES$ID
  
  ES <- get_ES(data, geneset, labels, miss_increment, is_hit, rank, absolute)
  better <- foreach (i=1:n_perm, .combine='+', .export= c("get_ES", "rank_genes", "ES_i"), .packages=c("matrixStats")) %dopar% {
    shuffled_labels <- transform( labels, Group = sample(Group) )
    es <- get_ES(data, geneset, shuffled_labels, miss_increment, is_hit, rank, absolute)
    better_one <- 0
    if (es > ES){
      # the stuff taking the sign into account, now only for absolute bigger
      better_one <- 1
    }
    better_one
  }
  stopCluster(cl)
  p_val <- better/n_perm
}