source("source/ranks.R")
library(foreach)
library(doParallel)

get_ES <- function(data, geneset, labels, miss_increment, rank = "s2n", absolute=TRUE){
  data$ranks <- rank_genes(data, labels, rank)
  if (absolute){
    data$ranks <- abs(data$ranks)
  }
  data = data[order(-data$ranks), ]
  is_hit <- row.names(data)  %in% geneset$GENES$ID
  data$ranks <- abs(data$ranks)
  N_R <- sum(data[is_hit, "ranks"])
  data[is_hit, "ranks"] <- data[is_hit, "ranks"]/N_R
  es <- if(is_hit[1]) data[1, "ranks"] else (-miss_increment)
  ES <- es
  for (i in 2:nrow(data)){
    es <- if(is_hit[i]) (es+data[i, "ranks"]) else (es-miss_increment)
    if(abs(es) > abs(ES)){
      ES <- es
    }
  }
  ES
}


single_gsea <- function(data, geneset, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
  N <- nrow(data)
  N_H <- length(geneset$GENES$ID)
  miss_increment <- 1/(N-N_H)
  
  ES <- get_ES(data, geneset, labels, miss_increment, rank, absolute)
  if (ES > 0){
    better <- foreach (i=1:n_perm, .combine='+', .export= c("get_ES", "rank_genes"), .packages=c("matrixStats")) %dopar% {
      shuffled_labels <- transform( labels, Group = sample(Group) )
      es <- get_ES(data, geneset, shuffled_labels, miss_increment, rank, absolute)
      better_one <- 0
      if (es > ES){
        better_one <- 1}
      better_one
    }}else{
      better <- foreach (i=1:n_perm, .combine='+', .export= c("get_ES", "rank_genes"), .packages=c("matrixStats")) %dopar% {
        shuffled_labels <- transform( labels, Group = sample(Group) )
        es <- get_ES(data, geneset, shuffled_labels, miss_increment, rank, absolute)
        better_one <- 0
        if (es < ES){
          better_one <- 1}
        better_one}}
  p_val <- better/n_perm
  data.frame(ES=ES, p_val = p_val)
}

gsea <- function(data, genesets, labels, rank="s2n", absolute=TRUE, n_perm=1000){
  results <- data.frame()
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  clusterExport(cl, c("get_ES", "rank_genes"))
  registerDoParallel(cl)
  for (i in 1:length(genesets)){
    result <- single_gsea(data, genesets[i, ], labels, rank, absolute, n_perm)
    results <- rbind(results, result)
  }
  stopCluster(cl)
  pvals_adjusted = p.adjust(results$p_val, method="BH")
  x = data.frame(ID=genesets$MODULES$ID, Title=genesets$MODULES$Title, corrected_pval=pvals_adjusted)
  results = cbind(results, x)
}
