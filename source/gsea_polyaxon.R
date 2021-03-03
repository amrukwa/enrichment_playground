source("source/ranks.R")
library(foreach)
library(doParallel)

get_ES <- function(geneset, data, labels){
  data = data[order(-data$ranks), ]
  is_hit <- row.names(data)  %in% unlist(strsplit(geneset["features"], ","))
  N_R <- sum(abs(data[is_hit, "ranks"]))
  P_hit <- rep(0, each=nrow(data))
  P_hit[is_hit] <- (abs(data[is_hit, "ranks"])/N_R)
  m_inc <- as.numeric(geneset["miss_increment"])
  P_miss <- rep(m_inc, each=nrow(data))
  P_miss[is_hit] <- 0
  ES <- cumsum(P_hit - P_miss)
  ES[which.max(abs(ES))]
}

get_miss_inc <- function(geneset, N){
  increment <- 1/(N- length(unlist(strsplit(geneset["features"], ","))))
  }

get_miss_increments <- function(genesets, N){
  increments <- apply(genesets, 1, get_miss_inc, N=N)
}

single_gsea <- function(data, genesets, labels, rank = "s2n", absolute=TRUE){
  data$ranks <- rank_genes(data, labels, rank)
  if (absolute){
    data$ranks <- abs(data$ranks)
  }
  ES <- apply(genesets, 1, get_ES, data=data, labels=labels)
  ES
}

gsea <- function(data, genesets, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
  genesets$miss_increment <- get_miss_increments(genesets, nrow(data))
  results <- data.frame(ID=genesets$ID, Title=genesets$Title) 
  results$ES <- single_gsea(data, genesets, labels, rank, absolute)
  pvals <- foreach(i=1:n_perm, .combine='+') %do% { 
    shuffled_labels <- transform( labels, Group = sample(Group) )
    es <- single_gsea(data, genesets, shuffled_labels, rank, absolute)
    better <- rep(0, nrow(genesets))
    for (i in 1:nrow(genesets)){
      if(results[i, "ES"] > 0 & results[i, "ES"] < es[i]){
        better[i] <- 1
      }else if(results[i, "ES"] < 0 & results[i, "ES"] > es[i]){
        better[i] <- 1
      }
    }
    better
  }
  results$pval <- pvals/n_perm
  results$corrected_pval <- p.adjust(results$pval, method="BH")
  results
}

# for single permutation, takes already shuffled labels
# single_gsea_gs <- function(data, geneset, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
#   N <- nrow(data)
#   N_H <- length(unlist(strsplit(geneset[, 3], ",")))
#   print(N_H)
#   miss_increment <- 1/(N-N_H)
#   data$ranks <- rank_genes(data, labels, rank)
#   if (absolute){
#     data$ranks <- abs(data$ranks)
#   }
#   ES <- get_ES(data, geneset, labels, miss_increment, rank, absolute)
#   if (ES > 0){
#     better <- foreach (i=1:n_perm, .combine='+', .export= c("get_ES", "rank_genes"), .packages=c("matrixStats")) %dopar% {
#       shuffled_labels <- transform( labels, Group = sample(Group) )
#       data$ranks <- rank_genes(data, labels, rank)
#       if (absolute){
#         data$ranks <- abs(data$ranks)
#       }
#       es <- get_ES(data, geneset, miss_increment, rank)
#       better_one <- 0
#       if (es > ES){
#         better_one <- 1}
#       better_one
#     }}else{
#       better <- foreach (i=1:n_perm, .combine='+', .export= c("get_ES", "rank_genes"), .packages=c("matrixStats")) %dopar% {
#         shuffled_labels <- transform( labels, Group = sample(Group) )
#         es <- get_ES(data, geneset, shuffled_labels, miss_increment, rank, absolute)
#         better_one <- 0
#         if (es < ES){
#           better_one <- 1}
#         better_one}}
#   p_val <- better/n_perm
#   data.frame(ES=ES, p_val = p_val)
# }
# 
# # pvals <- foreach(i=1:n_perm, .combine='+') %do% { shuffled_labels <- transform( labels, Group = sample(Group) )
# # c(1:5) gene set stuff} get p-values for all permutations at once
# 
# gsea_gs <- function(data, genesets, labels, rank="s2n", absolute=TRUE, n_perm=1000){
#   results <- data.frame()
#   cores=detectCores()
#   cl <- makeCluster(cores[1]-1)
#   clusterExport(cl, c("get_ES", "rank_genes"))
#   registerDoParallel(cl)
#   for (i in 1:nrow(genesets)){
#     result <- single_gsea_gs(data, genesets[i, ], labels, rank, absolute, n_perm)
#     results <- rbind(results, result)
#   }
#   stopCluster(cl)
#   pvals_adjusted = p.adjust(results$p_val, method="BH")
#   x = data.frame(ID=genesets$ID, Title=genesets$Title, corrected_pval=pvals_adjusted)
#   results = cbind(results, x)
# }