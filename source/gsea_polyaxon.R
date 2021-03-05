source("source/ranks.R")

get_ES <- function(geneset, ranks, genes, labels){
  is_hit <- genes  %in% unlist(strsplit(geneset["features"], ","))
  N_R <- sum(ranks[is_hit]) 
  ranks[is_hit] = ranks[is_hit]/N_R
  es <- 0
  ES <- es
  m_inc <- as.numeric(geneset["miss_increment"])
  for (i in 1:length(ranks)){
    es <- if(is_hit[i]) (es+ranks[i]) else (es-m_inc)
    if(abs(es) > abs(ES)){
      ES <- es
    }
  }
  ES
}

get_miss_inc <- function(geneset, genes){
  N = length(genes)
  in_gs <- which(genes %in% unlist(strsplit(geneset["features"], ",")))
  increment <- 1/(N- length(in_gs))
  }

get_miss_increments <- function(genesets, genes){
  increments <- apply(genesets, 1, get_miss_inc, genes=genes)
}

single_gsea <- function(data, genesets, labels, rank = "s2n", absolute=TRUE){
  data$ranks <- rank_genes(data, labels, rank)
  if (absolute){
    data$ranks <- abs(data$ranks)
  }
  data = data[order(-data$ranks), ]
  genes = row.names(data)
  ranks = abs(data$ranks)
  ES <- apply(genesets, 1, get_ES, ranks=ranks, genes=genes, labels=labels)
  ES
}

gsea <- function(data, gs, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
  genes <- row.names(data)
  gs$miss_increment <- get_miss_increments(gs, genes)
  results <- data.frame(ID=gs$ID, Title=gs$Title) 
  ES <- single_gsea(data, gs, labels, rank, absolute)
  
  pvals <- foreach(i=1:n_perm, .combine='+') %do% { 
    shuffled_labels <- transform( labels, Group = sample(Group) )
    es <- single_gsea(data, gs, shuffled_labels, rank, absolute)
    better <- rep(0, nrow(gs))
    for (i in 1:nrow(gs)){
      if(ES[i] > 0 & ES[i] < es[i]){
        better[i] <- 1
      }else if(ES[i] < 0 & ES[i] > es[i]){
        better[i] <- 1
      }
    }
    better
  }
  results$ES <- ES
  results$pval <- pvals/n_perm
  #results$corrected_pval <- p.adjust(results$pval, method="BH")
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