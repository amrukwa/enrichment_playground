source("source/ranks.R")

get_ES <- function(geneset, ranks, genes){
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
  ES <- apply(genesets, 1, get_ES, ranks=ranks, genes=genes)
  ES
}

gsea <- function(data, gs, labels, rank = "s2n", absolute=TRUE, n_perm=1000){
  genes <- row.names(data)
  gs$miss_increment <- get_miss_increments(gs, genes)
  results <- data.frame(ID=gs$ID, Title=gs$Title) 
  ES <- single_gsea(data, gs, labels, rank, absolute)
  es <- foreach(i=1:n_perm, .combine=rbind) %dopar% { 
    shuffled_labels <- transform(labels, Group = sample(Group))
    es <- single_gsea(data, gs, shuffled_labels, rank, absolute)
    es
  }
  NES <- foreach(i=1:length(ES), .combine=c) %dopar% {
    NES = ES[i]/abs(mean(es[, i]))
    NES
  }
  res = foreach(i=1:length(ES),.combine=rbind)%dopar%{
    if (ES[i] >= 0){
      better <- (es[, i] > ES[i])
      pval <- length(which(es[, i] >= ES[i]))
      NES_pval <- max(sum(es[better, i])/n_perm, 1/n_perm)
      NES_qval <- min(1, NES_pval/(sum(NES[better]>NES[i])/sum(better)))
      }else{
      better <- (es[, i] < ES[i])
      pval <- length(which(es[, i] <= ES[i]))
      NES_pval <- max(sum(es[better, i])/n_perm, 1/n_perm)
      NES_qval <- min(1, NES_pval[i]/(sum(NES[better] < NES[i])/sum(better)))
      }
    data.frame(pval=pval,NES_pval=NES_pval, NES_qval=NES_qval)
  }
  results$NES <- NES
  results <- cbind(results, res)
  results
}
