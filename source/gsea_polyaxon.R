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
  es <- foreach(i=1:n_perm, .combine=rbind, .packages = "matrixStats") %dopar% { 
    shuffled_labels <- transform(labels, Group = sample(Group))
    es <- single_gsea(data, gs, shuffled_labels, rank, absolute)
    es
  }
  ES_is_pos <- (ES >= 0)
  
  NES <- foreach(i=1:length(ES), .combine=c) %dopar% {
    if (ES_is_pos[i]){
      es_ <- es[es[, i] > 0, i]
    }else{
      es_ <- es[es[, i] < 0, i]
    }
    NES = ES[i]/abs(mean(es_))
    NES
  }
  
  nes <- foreach(i=1:length(ES), .combine=cbind) %dopar% {
    nes <- rep(0, n_perm)
    if (ES_is_pos[i]){
      sign <- es[, i] > 0
    }else{
      sign <- es[, i] < 0
    }
    es_ <- es[sign, i]
    nes[sign] <- es_/abs(mean(es_))
    nes
  }

  nes_pos <- nes[nes>0]
  nes_neg <- nes[nes<0]
  
  res = foreach(i=1:length(ES), .combine=rbind)%dopar%{
    if (ES_is_pos[i]){
      better <- (nes_pos > NES[i])
      pval <- length(which(es[, i] >= ES[i]))/length(which(es[ , i] > 0))
      
      NES_pval <- max(sum(better)/length(nes_pos), 1/length(nes_pos))
      NES_qval <- min(1, NES_pval/((sum(NES >= NES[i]))/sum(ES_is_pos)))
      }else{
      better <- (nes_neg < NES[i])
      pval <- length(which(es[, i] <= ES[i]))/length(which(es[ , i] < 0))
      
      NES_pval <- max(sum(better)/length(nes_neg), 1/length(nes_neg))
      NES_qval <- min(1, NES_pval/((sum(NES <= NES[i]))/(sum(!ES_is_pos))))
      }
    data.frame(pval=pval,NES_pval=NES_pval, NES_qval=NES_qval)
  }
  results$ES <- ES
  results$NES <- NES
  results <- cbind(results, res)
  results$pval_corrected <- p.adjust(results$pval, method="BH")
  results
  results[,c('ID', 'Title','ES','pval','pval_corrected', 'NES', 'NES_pval', 'NES_qval')]
}
