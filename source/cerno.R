# as in tmod the function takes already sorted gene list, my implementation will do the same

single_cerno <- function(ordered_labels, geneset, combining_method="fisher"){
  N <- length(ordered_labels) # all analyzed genes
  all_ranks <- c(1:N)
  ranked_labels <- data.frame(id=ordered_labels, rank = all_ranks) # all genes and their ranks
  # choose by ids that are in GS
  in_gs <- ranked_labels[ranked_labels$id %in% geneset$GENES$ID,]
  N_gs <- nrow(in_gs)
  sum_of_ranks <- sum(in_gs$rank)
  AUC <- (N_gs*(N-N_gs) + N_gs*(N_gs+1)/2 - sum_of_ranks)/(N_gs*(N-N_gs))
  # the values to be combined
  in_gs$P <- in_gs$rank/N
  if(combining_method=="fisher"){
    stat_val <- (-2)*sum(log(in_gs$P))
    pval <- pchisq(stat_val, 2*N_gs, lower.tail = FALSE)
    cES <- stat_val/(N_gs*2)
    return (data.frame(F=stat_val, N_gs = N_gs, AUC=AUC, cES=cES, pval=pval))
  }
  # z transform version
  if(combining_method=="z-transform"){
    stat_val <- sum(qnorm(in_gs$P))/sqrt(N_gs)
    pval <- pnorm(stat_val)
    return (data.frame(z=stat_val, N_gs = N_gs, AUC=AUC, pval=pval))
  }
}

cerno <- function(ordered_labels, genesets, combining_method="fisher"){
  results <- data.frame()
  for (i in 1:length(genesets)){
    result <- single_cerno(ordered_labels, genesets[i], combining_method)
    result <- cbind(ID=genesets[i]$MODULES$ID, Title=genesets[i]$MODULES$Title, result)
    results <- rbind(results, result)
  }
  results$corrected_pvals <- p.adjust(results$pval, method= "BH")
  results
}