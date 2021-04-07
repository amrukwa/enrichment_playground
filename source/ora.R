# Packages
library(tmod)

# FOR UNIQUE GENE SETS
single_ora <- function(dataset, expressed_indices, geneset){
  N <- nrow(dataset)
  K <-  length(expressed_indices) # all significant genes
  M <- length(intersect(geneset$GENES$ID, rownames(dataset))) # all genes in gene set
  x <- length(intersect(as.numeric(geneset$GENES$ID), expressed_indices)) # significant genes in gene set
  # contingency table
  c_table = matrix(c(x, K-x, M-x, (N-M)-(K-x)), nrow=2, ncol=2)
  or <- (c_table[1,1]*c_table[2,2])/(c_table[2,1]*c_table[1,2])
  # Hypergeometric test
  pval = fisher.test(c_table, alternative="greater")$p.value
  data.frame(x=x, M=M, K=K, N=N, odds_ratio=or, pval=pval)
}

# FOR WHOLE TMOD OBJECT
ora <- function(dataset, expressed_indices, genesets){
  results <- data.frame()
  for (i in 1:length(genesets$MODULES$Title)){
    result <- single_ora(dataset, expressed_indices, genesets[i])
    result <- cbind(ID=genesets[i]$MODULES$ID, Title=genesets[i]$MODULES$Title, result)
    results <- rbind(results, result)
  }
  results$corrected_pvals <- p.adjust(results$pval, method= "BH")
  results
}