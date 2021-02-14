# Packages
library(tmod)

# FOR UNIQUE GENE SETS
single_ora <- function(dataset, expressed_indices, geneset){
  tot <- unique(c(rownames(dataset), as.character(expressed_indices)))
  N  <- length(tot)
  K =  length(expressed_indices) # all significant genes
  M <- length(intersect(geneset$GENES$ID, tot)) # all genes in gene set and in dataset
  x <- length(intersect(as.numeric(geneset$GENES$ID), expressed_indices)) # significant genes in gene set
  # contingency table
  c_table = matrix(c(x, K-x, M-x, (N-M)-(K-x)), nrow=2, ncol=2)
  # Hypergeometric test
  pval = fisher.test(c_table, alternative="greater")$p.value
  pval
}

# FOR WHOLE TMOD OBJECT
ora <- function(dataset, expressed_indices, genesets){
  pvals <- c()
  for (i in 1:length(genesets)){
    pval <- single_ora(dataset, expressed_indices, genesets[i])
    pvals <- c(pvals, pval)
  }
  corrected_pvals <- p.adjust(pvals, method= "BH")
  for (i in 1:length(genesets)){
    if(corrected_pvals[i] < 0.05){
      print(genesets[i]$MODULES$Title)
    }
  }
  corrected_pvals
}