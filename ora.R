# Packages
library(tmod)

# FOR UNIQUE GENE SETS
single_ora <- function(N, expressed_indices, geneset){
  # N - number of genes in dataset
  K =  length(expressed_indices) # all significant genes
  M <- length(geneset$GENES$ID) # all genes in gene set
  x <- length(intersect(as.numeric(geneset$GENES$ID), expressed_indices)) # significant genes in gene set
  # contingency table
  c_table = matrix(c(x, K-x, M-x, (N-M)-(K-x)), nrow=2, ncol=2)
  # Hypergeometric test
  pval = fisher.test(c_table, alternative="greater")$p.value
  pval
}

# FOR WHOLE TMOD OBJECT
