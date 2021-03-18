library(tmod)

decomposition <- function(dataset, geneset){
  genes <- rownames(dataset) %in% geneset$GENES$ID
  x <- t(dataset[genes,])
  values <- prcomp(x, center=TRUE, scale.=TRUE, rank.=1)$x
  t(values)
}

plage <- function(dataset, genesets, labels){
  
}