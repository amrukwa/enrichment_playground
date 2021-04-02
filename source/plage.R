library(tmod)

decomposition <- function(dataset, geneset){
  genes <- rownames(dataset) %in% geneset$GENES$ID
  x <- t(dataset[genes,])
  values <- prcomp(x, center=TRUE, scale.=TRUE, rank.=1)$x
  t(values)
}

plage <- function(dataset, genesets, labels){
  pvals <- c()
  stat <- c()
  for (i in 1:length(genesets$MODULES$Title)){
    dec <- decomposition(dataset, genesets[i])
    res <- t.test(dec[labels$Group=="d"], dec[labels$Group=="c"])
    pvals <- append(pvals, res$p.value)
    stat <- append(stat, res$statistic)
  }
  corrected_pvals <- p.adjust(pvals, method= "BH")
  data.frame(ID=genesets$MODULES$ID, Title=genesets$MODULES$Title, t=stat, pval=pvals, corrected_pvals=corrected_pvals)
}