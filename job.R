load(".RData")
colnames(genesets) = c("ID", "Title", "features")
source("source/gsea_polyaxon.R")
ES <- gsea(data, genesets, metaInfo, rank = 's2n')
ES

