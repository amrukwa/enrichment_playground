load(".RData")
colnames(genesets) = c("ID", "Title", "features")
source("source/gsea_polyaxon.R")
ES <- single_gsea(data, genesets[1,], metaInfo, rank = 's2n')
ES

