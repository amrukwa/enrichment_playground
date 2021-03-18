library(foreach)
library(doParallel)
source("source/gsea_polyaxon.R")

load("polyaxon.RData")
colnames(genesets) = c("ID", "Title", "features")

cores=detectCores()
cl <- makeCluster(cores[1])
clusterExport(cl, c("get_ES", "rank_genes", "single_gsea", "get_miss_increments", "get_miss_inc"))
registerDoParallel(cl)

# s2n_abs <- gsea(data, genesets, metaInfo, rank = 's2n')
# s2n <- gsea(data, genesets, metaInfo, rank = 's2n', absolute=FALSE)
# lfc_abs <- gsea(data, genesets, metaInfo, rank = 'lfc')
lfc <- gsea(data, genesets, metaInfo, rank = 'lfc', absolute=FALSE)

stopCluster(cl)

#save(lfc_abs, lfc, s2n_abs, s2n, file = "results.RData")
