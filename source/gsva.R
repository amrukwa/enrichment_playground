library(matrixStats)
library(foreach)
library(doParallel)

transform_expression <- function(gene_expression, gene_bw){
  gene_expression <- gene_expression/gene_bw
  expressions <- c()
  for (patient in gene_expression){
    transformed_values <- (gene_expression)*(-1) + patient
    transformed_expression <- mean(pnorm(transformed_values))
    expressions <- c(expressions, transformed_expression)
  }
  expressions
}

get_m_inc <- function(geneset, genes){
  N = length(genes)
  in_gs <- which(genes %in% geneset$GENES)
  increment <- 1/(N- length(in_gs))
}

get_m_incs <- function(genesets, genes){
  N <- length(genesets)
  increments <- rep(0, N)
  for (i in 1:N){
    increments[i] <- get_m_inc(genesets[i], genes)
  }
  increments
}

calculate_ES <- function(ranks, is_hit, miss_increment){
  N_R <- sum(ranks[is_hit])
  ranks[is_hit] = ranks[is_hit]/N_R
  es <- ES_max <- ES_min <- 0
  for (i in 1:length(ranks)){
    es <- if(is_hit[i]) (es+ranks[i]) else (es-miss_increment)
    if (es < ES_min){
      ES_min = es
    }else if(es > ES_max){
      ES_max = es}
  }
  ES <- ES_max - ES_min
}

gsva <- function(dataset, genesets){
  cores=detectCores()-1
  cl <- makeCluster(cores[1])
  clusterExport(cl, c("calculate_ES", "get_m_incs", "transform_expression", "get_m_inc"))
  registerDoParallel(cl)
  
  N <- nrow(dataset)
  genes = row.names(dataset)
  # for each gene in the dataset calculate bandwidth
  d_mat <- data.matrix(dataset)
  bandwidths <- rowSds(d_mat)/4
  
  # transform gene expressions
  transformed <- foreach(i = 1:N, .combine=rbind) %do% {
    transformed_gene <- transform_expression(d_mat[i, ], bandwidths[i])
    transformed_gene
  }
  
  # normalize the results
  transformed <- abs(transformed - N/2)
  # miss increments for all genesets
  m_inc <- get_m_incs(genesets, genes)
  
  ES <- foreach(i = 1:ncol(dataset), .combine=cbind) %do% {
    # rank the genes for this patient
    ranks <- transformed[,i]
    ranks_order <- order(ranks, decreasing=TRUE)
    ranks <- ranks[ranks_order]
    gene_names <- row.names(dataset)[ranks_order]
    
    # calculate ES for each geneset for this patient
    foreach(j = 1:length(genesets), .combine=c) %do% {
      is_hit <- gene_names  %in% genesets[j]$GENES$ID
      es <- calculate_ES(ranks, is_hit, m_inc[j])
      es
      }
    }
  stopCluster(cl)
  res <- data.frame(ES)
  colnames(res) <- NULL
  rownames(res) <- genesets$MODULES$ID
  res
}
