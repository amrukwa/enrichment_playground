library(matrixStats)
library(foreach)
library(doParallel)
source("source/expression_testing.R")

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
  miss_increment <- miss_increment*(-1)
  N_R <- sum(ranks[is_hit])
  ranks[is_hit] = ranks[is_hit]/N_R
  ranks[!is_hit] = miss_increment
  es <- ES_max <- ES_min <- 0
  for (i in 1:length(ranks)){
    es <- es + ranks[i]
    if (es < ES_min){
      ES_min = es
    }else if(es > ES_max){
      ES_max = es}
  }
  ES <- ES_max - abs(ES_min)
}

gsva <- function(dataset, genesets, labels, rownames_title=TRUE){
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  clusterExport(cl, c("calculate_ES", "get_m_incs", "transform_expression", "get_m_inc"))
  registerDoParallel(cl)
  N <- nrow(dataset)
  genes = row.names(dataset)
  # for each gene in the dataset calculate bandwidth
  d_mat <- data.matrix(dataset)
  bandwidths <- rowSds(d_mat)/4
  
  # transform gene expressions
  transformed <- foreach(i = 1:N, .combine=rbind) %dopar% {
    transformed_gene <- transform_expression(d_mat[i, ], bandwidths[i])
    transformed_gene
  }
  # miss increments for all genesets
  m_inc <- get_m_incs(genesets, genes)
  ES <- foreach(i = 1:ncol(dataset), .combine=cbind) %do% {
    # rank the genes for this patient
    ranks <- transformed[,i]
    ranks_order <- order(ranks, decreasing=TRUE)
    ranks <- ranks[ranks_order]
    gene_names <- row.names(dataset)[ranks_order]
    # normalize the ordered genes
    ranks <- abs(seq(from=N, to=1) - N/2)
    # calculate ES for each geneset for this patient
    patient_es <- c()
    for (j in 1:length(genesets)){
      is_hit <- gene_names  %in% genesets[j]$GENES$ID
      es <- calculate_ES(ranks, is_hit, m_inc[j])
      patient_es <- c(patient_es, es)
      }
    patient_es
  }
  stopCluster(cl)
  ES <- data.frame(ES)
  colnames(ES) <- NULL
  if(rownames_title){
    rownames(ES) <- genesets$MODULES$Title}else{
    rownames(ES) <- genesets$MODULES$ID
    }
  ES$pval <- apply(ES, 1, do_ftest, labels)
  ES$pval <- p.adjust(ES$pval, method= "BH")
  ES$pval <- apply(ES, 1, do_ttest, labels, colname="pval")
  ES$corrected_pval <- p.adjust(ES$pval)
  pvals <- ES[, (ncol(ES)-1):ncol(ES)]
  ES <- ES[1:(length(ES)-2)]
  list(ES = ES, pvals= pvals)
}
