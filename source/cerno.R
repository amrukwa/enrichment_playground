library(plotly)
library(foreach)
library(doParallel)
library(heatmaply)

# as in tmod the function takes already sorted gene list, my implementation will do the same

rank_patient_genes <- function(gene_expressions, sort_type="abs", gene_names){
  if(sort_type=="abs"){
    ordered <- gene_names[order(abs(gene_expressions), decreasing = TRUE)]
  }else if(sort_type=="htl"){
    ordered <- gene_names[order(gene_expressions, decreasing = TRUE)]
  }else{
    ordered <- gene_names[order(gene_expressions)]
  }
  ordered
}

single_patient_cerno <- function(ordered_labels, pathway, N, all_ranks){
  # choose by ids that are in GS
  in_gs <- ordered_labels %in% pathway$GENES$ID
  significant_ranks <- all_ranks[in_gs]
  sum_of_ranks <- sum(significant_ranks)
  N_gs <- length(significant_ranks)
  # the values to be combined
  P <- significant_ranks/N
  stat_val <- (-2)*sum(log(P))
  pval <- pchisq(stat_val, 2*N_gs, lower.tail = FALSE)
  pval
}

patients_pathway_cerno <- function(ordered_labels, pathway, N, all_ranks, patients_n) {
  results <- foreach (j = 1:patients_n, .combine=rbind) %dopar% {
    res <- single_patient_cerno(ordered_labels[,j], pathway, N, all_ranks)
  }
  results
}

cerno_heatmaps <- function(dataset, pathways, color_labels=NULL, sort_type="abs", with_dendro='both'){
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  clusterExport(cl, c("rank_patient_genes", "patients_pathway_cerno", "single_patient_cerno"))
  registerDoParallel(cl)
  pathways_number <- length(pathways$MODULES$Title)
  gene_names <- rownames(dataset)
  patient_n <- length(dataset)
  N <- nrow(dataset)
  all_ranks <- c(1:N)
  
  # get ranks for every patient
  x_normalized <- t(apply(dataset, 1, function(x) pnorm(x, mean=mean(x[color_labels=='c']), sd=sd(x[color_labels=='c']), lower.tail = TRUE)))
  labels <- apply(x_normalized, 2, rank_patient_genes, sort_type=sort_type, gene_names=gene_names)
  # for each pathway calculate values
  pvalues <- matrix(NA, nrow=pathways_number, ncol=patient_n)
  for (i in 1:pathways_number){
    pathway <- pathways[i]
    values <- patients_pathway_cerno(labels, pathway, N, all_ranks, patient_n)
    pvalues[i,] <- values
  }
  stopCluster(cl)
  rownames(pvalues) <- pathways$MODULES$Title
  colnames(pvalues) <- colnames(dataset)
  ByPal <- colorRampPalette(c('blue', 'red'))
  fig <- heatmaply(pvalues, k_row = 2, k_col = 2, dendrogram = with_dendro,
                   col_side_colors = data.frame("Group" = color_labels, check.names=FALSE),
                   col_side_palette = ByPal,   showticklabels = c(FALSE, TRUE),
                   width=1000, height=625, fontsize_col=8, fontsize_row=8)
  fig
}

single_cerno <- function(ordered_labels, geneset, combining_method="fisher"){
  N <- length(ordered_labels) # all analyzed genes
  all_ranks <- c(1:N)
  ranked_labels <- data.frame(id=ordered_labels, rank = all_ranks) # all genes and their ranks
  # choose by ids that are in GS
  in_gs <- ranked_labels[ranked_labels$id %in% geneset$GENES$ID,]
  N_gs <- nrow(in_gs)
  sum_of_ranks <- sum(in_gs$rank)
  AUC <- (N_gs*(N-N_gs) + N_gs*(N_gs+1)/2 - sum_of_ranks)/(N_gs*(N-N_gs))
  # the values to be combined
  in_gs$P <- in_gs$rank/N
  if(combining_method=="fisher"){
    stat_val <- (-2)*sum(log(in_gs$P))
    pval <- pchisq(stat_val, 2*N_gs, lower.tail = FALSE)
    cES <- stat_val/(N_gs*2)
    return (data.frame(F=stat_val, N_gs = N_gs, AUC=AUC, cES=cES, pval=pval))
  }
  # z transform version
  if(combining_method=="z-transform"){
    stat_val <- sum(qnorm(in_gs$P))/sqrt(N_gs)
    pval <- pnorm(stat_val)
    return (data.frame(z=stat_val, N_gs = N_gs, AUC=AUC, pval=pval))
  }
}

cerno <- function(ordered_labels, genesets, combining_method="fisher"){
  results <- data.frame()
  for (i in 1:length(genesets$MODULES$Title)){
    result <- single_cerno(ordered_labels, genesets[i], combining_method)
    result <- cbind(ID=genesets[i]$MODULES$ID, Title=genesets[i]$MODULES$Title, result)
    results <- rbind(results, result)
  }
  results$corrected_pvals <- p.adjust(results$pval, method= "BH")
  results
}