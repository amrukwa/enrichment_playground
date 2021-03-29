library(plotly)

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

cerno_heatmaps <- function(dataset, genesets){
  gene_names <- rownames(dataset)
  patient_n <- length(dataset)
  plots <- list("CERNO_ABS", "CERNO_HTL", "CERNO_LTH")
  
  # get ranks for every patient
  x_normalized <- t(apply(dataset, 1, function(x) pnorm(x, mean=mean(x), sd=sd(x), lower.tail = TRUE)))
  abs_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="abs", gene_names=gene_names)
  htl_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="htl", gene_names=gene_names)
  lth_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="lth", gene_names=gene_names)
  # for each pathway calculate values
  for (i in 1:length(pathways)){
    pathway <- pathways[i]
    title <- pathway$MODULES$Title
    genes <- rownames(dataset) %in% pathway$GENES$ID
    CERNO_abs <- c()
    CERNO_htl <- c()
    CERNO_lth <- c()
    # parallelize it
    for (j in 1:patient_n){
      abs_res <- single_cerno(abs_labels[,j], pathway)
      htl_res <- single_cerno(htl_labels[,j], pathway)
      lth_res <- single_cerno(lth_labels[,j], pathway)
      CERNO_abs <- c(CERNO_abs, abs_res[, "F"])
      CERNO_htl <- c(CERNO_htl, htl_res[, "F"])
      CERNO_lth <- c(CERNO_lth, lth_res[, "F"])
      # save pvals as well
    }
  }
  # save the results as 6 matrixes <- pvals and Fs x6
  # x, y for heatmaps <- pathway names, patient numbers
}

cerno <- function(ordered_labels, genesets, combining_method="fisher"){
  results <- data.frame()
  for (i in 1:length(genesets)){
    result <- single_cerno(ordered_labels, genesets[i], combining_method)
    result <- cbind(ID=genesets[i]$MODULES$ID, Title=genesets[i]$MODULES$Title, result)
    results <- rbind(results, result)
  }
  results$corrected_pvals <- p.adjust(results$pval, method= "BH")
  results
}