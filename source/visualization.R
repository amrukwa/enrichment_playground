library(umap)
library(plotly)
source("source/plage.R")
source("source/cerno.R")

annotation_style <- list(
  xref = "paper",
  yref = "paper",
  x = 0.5,
  y = 1.05,
  yanchor = "top",
  xanchor = "middle",
  align = "center",
  showarrow = FALSE,
  font = list(size = 15)
)

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

visualize_pathways <- function(dataset, labels, pathways, GSVA_ES){
  # prepare for visualizations
  x_umap <- umap(t(dataset))
  x_reduced <- x_umap$layout
  df <- as.data.frame(x_reduced)
  colnames(df) <- c("UMAP1", "UMAP2")
  # label each patient
  groups <- gsub("d", "disease", gsub("c", "control", labels$Group))
  df$Group <- groups
  patient_n <- length(dataset)
  plots1 <- list("Group", "GSVA", "Z", "SVD")
  plots2 <- list("Group", "CERNO_ABS", "CERNO_HTL", "CERNO_LTH")
  gradients <- list(plots1, plots2)
  
  # normalization of each row (Z scores for Z and CERNO)
  x_normalized <- t(apply(dataset, 1, function(x) pnorm(x, mean=mean(x), sd=sd(x), lower.tail = TRUE)))
  # prepare ranks for CERNO
  gene_names <- rownames(dataset)
  abs_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="abs", gene_names=gene_names)
  htl_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="htl", gene_names=gene_names)
  lth_labels <- apply(x_normalized, 2, rank_patient_genes, sort_type="lth", gene_names=gene_names)
  # prepare for subplots
  subplots <- vector(mode = "list", length = length(gradients))
  
  for (i in 1:length(pathways)){
    pathway <- pathways[i]
    title <- pathway$MODULES$Title
    genes <- rownames(dataset) %in% pathway$GENES$ID
    N <- length(genes)
    # choose ES from GSVA_ES
    df$GSVA <- t(GSVA_ES[title, ])
    # SVD
    df$SVD <- t(decomposition(dataset, pathway))
    # rows specific for the pathway
    x_norm <- x_normalized[genes,]
    # Z
    df$Z <- colSums(x_norm)/sqrt(N)
    # CERNO
    CERNO_abs <- c()
    CERNO_htl <- c()
    CERNO_lth <- c()
    for (j in 1:patient_n){
      abs_res <- single_cerno(abs_labels[,j], pathway)
      htl_res <- single_cerno(htl_labels[,j], pathway)
      lth_res <- single_cerno(lth_labels[,j], pathway)
      CERNO_abs <- c(CERNO_abs, abs_res[, "F"])
      CERNO_htl <- c(CERNO_htl, htl_res[, "F"])
      CERNO_lth <- c(CERNO_lth, lth_res[, "F"])
    }
    df$CERNO_ABS <- CERNO_abs
    df$CERNO_HTL <- CERNO_htl
    df$CERNO_LTH <- CERNO_lth
    fig <- group_comparison(df, gradients)
    subplots[[i]] <- fig
  }
  subplots
}

group_comparison <- function(df, gradients){
  pos <- 1
  plots <- vector(mode = "list", length = length(gradients))
  for (group in gradients){
    figs <- vector(mode = "list", length = length(group))
    for (i in 1:length(group)){
      gradient_name <- group[[i]]
      figs[[i]] <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', color= as.formula(paste0('~', gradient_name)), 
                           name=gradient_name, text = as.formula(paste0('~', gradient_name)),
                           mode = 'markers', showlegend=FALSE) %>% layout(annotations = c(annotation_style, text=gradient_name))
    }
    fig <- subplot(figs, nrows = 2, shareX = TRUE, shareY = TRUE)
    fig <- fig %>% layout(autosize = F, width = 725, height = 500)
    plots[[pos]] <- fig
    pos = pos+1
  }
  plots
}