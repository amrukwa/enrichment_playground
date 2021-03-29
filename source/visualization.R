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
  
  # normalization of each row (Z scores)
  x_normalized <- t(apply(dataset, 1, function(x) pnorm(x, mean=mean(x), sd=sd(x), lower.tail = TRUE)))
  # prepare for subplots
  subplots <- vector(mode = "list", length = length(pathways))
  
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

    fig <- group_comparison(df, plots1)
    subplots[[i]] <- fig
  }
  subplots
}

group_comparison <- function(df, gradients){
  figs <- vector(mode = "list", length = length(gradients))
  for (i in 1:length(gradients)){
    gradient_name <- gradients[[i]]
    figs[[i]] <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', color= as.formula(paste0('~', gradient_name)),
                         name=gradient_name, text = as.formula(paste0('~', gradient_name)),
                         mode = 'markers', showlegend=FALSE) %>% layout(annotations = c(annotation_style, text=gradient_name))
  }
  fig <- subplot(figs, nrows = 2, shareX = TRUE, shareY = TRUE)
  fig <- fig %>% layout(autosize = F, width = 725, height = 500)
  fig
}