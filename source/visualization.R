library(umap)
library(plotly)
source("source/gsva.R")

group_comparison <- function(dataset, labels, pathway){
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
  
  
  title <- pathway$MODULES$Title
  genes <- rownames(dataset) %in% pathway$GENES$ID
  x <- dataset[genes,]
  x_umap <- umap(t(x))
  x_reduced <- x_umap$layout
  df <- as.data.frame(x_reduced)
  colnames(df) <- c("UMAP1", "UMAP2")
  
  # normalization of each row
  x_normalized <- t(apply(x, 1, function(x) pnorm(x, mean=mean(x), sd=sd(x), lower.tail = TRUE)))
  # integrations
  N <- length(genes)
  # from gsva
  df$GSVA <- t(gsva(dataset, pathway))
  # Z
  df$Z <- colSums(x)/sqrt(N)
  # SVD
  df$SVD <- t(decomposition(data, pathway))
  
  labels <- gsub("d", "disease", gsub("c", "control", labels$Group))
  df$Group <- labels
  
  fig1 <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, color= ~Group, type = 'scatter', mode = 'markers',
                  colors = "Dark2") %>% layout(annotations = c(annotation_style, text="Group"))
  fig2 <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', color= ~SVD, 
                  name="SVD", text = ~paste('Value: ', SVD),
                  mode = 'markers', showlegend=FALSE) %>% layout(annotations = c(annotation_style, text="SVD"))
  
  fig3 <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', color= ~GSVA, 
                  name="GSVA", text = ~paste('Value: ', GSVA),
                  mode = 'markers', showlegend=FALSE) %>% layout(annotations = c(annotation_style, text="GSVA"))
  
  fig4 <- plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', color= ~Z, 
                  name="Z", text = ~paste('Value: ', Z),
                  mode = 'markers', showlegend=FALSE) %>% layout(annotations = c(annotation_style, text="Z"))
  fig <- subplot(fig1, fig2, fig3, fig4, nrows = 2, shareX = TRUE, shareY = TRUE)
  fig <- fig %>% layout(autosize = F, width = 725, height = 500)
  fig
}