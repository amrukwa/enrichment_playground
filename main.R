library(tmod)

# Load the data

load("data/data_lung_cancer.RData")


# AUCell
library(AUCell)
cells_rankings <- AUCell_buildRankings(data.matrix(data))
geneSets <- list()
for (i in 1:length(KEGGhsa)){
  title = KEGGhsa[i]$MODULES$Title
  genes <- list(KEGGhsa[i]$GENES$ID)
  geneSets <- c(geneSets, genes)
}
names(geneSets) <- KEGGhsa$MODULES$Title
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

set.seed(123)
miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]
library(NMF)
pdf("heatmap.pdf")
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)
dev.off()

set.seed(NULL)
library(umap)
x_umap <- umap(t(data))
x_reduced <- x_umap$layout
df <- as.data.frame(x_reduced)
colnames(df) <- c("UMAP1", "UMAP2")
labels = gsub("d", "red", gsub("c", "blue", metaInfo$Group))
plot(df, col=labels, bg=labels, pch=16)
legend("topright", legend=c("disease", "control"),
       col=c("red", "blue"), pt.cex = 1, pch=16)

set.seed(123)
selectedThresholds <- getThresholdSelected(cells_assignment)
par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
for(geneSetName in names(selectedThresholds))
{
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    plot(df, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(df)], pch=16)
  }
}
