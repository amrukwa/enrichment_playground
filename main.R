# Load source files
source("source/expression_testing.R")
source("source/ora.R")
source("source/cerno.R")

library(tmod)

# Load the data

load("data_lung_cancer.RData")

heatmap_all_cerno <- cerno_heatmaps(data, KEGGhsa, color_labels=metaInfo$Group, sort_type="abs", with_dendro='none')

dendro_all_cerno <- cerno_heatmaps(data, KEGGhsa, color_labels=metaInfo$Group, sort_type="abs", with_dendro='both')

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
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)

library(umap)
x_umap <- umap(t(data))
x_reduced <- x_umap$layout
df <- as.data.frame(x_reduced)
colnames(df) <- c("UMAP1", "UMAP2")
labels = gsub("d", "red", gsub("c", "blue", metaInfo$Group))
plot(df)

plot(df, col=labels)
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

# My GSVA
source("source/gsva.R")
gsva_implementation <- gsva
gsva_result_implementation <- gsva_implementation(data, KEGGhsa[1], metaInfo$Group)


# TMOD PLAGE
plage_result_tmod <- tmodPLAGEtest(mset=KEGGhsa, x=data, group=metaInfo$Group, l=rownames(data), order.by="none")
plage_result_tmod

# My implementation of PLAGE

source("source/plage.R")
plage_result <- plage(data, KEGGhsa, metaInfo)

# full functions - differential expression
df <- means_tests(data, metaInfo)
hist(df$pval)
hist(df$corrected_pval)
de_genes <- rownames(df[df$corrected_pval < 0.05,])# labels cause i did it on the data and its original indices

# TMOD ORA

bg <- rownames(data) # all genes
fg <- as.character(de_genes) # foreground
ora_result_tmod <- tmodHGtest(fg, bg, mset = KEGGhsa)
ora_result_tmod

# My implementation of ORA

ora_result <- ora(data, de_genes, KEGGhsa)
significant_ora <- ora_result[ora_result$corrected_pvals < 0.05,]

# CERNO 
# as the tmod implementation of CERNO requires sorted list of the gene names, I sort the df by the p values

ordered_genes <- rownames(df[with(df, order(corrected_pval)), ])

# My implementation of CERNO

cerno_result <- cerno(ordered_genes, KEGGhsa)
significant_cerno <- cerno_result[cerno_result$corrected_pvals < 0.05,]

# TMOD CERNO

cerno_result_tmod <- tmodCERNOtest(ordered_genes, mset = KEGGhsa)
cerno_result_tmod

# My implementation of Z transform

z_result <- cerno(ordered_genes, KEGGhsa, combining_method = "z-transform")
significant_z <- z_result[z_result$corrected_pvals < 0.05,]

# TMOD Z transform

z_result_tmod <- tmodZtest(ordered_genes, mset = KEGGhsa)
z_result_tmod
