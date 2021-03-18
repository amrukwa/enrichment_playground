library(R.matlab)

# exporting geneset to matlab
for (i in 1:length(KEGGhsa)){
  x = append(KEGGhsa[i]$MODULES$ID, KEGGhsa[i]$MODULES$Title)
  x = append(x, KEGGhsa[i]$GENES$ID)
  x = matrix(x, nrow=1, byrow=TRUE)
  x = data.frame(x)
  write.table(x, file = "trial.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}

# exporting data to matlab

group <- metaInfo$Group == 'd'
group <- group*1
dataF <- as.matrix(data)
rownames(dataF) <- NULL
prob <- as.numeric(rownames(data))
writeMat("data.mat", group = group, dataF = dataF, prob = prob)
