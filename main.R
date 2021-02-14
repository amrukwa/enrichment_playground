# Load source files
source("expression_testing.R")
source("ora.R")

# Load the data

load("data_lung_cancer.RData")

# Look at the data

head(data, n=2L)
head(metaInfo, n=2L)

# Get the p-values for variances and adjust for multiple testing
vpvals <- apply(data, 1, do_ftest, labels=metaInfo)
data$vpvals <- p.adjust(vpvals, method= "BH")

# Get the p-values for means and adjust for multiple testing
diff <- apply(data, 1, do_ttest, labels=metaInfo)
hist(diff)
corrected_diff <- p.adjust(diff, method= "BH")
hist(corrected_diff)

# get the number of genes in dataset (N)
N =  nrow(data)

# get the de genes - labels
de_labels <- which(corrected_diff < 0.05)

# FOR UNIQUE GENE SETS

for (i in 1:length(KEGGhsa)){
  enrichment_pval <- single_ora(N, de_labels, KEGGhsa[i])
  # add to some list again and do bh
}


# TMOD
bg <- rownames(data) # all genes
fg <- as.character(de_labels) # foreground
result_tmod <- tmodHGtest(fg, bg, mset = KEGGhsa) 
result_tmod
