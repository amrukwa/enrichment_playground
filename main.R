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

# get the de genes - labels
de_labels <- which(corrected_diff < 0.05)
de_genes <- rownames(data)[de_labels]


# FOR UNIQUE GENE SETS
result <- ora(data, de_genes, KEGGhsa)
result
length(which(result < 0.05))

# TMOD
bg <- rownames(data) # all genes
fg <- as.character(de_genes) # foreground
result_tmod <- tmodHGtest(fg, bg, mset = KEGGhsa)
result_tmod

