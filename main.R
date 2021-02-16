# Load source files
source("source/expression_testing.R")
source("source/ora.R")

# Load the data

load("data_lung_cancer.RData")

# Look at the data

head(data, n=2L)
head(metaInfo, n=2L)

# FOR UNIQUE GENE SETS
result <- ora(data, de_genes, KEGGhsa)
result
length(which(result < 0.05))

# TMOD
bg <- rownames(data) # all genes
fg <- as.character(de_genes) # foreground
result_tmod <- tmodHGtest(fg, bg, mset = KEGGhsa)
result_tmod

# full functions - differential expression
df <- means_tests(data, metaInfo)
hist(df$pval)
hist(df$corrected_pval)
de_genes <- rownames(df[df$corrected_pval < 0.05,])# labels cause i did it on the data and its original indices

ora_result <- ora(data, de_genes, KEGGhsa)
ora_result[ora_result$corrected_pvals < 0.05,]
tmodHGtest(fg, bg, mset = KEGGhsa, order.by ="none", qval=1.1)
