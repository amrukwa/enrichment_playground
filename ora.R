# Packages
library(tmod)
library(testthat)
source("expression_testing.R")
# Load the data

load("data_lung_cancer.RData")

# Look at the data

head(data, n=2L)
head(metaInfo, n=2L)


# df <- tmod2DataFrame(KEGGhsa)
# head(df, n=2L)$feature_id #<- genes from data feature_col

# get the number of genes in dataset (N)

N =  nrow(data)

vpvals <- apply(data, 1, do_ftest, labels=metaInfo)
data$vpvals <- p.adjust(vpvals, method= "BH")

diff <- apply(data, 1, do_ttest, labels=metaInfo)
hist(diff)
corrected_diff <- p.adjust(diff, method= "BH")
hist(corrected_diff)

# get the de genes - labels and number (K)
de_labels <- which(corrected_diff < 0.05)
K =  length(de_labels) # significant genes

# UNIQUE GENE SETS
# find the de genes in GS (x) and M - number of genes in GS

M = 82 # genes in gene set (split or sth on the string)
x = 10 # significant in GS (compare how many on de_labels and list)

# Contingency table
c_table = matrix(c(x, K-x, M-x, (N-M)-(K-x)), nrow=2, ncol=2)

# Hypergeometric test
pval = fisher.test(c_table, alternative="greater")$p.value
pval

# TMOD

