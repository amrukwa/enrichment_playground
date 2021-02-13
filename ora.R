# Packages
library(tmod)
library(plotly)
library(testthat)

# Load the data

load("data_lung_cancer.RData")

# Look at the data

head(data, n=2L)
head(metaInfo, n=2L)

# df <- tmod2DataFrame(KEGGhsa)
# head(df, n=2L)$feature_id #<- genes from data feature_col

# get the number of genes in dataset (N)

N =  nrow(data)

do_ftest <- function(row){
  control_row = row[metaInfo$Group=="c"]
  disease_row = row[metaInfo$Group=="d"]
  var_pval <- var.test(control_row, disease_row)$p.value
  var_pval
}

do_ttest <- function(row){
  gene_vals = head(row, -1)
  control_row =gene_vals[metaInfo$Group=="c"]
  disease_row = gene_vals[metaInfo$Group=="d"]
  var_are_equal = TRUE
  if (row["vpvals"] < 0.05) {
    var_are_equal = FALSE
  }
  mean_pval <- t.test(control_row, disease_row, var.equal=var_are_equal)$p.value
  mean_pval
}

vpvals <- apply(data, 1, do_ftest)
data$vpvals <- p.adjust(vpvals, method= "BH")

diff <- apply(data, 1, do_ttest)
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

