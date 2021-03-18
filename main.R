# Load source files
source("source/expression_testing.R")
source("source/ora.R")
source("source/cerno.R")

library(tmod)

# Load the data

load("data_lung_cancer.RData")

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
