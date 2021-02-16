# check null hypothesis: distribution is normal for both control and disease, return both pvals, starting with c, then d
check_distributions <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  c_dist_pval <- shapiro.test(control_row)$p.value
  d_dist_pval <- shapiro.test(disease_row)$p.value
  c(c_dist_pval, d_dist_pval)
}

# check if both distributions are normal
check_dist_pvals <- function(row){
  if (row[1] > 0.05 && row[2] > 0.05){
    TRUE
  } else {
    FALSE
  }
}

# check all pairs of distributions
check_all_distr <- function(dataset, labels){
  pvals <- apply(dataset, 1, check_distributions, labels)
  corrected_pvals <- p.adjust(pvals)
  by_genes <- matrix(corrected_pvals, ncol=2, byrow = TRUE)
  are_normal <- apply(by_genes, 1, check_dist_pvals)
  are_normal
}

# variance test for single row
do_ftest <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  var_pval <- var.test(control_row, disease_row)$p.value
  var_pval
}

# wilcoxon test for single row
do_wilcoxon <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  mean_pval <- wilcox.test(control_row, disease_row)$p.value
  mean_pval
}

# t-test for single row
do_ttest <- function(row, labels, colname="vpvals"){
  gene_vals = head(row, -1)
  control_row =gene_vals[labels$Group=="c"]
  disease_row = gene_vals[labels$Group=="d"]
  var_are_equal = TRUE
  if (row[colname] < 0.05) {
    var_are_equal = FALSE
  }
  mean_pval <- t.test(control_row, disease_row, var.equal=var_are_equal)$p.value
  mean_pval
}

# test means of the whole dataset
means_tests <- function(dataset, labels){
  dataset$pval <- 0
  dataset$are_normal <- check_all_distr(dataset, labels)

  dataset[which(dataset$are_normal == 1), ]$pval <- apply(dataset[which(dataset$are_normal == 1), -(1:2)], 1, do_ftest, labels) # get rid of pval and are_normal
  dataset[which(dataset$are_normal == 1), ]$pval <- p.adjust(dataset[which(dataset$are_normal == 1), ]$pval)
  dataset[which(dataset$are_normal == 1), ]$pval <- apply(dataset[which(dataset$are_normal == 1), -1], 1, do_ttest, labels, colname="pval") # get rid of are_normal
  
  dataset[which(dataset$are_normal == 0), ]$pval <- apply(dataset[which(dataset$are_normal == 0), -(1:2)], 1, do_wilcoxon, labels) # get rid of pval and are_normal
  dataset$corrected_pval <- p.adjust(dataset$pval)
  
  dataset[, c("pval", "corrected_pval")]
}