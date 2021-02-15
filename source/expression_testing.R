# check null hypothesis: distribution is normal for both control and disease, return both pvals, starting with c, then d
check_distributions <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  c_dist_pval <- shapiro.test(control_row)$p.value
  d_dist_pval <- shapiro.test(disease_row)$p.value
  c(c_dist_pval, d_dist_pval)
}

check_all_distr <- function(dataset, labels){
  # pval correction
  # as_numeric will have to be converted to df with cols c and d and checked
}

# variance test for single row
do_ftest <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  var_pval <- var.test(control_row, disease_row)$p.value
  var_pval
}

# wilcoxon

do_ttest <- function(row, labels){
  gene_vals = head(row, -1) # will be -3 after Shapiro-Wilk
  control_row =gene_vals[labels$Group=="c"]
  disease_row = gene_vals[labels$Group=="d"]
  var_are_equal = TRUE
  if (row["vpvals"] < 0.05) {
    var_are_equal = FALSE
  }
  mean_pval <- t.test(control_row, disease_row, var.equal=var_are_equal)$p.value
  mean_pval
}

# compare means with if both bigger then for_normal -> ftest, ttest else wilcoxon

means_tests <- function(dataset, labels){
  # check all distrs, correct pval
  # think about correct way of storing wilcoxon and ttest -> correct asserting pvals to dataframe later on
  # then for all normal: variances pvals, correct them, calculate means pvals
  # for all non-normal: wilcoxon pvals
  # correction for multiple testing 
  # return both corrected and not corrected -> histograms
}