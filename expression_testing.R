do_ftest <- function(row, labels){
  control_row = row[labels$Group=="c"]
  disease_row = row[labels$Group=="d"]
  var_pval <- var.test(control_row, disease_row)$p.value
  var_pval
}

do_ttest <- function(row, labels){
  gene_vals = head(row, -1)
  control_row =gene_vals[labels$Group=="c"]
  disease_row = gene_vals[labels$Group=="d"]
  var_are_equal = TRUE
  if (row["vpvals"] < 0.05) {
    var_are_equal = FALSE
  }
  mean_pval <- t.test(control_row, disease_row, var.equal=var_are_equal)$p.value
  mean_pval
}