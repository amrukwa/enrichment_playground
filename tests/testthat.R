library(testthat)
source("source/expression_testing.R", chdir = TRUE)

test_that("variances", {
  X <- c(59, 54, 74, 35, 71, 46, 50, 31, 24, 42, 22, 37, 38, 52, 47, 48, 49, 64, 36, 67)
  x_labels <- rep("c", length(X))
  Y <- c(63, 43, 37, 57, 53, 40, 14, 69, 36, 44, 47, 42, 60, 33, 56, 25, 65, 41)
  y_labels <- rep("d", length(Y))
  labels <- data.frame(c(x_labels, y_labels))
  colnames(labels) <- c("Group")
  dataset <- data.frame(matrix(c(X, Y), nrow=1))
  pval <- apply(dataset, 1, do_ftest, labels=labels)
  expected_pval <- var.test(X, Y, alternative="two.sided")$p.value
  expect_equal(pval, expected_pval)
  expect_gte(pval, 0.05)
})

test_that("means", {
  X <- c(59, 54, 74, 35, 71, 46, 50, 31, 24, 42, 22, 37, 38, 52, 47, 48, 49, 64, 36, 67)
  x_labels <- rep("c", length(X))
  Y <- c(63, 43, 37, 57, 53, 40, 14, 69, 36, 44, 47, 42, 60, 33, 56, 25, 65, 41)
  y_labels <- rep("d", length(Y))
  labels <- data.frame(c(x_labels, y_labels))
  colnames(labels) <- c("Group")
  vpval <- var.test(X, Y, alternative="two.sided")$p.value
  dataset <- data.frame(matrix(c(X, Y), nrow=1))
  dataset["vpvals"] <- vpval
  pval <- apply(dataset, 1, do_ttest, labels=labels)
  expected_pval <- t.test(X, Y, alternative="two.sided", var.equal = TRUE)$p.value
  expect_equal(pval, expected_pval)
  expect_gte(pval, 0.05)
})

test_that("distribution_row", {
  X <- c(59, 54, 74, 35, 71, 46, 50, 31, 24, 42, 22, 37, 38, 52, 47, 48, 49, 64, 36, 67)
  x_labels <- rep("c", length(X))
  Y <- c(63, 43, 37, 57, 53, 40, 14, 69, 36, 44, 47, 42, 60, 33, 56, 25, 65, 41)
  y_labels <- rep("d", length(Y))
  labels <- data.frame(c(x_labels, y_labels))
  colnames(labels) <- c("Group")
  dataset <- data.frame(matrix(c(X, Y), nrow=1))
  
  pvals <- apply(dataset, 1, check_distributions, labels=labels)
  expected_pval_c <- shapiro.test(X)$p.value
  expected_pval_d <- shapiro.test(Y)$p.value
  
  expect_equal(pvals[1], expected_pval_c)
  expect_equal(pvals[2], expected_pval_d)
})

test_that("distribution_dataset", {
  
})