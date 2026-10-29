# Run from the repository root; evaluates the actual wrapper on saved results.
library(testthat)
lca_wrapper_env <- new.env(parent = globalenv())
lca_wrapper_env$packages_dir <- "equivalence/lca-packages"
lca_wrapper_env$compare_all_lca <- function() read.csv("equivalence/lca-packages/comparison.csv")
expression <- tail(parse("equivalence/test-lca-packages.R"), 1)[[1]]
eval(expression, lca_wrapper_env)
