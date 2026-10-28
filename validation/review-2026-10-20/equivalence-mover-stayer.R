pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-expected-warnings.R")
expressions <- parse("tests/equivalence/test-lta-general.R")
selected <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("test_that")) &&
  grepl("mover-stayer: the stayer", as.character(e[[2]]), fixed = TRUE), expressions)
stopifnot(length(selected) == 1L)
frame <- new.env(parent = globalenv())
frame$expect_identical <- function(...) invisible(NULL)
frame$expect_true <- function(...) invisible(NULL)
frame$expect_equal <- function(...) invisible(NULL)
eval(selected[[1]][[3]], envir = frame)
cat("Initial logits used by reference:\n")
print(frame$p$initial)
cat("Original reference shares:\n")
print(frame$shares)
frame$shares <- lapply(1:2, function(h) {
  logits <- c(frame$p$initial[1, 1, h], 0)
  e <- exp(logits - max(logits))
  e / sum(e)
})
cat("Stable reference shares:\n")
print(frame$shares)
assignments <- as.list(selected[[1]][[3]])[-1L]
exact_expression <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
  identical(e[[2]], as.name("exact")), assignments)
original_exact <- frame$exact
eval(exact_expression[[1]], envir = frame)
cat(sprintf("Production LL %.14f\nOriginal reference %.14f\nStable independent reference %.14f\nDifference %.16g\n",
  frame$fit$log_likelihood, original_exact, frame$exact,
  frame$fit$log_likelihood - frame$exact))
stopifnot(is.nan(original_exact), is.finite(frame$exact),
          abs(frame$fit$log_likelihood - frame$exact) < 1e-8)
