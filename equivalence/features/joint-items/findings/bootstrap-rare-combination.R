# FINDING 2: bootstrap inference stops with an internal, unclassed error when
# a resample misses a rare observed combination of a joined pair.
#
# Output of this script (2026-10-20): both the joined pair (levels 0:0 1:0
# 1:1, the combination 1:0 observed twice in 300 rows) and a single item with
# a category observed twice stop with simpleError "the two signatures must
# have the same shape" (from .multilpa_match_order()). A resample that misses
# the category refits an item with one category fewer, whose profiles cannot
# be aligned with the original fit's. So it is the bootstrap's handling of a
# category absent from a resample, not the join; a joined pair meets it far
# more often, since its categories are products of the two items'. Expected:
# the replicate estimated with that category at its bound, or skipped and
# counted as a failed replicate (latents_failed_replicates), not a stop with
# an unclassed error. Failing assertion: tests/testthat/test-features-joint-
# items.R, "bootstrap inference survives a resample without a rare
# combination" (skipped on CRAN).
# Run from the project root:
#   Rscript equivalence/features/joint-items/findings/bootstrap-rare-combination.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
quiet <- function(expr) withCallingHandlers(expr, message = \(m) invokeRestart("muffleMessage"))
set.seed(3)
n <- 300
class <- sample(1:2, n, TRUE)
a <- stats::rbinom(n, 1, c(0.8, 0.2)[class])
b <- a
b[which(a == 1)[1:2]] <- 0   # the combination a = 1, b = 0 occurs twice
data <- data.frame(a = a, b = b, c = stats::rbinom(n, 1, c(0.8, 0.3)[class]),
                   e = stats::rbinom(n, 1, c(0.7, 0.2)[class]))
joined <- quiet(lca(data, c("a", "b", "c", "e"), 2, local_dependence = "a:b", n_starts = 5,
                    seed = 1))
print(joined$categorical_levels[["a:b"]])
outcome <- tryCatch(quiet(parameter_inference(joined, method = "bootstrap", iter = 30,
                                              seed = 1, n_starts = 3)),
                    error = function(e) e)
cat("joined pair:", class(outcome)[1], "-", conditionMessage(outcome), "\n")
# The same rare response in a single categorical item (a three-category item
# whose third category occurs twice):
single <- transform(data, ab = a + b)
print(table(single$ab))
plain <- quiet(lca(single, c("ab", "c", "e"), 2, n_starts = 5, seed = 1))
outcome <- tryCatch(quiet(parameter_inference(plain, method = "bootstrap", iter = 30, seed = 1,
                                              n_starts = 3)),
                    error = function(e) e)
cat("single item:", class(outcome)[1], "-",
    if (inherits(outcome, "error")) conditionMessage(outcome) else "ran", "\n")
