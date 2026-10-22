# The covariance matrices of a pair-covariance fit (and of a full VVV fit)
# are not exactly symmetric: the two triangles can differ in the last bit.
# Run from the project root:
#   Rscript equivalence/features/pair-covariance/findings/asymmetric-covariance.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(11)
n <- 300
class <- sample.int(2L, n, replace = TRUE, prob = c(0.55, 0.45))
rho <- c(0.6, -0.5)[class]
first <- stats::rnorm(n)
second <- rho * first + sqrt(1 - rho^2) * stats::rnorm(n)
third <- stats::rnorm(n)
fourth <- -0.4 * third + sqrt(1 - 0.16) * stats::rnorm(n)
data <- data.frame(y1 = c(0, 2.2)[class] + first, y2 = c(0, 1.8)[class] + 1.2 * second,
                   y3 = c(0, 1.5)[class] + third, y4 = c(0, -1.4)[class] + 0.8 * fourth)
asymmetry <- function(fit) {
  vapply(seq_len(fit$n_profiles), function(j) {
    s <- unname(fit$covariances[, , j])
    max(abs(s - t(s)))
  }, numeric(1))
}
paired <- lpa(data, paste0("y", 1:4), 2, local_dependence = c("y1:y2", "y3:y4"),
              n_starts = 4, seed = 1, tol = 1e-10, max_iter = 5000)
full <- lpa(data, paste0("y", 1:4), 2, model = "VVV", n_starts = 4, seed = 1,
            tol = 1e-10, max_iter = 5000)
cat("pairs: max |S - t(S)| per profile:", format(asymmetry(paired)), "\n")
cat("VVV:   max |S - t(S)| per profile:", format(asymmetry(full)), "\n")
# The pair's covariance read from either triangle: get_results() reads
# covariance[item_1, item_2], so writing the pair the other way round reports
# the other triangle.
forward <- get_results(paired, "local_dependence")$covariance
reversed <- get_results(lpa(data, paste0("y", 1:4), 2, local_dependence = c("y2:y1", "y4:y3"),
                            n_starts = 4, seed = 1, tol = 1e-10, max_iter = 5000),
                        "local_dependence")$covariance
cat("y1:y2 vs y2:y1 reported covariances differ by:", format(forward - reversed), "\n")
