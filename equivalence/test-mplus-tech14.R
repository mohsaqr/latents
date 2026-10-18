# bootstrap_lrt() against Mplus TECH14, the parametric bootstrapped likelihood
# ratio test (McLachlan, 1987; Nylund, Asparouhov and Muthen, 2007).
#
# - The statistic (2 times the log likelihood difference), the parameter
#   difference and the null model's log likelihood are deterministic, and are
#   compared with TECH14 on four retained Mplus runs (jstats-mplus/). Their
#   p-values are all 0.0000 and test nothing.
# - The p-value is compared on a run where the null hypothesis is true
#   (mplus/tech14-null/: 3 against 4 classes on a Monte Carlo dataset
#   generated from 3, 100 bootstrap draws, p = 0.23). Both p-values are
#   simulation estimates, so they must agree within the Monte Carlo error of
#   their difference.

tech14_reference <- function(path) {
  lines <- read_mplus_lines(path)
  section <- lines[grep("TECHNICAL 14 OUTPUT", lines)[1L]:length(lines)]
  list(h0 = printed_value(section, "H0 Loglikelihood Value"),
       statistic = printed_value(section, "2 Times the Loglikelihood Difference"),
       df = as.integer(printed_value(section, "Difference in the Number of Parameters")),
       p_value = printed_value(section, "Approximate P-Value"),
       draws = as.integer(printed_value(section, "Successful Bootstrap Draws")))
}

tech14_fits <- function(data, vars, k, categorical) {
  fit <- function(classes) {
    if (categorical) {
      return(suppressWarnings(lca(data, vars, n_classes = classes, n_starts = 30,
                                  seed = 1, tol = 1e-10, max_iter = 5000)))
    }
    suppressWarnings(lpa(data, vars, n_profiles = classes, model = "VVI", n_starts = 30,
                         seed = 1, tol = 1e-10, max_iter = 5000))
  }
  list(null = fit(k - 1L), alternative = fit(k))
}

test_that("the bootstrap LRT statistic and parameter difference match TECH14", {
  directory <- "jstats-mplus"
  runs <- list(
    list(stem = "mplus-demo-lca-tech14", vars = paste0("u", 1:5), k = 2L, categorical = TRUE),
    list(stem = "mplus-demo-lca-k3-tech14", vars = NULL, k = 3L, categorical = TRUE),
    list(stem = "mplus-demo-lpa-k3-vvi-tech14", vars = paste0("y", 1:4), k = 3L,
         categorical = FALSE))
  invisible(lapply(runs, function(run) {
    lines <- read_mplus_lines(file.path(directory, paste0(run$stem, ".out.txt")))
    data <- utils::read.table(file.path(directory, paste0(run$stem, ".dat")))
    vars <- run$vars %||% paste0("u", seq_len(ncol(data)))
    names(data) <- vars
    reference <- tech14_reference(file.path(directory, paste0(run$stem, ".out.txt")))
    fits <- tech14_fits(data, vars, run$k, run$categorical)
    statistic <- 2 * (fits$alternative$log_likelihood - fits$null$log_likelihood)
    expect_lt(abs(fits$null$log_likelihood - reference$h0), 5e-4 + 1e-4, label = run$stem)
    expect_lt(abs(statistic - reference$statistic), 2e-3, label = run$stem)
    expect_identical(as.integer(fits$alternative$n_parameters - fits$null$n_parameters),
                     reference$df, label = run$stem)
  }))
})

test_that("bootstrap_lrt()'s p-value agrees with TECH14 where the null is true", {
  directory <- file.path("mplus", "tech14-null")
  reference <- tech14_reference(file.path(directory, "tech14.out"))
  expect_identical(reference$draws, 100L)
  vars <- paste0("u", 1:6)
  data <- utils::read.table(file.path(directory, "tech14.dat"), col.names = c("id", vars))
  fits <- tech14_fits(data, vars, 4L, categorical = TRUE)
  expect_lt(abs(2 * (fits$alternative$log_likelihood - fits$null$log_likelihood) -
                  reference$statistic), 2e-3)
  test <- bootstrap_lrt(fits$null, fits$alternative, data = data, iter = 199,
                        n_starts = 10, seed = 1)
  ours <- as.data.frame(test)
  # Both p-values estimate the same quantity from independent simulations.
  p <- (ours$p_value + reference$p_value) / 2
  margin <- 3 * sqrt(p * (1 - p) * (1 / 200 + 1 / reference$draws))
  expect_lt(abs(ours$p_value - reference$p_value), margin)
  expect_equal(ours$statistic, reference$statistic, tolerance = 2e-3 / reference$statistic)
})
