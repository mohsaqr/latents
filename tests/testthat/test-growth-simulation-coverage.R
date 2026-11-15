# The retained growth simulation (validation/growth-simulation.R, its
# registry and outputs; validation/GROWTH_SIMULATION.md): the declared rows
# never change, the generator draws the residual process the model defines,
# and the retained summaries of the correlated-residual and bootstrap cells
# are what their replicates give under the predeclared recovery criteria
# (bias within 3.5 Monte Carlo SE; 95% coverage in the central 99.9%
# binomial range).

skip_on_cran()

growth_sim_path <- function(...) {
  testthat::test_path("..", "..", "validation", ...)
}

growth_sim_read <- function(name) {
  path <- growth_sim_path(name)
  skip_if_not(file.exists(path), "supporting study is outside the installed package")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# One top-level assignment of the study script, evaluated without running it.
growth_sim_definition <- function(name) {
  path <- growth_sim_path("growth-simulation.R")
  skip_if_not(file.exists(path), "supporting study is outside the installed package")
  definition <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name(name))
  }, as.list(parse(path)))[[1L]]
  eval(definition[[3L]], environment())
}

test_that("the growth registry keeps the rows declared on 2026-10-17", {
  registry <- growth_sim_read("growth-simulation-registry.csv")
  declared <- data.frame(
    cell = c("balanced", "unbalanced", "weak", "small_variance"),
    persons = c(300L, 300L, 500L, 300L), times = "0,1,2,3,4,5",
    jitter = c(0, 0.3, 0, 0), missing = c(0, 0.25, 0, 0),
    intercept = c("2,5", "2,5", "2,3.2", "2,5"),
    slope = c("0.6,-0.4", "0.6,-0.4", "0.6,0.1", "0.6,-0.4"),
    # small_variance ran with the base variances (GROWTH_SIMULATION.md).
    variance = "1,0.09,0.3 | 0.6,0.05,-0.2",
    sigma = "0.7,0.8", membership = "0.2,0.8", n_starts = 5L,
    base_seed = c(8e6L, 9e6L, 10e6L, 11e6L), stringsAsFactors = FALSE)
  expect_identical(registry[1:4, names(declared)], declared)
  expect_true(all(registry$study[1:4] == "wald"))
  added <- registry[-(1:4), ]
  expect_setequal(added$cell, c("ar", "brownian", "bootstrap_balanced", "bootstrap_weak"))
  expect_identical(added$residual_correlation[added$study == "correlated"], c("ar", "brownian"))
  expect_true(all(added$inference[added$study == "bootstrap"] == "bootstrap"))
  expect_false(anyDuplicated(registry$base_seed) > 0L)
})

test_that("the study's generator draws the residual process the model defines", {
  draw <- growth_sim_definition("draw")
  quiet_cell <- function(process) {
    list(J = 6000L, times = 0:7, jitter = 0, missing = 0, process = process,
         classes = list(intercept = c(0, 0), slope = c(0, 0), variance = list(0, 0),
                        sigma = c(0, 0)))
  }
  empirical <- function(process) {
    frame <- draw(quiet_cell(process), 42L)
    stats::cov(matrix(frame$y, ncol = 8L, byrow = TRUE))
  }
  times <- 0:7
  # The package's own covariance (lcmm's AR(time) and BM(time)).
  spec_ar <- list(residual_correlation = "ar")
  spec_bm <- list(residual_correlation = "brownian")
  model_ar <- .growth_process_covariance(spec_ar, c(rate = 0.7, sd = 0.6), times)
  model_bm <- .growth_process_covariance(spec_bm, c(sd = 0.3), times)
  expect_equal(empirical(list(type = "ar", rate = 0.7, sd = 0.6)), model_ar,
               tolerance = 0.06, ignore_attr = TRUE)
  expect_equal(empirical(list(type = "brownian", sd = 0.3)), model_bm,
               tolerance = 0.06, ignore_attr = TRUE)
})

# The summary each retained replicate file gives under the declared criteria.
growth_sim_resummarise <- function(replicates, methods, primary) {
  keys <- unique(replicates[c("cell", "block", "class", "term")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(r) {
    rows <- merge(keys[r, ], replicates, by = c("cell", "block", "class", "term"))
    ok <- rows[rows$available, ]
    sd_empirical <- stats::sd(ok$estimate)
    bias <- mean(ok$estimate) - rows$truth[1L]
    range <- stats::qbinom(c(0.0005, 0.9995), nrow(ok), 0.95) / nrow(ok)
    coverage <- vapply(methods, function(m) mean(ok[[paste0("covered_", m)]]), numeric(1))
    names(coverage) <- paste0("coverage_", methods)
    data.frame(keys[r, ], availability = mean(rows$available), n_available = nrow(ok),
               bias_mcse = bias / (sd_empirical / sqrt(nrow(ok))),
               as.list(coverage), coverage_low = range[1L], coverage_high = range[2L],
               pass = mean(rows$available) >= 0.98 &&
                 abs(bias / (sd_empirical / sqrt(nrow(ok)))) <= 3.5 &&
                 coverage[[primary]] >= range[1L] && coverage[[primary]] <= range[2L],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

growth_sim_check_summary <- function(study, methods, primary) {
  replicates <- growth_sim_read(sprintf("growth-simulation-%s-replicates.csv", study))
  summary <- growth_sim_read(sprintf("growth-simulation-%s-summary.csv", study))
  again <- growth_sim_resummarise(replicates, methods, primary)
  summary <- summary[match(do.call(paste, again[c("cell", "block", "class", "term")]),
                           do.call(paste, summary[c("cell", "block", "class", "term")])), ]
  expect_false(anyNA(summary$cell))
  for (column in setdiff(names(again), c("cell", "block", "class", "term"))) {
    expect_equal(summary[[column]], again[[column]], label = column, ignore_attr = TRUE)
  }
  expect_true(all(summary$criterion == primary))
  summary
}

test_that("correlated residuals: the retained cells recover the truth and cover it", {
  summary <- growth_sim_check_summary("correlated", c("observed", "robust"),
                                      "coverage_observed")
  # Every parameter of both cells, the process included: bias within 3.5
  # Monte Carlo SE, observed-information coverage in its binomial range.
  expect_true(all(abs(summary$bias_mcse) <= 3.5))
  expect_true(all(summary$coverage_observed >= summary$coverage_low &
                    summary$coverage_observed <= summary$coverage_high))
  expect_true(all(c("rate", "sd") %in% summary$term[summary$cell == "ar"]))
  expect_true("sd" %in% summary$term[summary$cell == "brownian" &
                                        summary$block == "residual_process"])
  # Brownian motion meets every criterion. The autoregressive cell misses
  # availability only: 8 of 150 fits leave the process on its bound, and
  # inference refuses them (GROWTH_SIMULATION.md).
  expect_true(all(summary$pass[summary$cell == "brownian"]))
  ar <- summary[summary$cell == "ar", ]
  expect_true(all(ar$availability < 0.98 & ar$availability > 0.9))
  expect_false(any(ar$pass))
  replicates <- growth_sim_read("growth-simulation-correlated-replicates.csv")
  reasons <- unique(replicates$reason[!replicates$available])
  expect_identical(reasons, "residual process on its bound")
})
