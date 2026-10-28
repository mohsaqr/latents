# Measurement-unit equivariance applies only to transformations preserving
# the covariance family. All fourteen admit a common scale; eight admit a
# separate scale for each indicator. Counts and categories are never rescaled.
# The structure panel compares finite trajectories before numerical plateaus;
# the CORE-1, mixed and replication cases below run to convergence.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

.release_engine_fit <- function(data, vars = names(data), ...) {
  suppressWarnings(lpa(data, vars, n_profiles = 2, n_starts = 3,
                       seed = 3, ...))
}

# `converged`: grand scaling divides by a standard deviation, which is not
# exact in floating point, and an unconverged SQUAREM trajectory amplifies
# that rounding differently on each platform's BLAS (CI, 2026-10-28: even the
# best of three 12-step starts moved by 9e-7). Those fits run to convergence
# instead, and only what convergence fixes is compared: the best likelihood,
# the posteriors and the replication count, not each start's iterations.
.release_engine_compare <- function(raw, scaled, jacobian, converged = FALSE) {
  starts <- scaled$starts
  likelihoods <- intersect(c("log_likelihood", "screened_log_likelihood"), names(starts))
  invisible(lapply(likelihoods, function(name) {
    adjusted <- starts[[name]] - jacobian
    if (converged) {
      expect_equal(max(adjusted), max(raw$starts[[name]]), tolerance = 1e-10)
    } else {
      expect_equal(adjusted, raw$starts[[name]], tolerance = 1e-10)
    }
  }))
  if (!converged) {
    expect_identical(starts[setdiff(names(starts), likelihoods)],
                     raw$starts[setdiff(names(raw$starts), likelihoods)])
  }
  expect_equal(scaled$subject_posteriors, raw$subject_posteriors,
               tolerance = 1e-6)
  expect_identical(scaled$n_best_replicated, raw$n_best_replicated)
  expect_identical(scaled$replication_tolerance, raw$replication_tolerance)
}

invisible(lapply(c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
                   "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"), function(model) {
  test_that(paste("SQUAREM common-unit start tables agree for", model), {
    data <- srl[seq_len(150L), c("cognitive_strategies", "intrinsic_value",
                                "test_anxiety")]
    raw <- .release_engine_fit(data, model = model, search = "kmeans", max_iter = 12)
    changed <- .release_engine_fit(data * 1000, model = model,
                                   search = "kmeans", max_iter = 12)
    .release_engine_compare(raw, changed, -nrow(data) * ncol(data) * log(1000))
  })
}))

invisible(lapply(c("EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVV", "VVV"), function(model) {
  test_that(paste("SQUAREM grand-scaled start tables agree for", model), {
    data <- srl[seq_len(150L), c("cognitive_strategies", "intrinsic_value",
                                "test_anxiety")]
    data$test_anxiety <- data$test_anxiety * 1000
    raw <- .release_engine_fit(data, model = model, search = "kmeans", tol = 1e-10)
    changed <- .release_engine_fit(data, model = model, search = "kmeans", tol = 1e-10,
                                   centering = "grand", scaling = "grand")
    jacobian <- nrow(data) * sum(log(vapply(data, stats::sd, numeric(1))))
    .release_engine_compare(raw, changed, jacobian, converged = TRUE)
  })
}))

test_that("the CORE-1 competing VVV maxima no longer depend on units", {
  vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy",
            "self_regulation", "test_anxiety")
  data <- srl[vars]
  data$test_anxiety <- data$test_anxiety * 1000
  fit <- function(...) suppressWarnings(lpa(data, vars, 4, model = "VVV",
                                            n_starts = 5, seed = 3, ...))
  raw <- fit()
  changed <- fit(centering = "grand", scaling = "grand")
  jacobian <- nrow(data) * sum(log(vapply(data, stats::sd, numeric(1))))
  .release_engine_compare(raw, changed, jacobian)
})

test_that("mixed and categorical start tables retain unit-free replication", {
  data <- srl[seq_len(150L), c("cognitive_strategies", "intrinsic_value")]
  data$item <- as.integer(data$intrinsic_value > median(data$intrinsic_value))
  raw <- .release_engine_fit(data, categorical = "item")
  changed <- .release_engine_fit(data, categorical = "item",
                                 centering = "grand", scaling = "grand")
  jacobian <- nrow(data) * sum(log(vapply(data[1:2], stats::sd, numeric(1))))
  .release_engine_compare(raw, changed, jacobian)
  categorical <- data.frame(a = data$item,
                            b = as.integer(data$cognitive_strategies > median(data$cognitive_strategies)),
                            c = as.integer(seq_len(nrow(data)) %% 3L))
  fit <- .release_engine_fit(categorical, categorical = names(categorical))
  scaled <- .release_engine_fit(categorical, categorical = names(categorical),
                                scaling = "grand")
  .release_engine_compare(fit, scaled, 0)
  expect_identical(fit$replication_tolerance, 1e-6 * length(as.matrix(categorical)))
})

test_that("replication uses the same absolute gap at different units under plain EM", {
  data <- srl[c("cognitive_strategies", "intrinsic_value", "self_efficacy",
                "self_regulation", "test_anxiety")]
  fit <- function(data) suppressWarnings(lpa(data, names(data), 4, model = "VVI",
    n_starts = 5, seed = 2, tol = 1e-4, acceleration = "none"))
  raw <- fit(data)
  changed <- fit(data * 1e4)
  .release_engine_compare(raw, changed, -nrow(data) * ncol(data) * log(1e4))
  expect_identical(raw$replication_tolerance, 1e-6 * length(as.matrix(data)))
})

test_that("replication counts observed cells and excludes zero-weight rows", {
  data <- srl[seq_len(100L), c("cognitive_strategies", "intrinsic_value")]
  data$item <- as.integer(seq_len(nrow(data)) %% 2L)
  data$cognitive_strategies[seq_len(8L)] <- NA_real_
  data$intrinsic_value[11:15] <- NA_real_
  data$weight <- rep(c(0, 1, 2, 1), length.out = nrow(data))
  vars <- names(data)[1:3]
  fit <- .release_engine_fit(data, vars, categorical = "item", missing = "fiml",
                             weights = "weight")
  expect_identical(fit$replication_tolerance,
                   1e-6 * sum(!is.na(data[data$weight > 0, vars])))
  scaled <- .release_engine_fit(data, vars, categorical = "item", missing = "fiml",
    weights = "weight", centering = "grand", scaling = "grand")
  jacobian <- sum(colSums((!is.na(data[1:2])) * fit$sampling_weights[fit$group_index]) *
                   log(vapply(data[data$weight > 0, 1:2], stats::sd, numeric(1), na.rm = TRUE)))
  .release_engine_compare(fit, scaled, jacobian)
})

test_that("a start that converged slowly to the best maximum counts as replicating it", {
  skip_on_cran()
  activities <- c("time_with_friends", "on_social_media", "tv_video_games",
                  "listened_music", "sports", "walking", "reading",
                  "part_time_job")
  fit <- suppressMessages(lca(student_esm, activities, n_classes = 3, seed = 2))
  starts <- get_results(fit, "starts")
  gaps <- sort(fit$log_likelihood - starts$log_likelihood)
  # The second start stops 7.6e-4 below the best: the same maximum, farther
  # than one EM step at `tol` (2.1e-4) but well inside 1e-6 per value.
  expect_gt(gaps[[2]], fit$tol * length(as.matrix(student_esm[activities])))
  expect_lt(gaps[[2]], fit$replication_tolerance)
  expect_identical(fit$replication_tolerance,
                   1e-6 * length(as.matrix(student_esm[activities])))
  expect_identical(fit$n_best_replicated, 2L)
  coarse <- suppressMessages(lca(student_esm, activities, n_classes = 3,
                                 tol = 1e-6, seed = 2))
  expect_identical(coarse$replication_tolerance, fit$replication_tolerance)
})
