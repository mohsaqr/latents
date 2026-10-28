# Astra review: active support, regression degeneracy, and transition inference.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

astra_lta_data <- function() {
  .latents_local_seed(21)
  do.call(rbind, lapply(seq_len(35), function(i) {
    states <- sample(1:2, 3, TRUE)
    data.frame(id = i, time = 1:3, z = stats::rnorm(1), w = stats::rnorm(3),
               y1 = stats::rnorm(3, c(0, 3)[states], .5),
               y2 = stats::rnorm(3, c(0, 2)[states], .5), sampling = 1)
  }))
}

astra_lta_fit <- function(data) {
  suppressWarnings(lta(data, c("y1", "y2"), "id", 2, "time",
                        weights = "sampling", transition_covariates = "w",
                        initial_covariates = "z", n_starts = 1, seed = 1,
                        max_iter = 500))
}

astra_lta_pad <- function(data, value) {
  extra <- data[seq_len(15), ]
  extra$id <- extra$id + 35
  extra$sampling <- 0
  extra$w <- extra$z <- value
  rbind(data, extra)
}

test_that("LTA charts use only covariates contributing to the likelihood (A02)", {
  data <- astra_lta_data()
  original <- astra_lta_fit(data)
  padded <- astra_lta_fit(astra_lta_pad(data, 1e8))
  expect_true(original$converged)
  expect_true(padded$converged)
  expect_equal(padded$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  pack <- function(fit) .lta_pack(.lta_parameters(fit), fit$variance_model, fit$extra_data)
  expect_equal(pack(padded), pack(original), tolerance = 1e-8)
  expect_equal(vcov(padded), vcov(original), tolerance = 1e-7)
  # A trailing occasion with no indicator contributes no transition likelihood.
  tail <- data[data$time == 3, ]
  tail$time <- 4
  tail[c("y1", "y2")] <- NA_real_
  tail$w <- 1e8
  fit_missing <- function(d) suppressWarnings(lta(
    d, c("y1", "y2"), "id", 2, "time", missing = "fiml", weights = "sampling",
    transition_covariates = "w", initial_covariates = "z", n_starts = 1,
    seed = 1, max_iter = 500))
  trailing <- fit_missing(rbind(data, tail))
  expect_equal(trailing$log_likelihood, original$log_likelihood, tolerance = 1e-8)
  expect_equal(pack(trailing), pack(original), tolerance = 1e-7)
})

test_that("excluded LTA probabilities cannot veto inference (A06)", {
  data <- astra_lta_data()
  original <- astra_lta_fit(data)
  padded <- astra_lta_fit(astra_lta_pad(data, 10000))
  expect_false(original$boundary)
  expect_false(padded$boundary)
  expect_equal(vcov(padded), vcov(original), tolerance = 1e-7)
  # The boundary rule still catches an active zero probability.
  terms <- list(per_class = list(list(log_initial = log(matrix(c(0, 1), 1)),
                                     log_transition = list(), log_transition2 = list())))
  layout <- list(slot = matrix(1L, 1, 1), within = matrix(TRUE, 1, 1))
  expect_true(.lta_probability_boundary(terms, layout))
})

test_that("constant raw outcomes can have identified regression residuals (A03)", {
  .latents_local_seed(1027)
  data <- data.frame(y = rep(2, 60), x = seq(-2, 3, length.out = 60),
                     o = stats::rnorm(60), w = 1)
  invisible(lapply(list(y ~ 0 + x, y ~ 1 + offset(o)), function(formula) {
    reference <- stats::lm(formula, data)
    invisible(lapply(c(FALSE, TRUE), function(weighted) {
      fit <- mixture_regression(formula, data, 1, n_starts = 0, vcov_type = "none",
                                 weights = if (weighted) "w" else NULL)
      expect_true(fit$converged)
      expect_equal(as.numeric(fit$params$beta), unname(stats::coef(reference)),
                    tolerance = 1e-10)
      expect_equal(unname(fit$params$sigma2), mean(stats::residuals(reference)^2),
                    tolerance = 1e-10)
      expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
                    tolerance = 1e-10)
    }))
  }))
  expect_error(mixture_regression(y ~ 1, data, 1, n_starts = 0),
                class = "latents_bad_data")
  # The residual-scale check must follow the public empty-sample refusal.
  expect_error(mixture_regression(y ~ x, data[0, ], 1, n_starts = 0),
                class = "latents_bad_data")
  # Varying excluded values cannot hide a constant carrying response.
  extra <- data[1:2, ]; extra$y <- c(-3, 4); extra$w <- 0
  expect_error(mixture_regression(y ~ 1, rbind(data, extra), 1,
                                   weights = "w", n_starts = 0),
                class = "latents_bad_data")
})

test_that("zero-weight-only categories do not change the fitted model (A07)", {
  invisible(lapply(c("character", "factor", "numeric"), function(kind) {
    data <- data.frame(a = rep(c("a", "b"), 30), b = rep(c("x", "y", "z"), 20), w = 1)
    padded <- rbind(data, data.frame(a = "absent-category", b = "x", w = 0))
    if (kind == "factor") {
      data$a <- factor(data$a); padded$a <- factor(padded$a)
    } else if (kind == "numeric") {
      data$a <- match(data$a, c("a", "b"))
      padded$a <- match(padded$a, c("a", "b", "absent-category"))
    }
    fit <- function(d) lca(d, c("a", "b"), 1, n_starts = 1, seed = 1, weights = "w")
    original <- fit(data); actual <- fit(padded)
    expect_equal(actual$n_parameters, original$n_parameters)
    expect_equal(actual$bic, original$bic, tolerance = 1e-12)
    expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-12)
    expect_equal(actual$categorical_levels, original$categorical_levels)
    expect_equal(vcov(actual), vcov(original), tolerance = 1e-8)
    expect_equal(vcov(actual, data = padded), vcov(original), tolerance = 1e-8)
    expect_equal(nrow(predict(actual)), nrow(padded))
    expect_error(predict(actual, padded), class = "latents_bad_data")
    rebuilt <- get_results(actual, "data")
    expect_equal(rebuilt$a, padded$a)
    expect_equal(fit(rebuilt)$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  }))
  # Genuine paired items retain the full declared Cartesian grid, including
  # unobserved combinations; trimming ordinary support must not drop these.
  pair <- data.frame(a = c("a", "a", "b"), b = c("x", "y", "x"))
  joint <- suppressWarnings(lca(pair, c("a", "b"), 1, local_dependence = "a:b",
                                n_starts = 1, seed = 1))
  expect_length(joint$categorical_levels[["a:b"]], 4L)
})

test_that("unidentified LTA bootstrap covariance rows are unavailable (A08)", {
  .latents_local_seed(21)
  data <- do.call(rbind, lapply(seq_len(80), function(i) {
    states <- Reduce(function(previous, t) {
      if (stats::runif(1) < .7) previous else 3L - previous
    }, 1:3, init = sample.int(2, 1), accumulate = TRUE)
    data.frame(id = i, time = 1:4,
               y1 = stats::rnorm(4, c(0, 2)[states], .7),
               y2 = stats::rnorm(4, c(0, 1.5)[states], .7))
  }))
  tail <- data[data$time == 4 & data$id <= 20, ]
  tail$time <- 5; tail[c("y1", "y2")] <- NA_real_
  data <- rbind(data, tail)
  fit <- lta(data, c("y1", "y2"), "id", 2, "time", n_starts = 2, seed = 1,
              missing = "fiml", transitions = "occasion")
  table <- parameter_inference(fit, data = data, method = "bootstrap", iter = 4,
                                n_starts = 1, seed = 45)
  covariance <- attr(table, "covariance")
  expect_true(any(!table$estimated))
  expect_true(all(is.na(covariance[!table$estimated, , drop = FALSE])))
  expect_true(all(is.na(covariance[, !table$estimated, drop = FALSE])))
  expect_equal(unname(diag(covariance)[table$estimated]),
                table$standard_error[table$estimated]^2, tolerance = 1e-12)
})

test_that("absent growth persons cannot contaminate starting values", {
  .latents_local_seed(1027)
  data <- data.frame(id = rep(1:30, each = 4), time = rep(0:3, 30), w = 1)
  data$y <- rep(stats::rnorm(30, sd = 1), each = 4) + .5 * data$time +
    stats::rnorm(120, sd = .5)
  fit <- function(d) suppressWarnings(mixture_regression(
    y ~ time, d, 1, id = "id", class_level = "group", random = ~1,
    n_starts = 0, seed = 5, weights = "w", vcov_type = "none"))
  original <- fit(data)
  extra <- data[seq_len(16), ]; extra$id <- rep(31:34, each = 4); extra$w <- 0
  extra$y <- rep(c(-1e100, 1e100), each = 8)
  padded <- fit(rbind(data, extra))
  expect_true(padded$converged)
  expect_equal(padded$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  expect_equal(padded$params, original$params, tolerance = 1e-10)
  expect_equal(padded$starts, original$starts, tolerance = 1e-10)
  expect_equal(nrow(padded$expectation$posterior), 34L)
})

test_that("covariate and transition category support follows positive weights", {
  .latents_local_seed(731)
  data <- data.frame(x = c(stats::rnorm(80, -2), stats::rnorm(80, 2)),
                     item = rep(c("a", "b"), 80), z = stats::rnorm(160),
                     id = seq_len(160), w = 1)
  extra <- data[seq_len(8), ]; extra$id <- extra$id + 160
  extra$w <- 0; extra$item <- "excluded"
  padded <- rbind(data, extra)
  fit_covariate <- function(d) multilpa(
    d, c("x", "item"), "id", 2, 1, categorical = "item", profile_covariates = "z",
    weights = "w", n_starts = 2, seed = 2)
  original <- fit_covariate(data); actual <- fit_covariate(padded)
  expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(actual$n_parameters, original$n_parameters)
  expect_equal(vcov(actual), vcov(original), tolerance = 1e-8)
  expect_equal(vcov(actual, data = padded), vcov(original), tolerance = 1e-8)
  expect_equal(get_results(actual, "data")$item, padded$item)

  data <- astra_lta_data()
  data$item <- rep(c("a", "b"), length.out = nrow(data))
  extra <- data[seq_len(15), ]; extra$id <- extra$id + 35
  extra$sampling <- 0; extra$item <- "excluded"
  padded <- rbind(data, extra)
  fit_lta <- function(d) lta(
    d, c("y1", "y2", "item"), "id", 2, "time", categorical = "item",
    weights = "sampling", transition_covariates = "w", initial_covariates = "z",
    n_starts = 1, seed = 1, max_iter = 500)
  original <- fit_lta(data); actual <- fit_lta(padded)
  expect_equal(actual$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(actual$n_parameters, original$n_parameters)
  expect_equal(vcov(actual), vcov(original), tolerance = 1e-8)
  expect_no_error(.lta_check_data(actual, padded))
})
