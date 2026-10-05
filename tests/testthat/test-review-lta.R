# Independent review: finite path sums and central derivatives do not call the
# package's density, softmax, forward-backward or score implementations.
review_lta_data <- function() {
  data <- data.frame(id = rep(1:6, each = 4L), time = rep(1:4, 6L),
                     y = sin(seq_len(24L)), z = cos(seq_len(24L)),
                     w = rep(seq_len(6L) / 7, each = 4L),
                     c = factor(rep(c("no", "yes"), 12L)),
                     o = ordered(rep(c("low", "mid", "high"), 8L),
                                 levels = c("low", "mid", "high")),
                     k = rep(c(0, 1, 3, 5), 6L), weight = rep(1:3, each = 8L))
  data[c(3L, 18L), "y"] <- NA_real_
  data[c(7L, 14L), "o"] <- NA
  data[c(8L, 17L), "k"] <- NA_real_
  data[-c(2L, 11L, 23L, 24L), ]
}

review_lta_fit <- function(...) {
  quietly(lta(review_lta_data(), c("y", "c", "o", "k"), "id", 3L,
    time = "time", categorical = "c", ordinal = "o", count = "k", missing = "fiml",
    transition_covariates = "z", initial_covariates = "w", occasions = "grid",
    n_starts = 1L, max_iter = 0L, seed = 51L, ...), c("latents_unconverged", "latents_boundary"))
}

review_lta_exact <- function(theta, fit) {
  p <- .lta_unpack(theta, fit)
  softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
  K <- fit$n_profiles
  vapply(seq_len(fit$n_groups), function(g) {
    T <- fit$layout$span[g]
    paths <- as.matrix(expand.grid(rep(list(seq_len(K)), T)))
    class_likelihood <- vapply(seq_len(fit$n_group_classes), function(h) {
      initial <- softmax(c(drop(fit$designs$initial[g, , drop = FALSE] %*%
        matrix(p$initial[, , h], ncol(fit$designs$initial))), 0))
      sum(apply(paths, 1L, function(path) {
        moves <- vapply(seq_len(T)[-1L], function(t) {
          previous <- path[t - 1L]
          if (fit$stayer[h]) return(as.numeric(path[t] == previous))
          pair <- if (fit$order >= 2L && t >= 3L)
            (path[t - 2L] - 1L) * K + previous else previous
          coefficients <- if (fit$order >= 2L && t >= 3L) p$transition2 else p$transition
          logits <- drop(fit$designs$transition[[t - 1L]][g, , drop = FALSE] %*%
            matrix(coefficients[, , pair, h], ncol(fit$designs$transition[[t - 1L]])))
          destinations <- c(setdiff(seq_len(K), previous), previous)
          softmax(c(logits, 0))[match(path[t], destinations)]
        }, numeric(1))
        density <- vapply(seq_len(T), function(t) {
          row <- fit$layout$slot[g, t]
          if (is.na(row)) return(1)
          b <- p$measurement[[if (length(p$measurement) == 1L) 1L else t]]
          state <- path[t]
          y <- fit$x[row, 1L]
          d <- if (is.na(y)) 1 else dnorm(y, b$means[state, 1L], sqrt(b$variances[state, 1L]))
          c <- fit$codes[row, 1L]
          if (!is.na(c)) d <- d * b$response_probabilities[[1L]][state, c]
          o <- fit$extra_data$ordinal[row, 1L]
          if (!is.na(o)) {
            logits <- c(0, b$ordinal_intercepts[[1L]]) +
              (seq_len(3L) - 1L) * b$ordinal_locations[state, 1L]
            d <- d * softmax(logits)[o]
          }
          k <- fit$extra_data$count[row, 1L]
          if (!is.na(k)) d <- d * if (.latents_negative_binomial(fit$extra_data))
            dnbinom(k, size = 1 / b$count_dispersion[state, 1L], mu = b$count_means[state, 1L]) else
              dpois(k, b$count_means[state, 1L])
          d
        }, numeric(1))
        initial[path[1L]] * prod(moves) * prod(density)
      }))
    }, numeric(1))
    log(sum(class_likelihood * p$group_probabilities))
  }, numeric(1))
}

test_that("all general LTA extensions match independent path likelihoods and scores", {
  skip_on_cran()
  variants <- list(
    list(n_group_classes = 2L),
    list(transitions = "occasion", measurement = "occasion"),
    list(order = 2L, transitions = "occasion", mover_stayer = TRUE),
    list(count_model = "negative_binomial", count_dispersion = "equal"),
    list(count_model = "negative_binomial", count_dispersion = "varying", weights = "weight"))
  lapply(variants, function(arguments) {
    fit <- do.call(review_lta_fit, arguments)
    theta <- coef(fit)
    theta[] <- sin(seq_along(theta)) / 3
    # Evaluate away from an optimized point: no stationary score can hide a
    # missing block. Finite coordinates keep all exact-path masses positive.
    exact <- review_lta_exact(theta, fit)
    expectation <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
      .lta_unpack(theta, fit), fit$occasion_of_row, extra = fit$extra_data)
    expect_equal(unname(expectation$group_log_likelihood), exact, tolerance = 1e-11)
    numeric_scores <- vapply(seq_along(theta), function(k) {
      plus <- minus <- theta
      plus[k] <- plus[k] + 1e-5
      minus[k] <- minus[k] - 1e-5
      (review_lta_exact(plus, fit) - review_lta_exact(minus, fit)) / 2e-5
    }, numeric(fit$n_groups)) * (unname(fit$sampling_weights) %||% 1)
    expect_equal(unname(.lta_group_scores(theta, fit)), numeric_scores, tolerance = 2e-7)
    expect_equal(length(theta), fit$n_parameters)
    expect_equal(.lta_pack(.lta_unpack(theta, fit), fit$variance_model, fit$extra_data),
                 unname(theta), tolerance = 1e-12)
    expect_invisible(.lta_check_data(fit, review_lta_data()))
    changed <- review_lta_data()
    changed$k[1L] <- changed$k[1L] + 1
    expect_error(.lta_check_data(fit, changed), class = "latents_bad_inference_data")
  })
})

test_that("second-order occasion models omit unused logits and report stayer identity", {
  fit <- review_lta_fit(order = 2L, transitions = "occasion", mover_stayer = TRUE)
  first <- names(coef(fit))[startsWith(names(coef(fit)), "transition.")]
  second <- names(coef(fit))[startsWith(names(coef(fit)), "transition2.")]
  expect_false(any(grepl("occasion_3", first, fixed = TRUE)))
  expect_false(any(grepl("occasion_4", first, fixed = TRUE)))
  expect_false(any(grepl("occasion_2", second, fixed = TRUE)))
  table <- subset(get_results(fit, "second_order_transitions"), group_class == "stayers")
  expect_equal(table$probability, as.numeric(table$from == table$to))
  expect_identical(unique(get_results(fit, "transitions")$occasion), 2L)
})

test_that("extended LTA refuses unidentified panels and invalid inference controls", {
  data <- review_lta_data()
  expect_error(lta(data[data$time <= 2, ], "y", "id", 2, time = "time",
    order = 2, missing = "fiml"), class = "latents_bad_transition")
  expect_error(review_lta_fit(n_group_classes = 7L), class = "latents_unidentified")
  fit <- review_lta_fit()
  fit$converged <- TRUE
  fit$boundary <- FALSE
  expect_error(.lta_inference(fit, step = 0), class = "latents_bad_argument")
  fit$measurement[[1L]]$count_dispersion <- matrix(.latents_min_dispersion, 3L, 1L)
  expect_error(.lta_inference(fit), class = "latents_boundary_fit")
})

test_that("transition nesting checks indicator roles, grid and covariance restrictions", {
  base <- list(vars = c("y", "c"), id = "id", time = "time", n_profiles = 2L,
    n_group_classes = 1L, variance_model = "varying", covariance_model = "diagonal",
    categorical = "c", transitions = "homogeneous", transition_covariates = character(),
    initial_covariates = character(), measurement = "invariant", order = 1L)
  fake <- function(...) structure(list(arguments = modifyList(base, list(...))), class = "multilpa_lta")
  expect_false(.lta_nested(fake(), fake(categorical = character())))
  expect_false(.lta_nested(fake(occasions = "observed"), fake(occasions = "grid")))
  expect_false(.lta_nested(fake(model = "VVV"), fake(model = "VVI")))
  expect_true(.lta_nested(fake(model = "EEI"), fake(model = "VVI")))
  expect_true(.lta_nested(fake(model = "VVI"), fake(model = "VVV")))
  expect_true(.lta_nested(fake(count_dispersion = "equal"), fake(count_dispersion = "varying")))
})

test_that("every covariance structure uses the correct observed Gaussian marginal in LTA", {
  data <- data.frame(id = rep(1:12, each = 3L), time = rep(1:3, 12L),
                     y = sin(seq_len(36L)) + seq_len(36L) / 40,
                     x = cos(seq_len(36L)) + sin(seq_len(36L) / 3),
                     z = sin(seq_len(36L) / 7))
  data$y[c(2L, 8L)] <- NA
  data$x[c(4L, 10L)] <- NA
  lapply(.multilpa_structures(), function(model) {
    fit <- quietly(lta(data, c("y", "x"), "id", 2, time = "time",
      transition_covariates = "z", missing = "fiml", model = model,
      n_starts = 1L, max_iter = 0L, seed = 62L), c("latents_unconverged", "latents_boundary"))
    block <- fit$measurement[[1L]]
    block$means <- sweep(block$means, 2L, fit$centers, "-")
    theta <- coef(fit)
    decoded <- .lta_unpack(theta, fit)
    expect_equal(length(theta), fit$n_parameters, info = model)
    expect_equal(unname(decoded$measurement[[1L]]$means), unname(block$means), info = model)
    expect_equal(unname(decoded$measurement[[1L]]$variances), unname(block$variances),
                 tolerance = 1e-10, info = model)
    expect_equal(unname(decoded$measurement[[1L]]$covariances), unname(block$covariances),
                 tolerance = 1e-10, info = model)
    expect_equal(.lta_pack(decoded, fit$variance_model), unname(theta),
                 tolerance = 1e-10, info = model)
    reconstructed <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
                                      decoded, fit$occasion_of_row)
    expect_equal(reconstructed$log_likelihood, fit$log_likelihood,
                 tolerance = 1e-10, info = model)
    actual <- .lta_log_density(fit$x, NULL, list(measurement = list(block)),
                              fit$occasion_of_row)
    expected <- vapply(seq_len(2L), function(k) {
      vapply(seq_len(nrow(data)), function(i) {
        observed <- which(!is.na(fit$x[i, ]))
        variance <- if (is.null(block$covariances)) diag(block$variances[k, ], 2L) else
          block$covariances[, , k]
        variance <- variance[observed, observed, drop = FALSE]
        residual <- fit$x[i, observed] - block$means[k, observed]
        -0.5 * (length(observed) * log(2 * pi) +
          as.numeric(determinant(variance, logarithm = TRUE)$modulus) +
          drop(t(residual) %*% solve(variance, residual)))
      }, numeric(1))
    }, numeric(nrow(data)))
    expect_equal(unname(actual), expected, ignore_attr = TRUE, tolerance = 1e-10,
                 info = model)
  })
})

test_that("LTA simulation preserves ordinal types, counts and indicator missingness", {
  data <- review_lta_data()
  lapply(c("poisson", "negative_binomial"), function(model) {
    fit <- review_lta_fit(count_model = model, order = 2L, mover_stayer = TRUE,
                          measurement = "occasion")
    set.seed(64L)
    simulated <- .lta_simulate(fit, data)
    expect_identical(simulated$id, data$id)
    expect_identical(simulated$time, data$time)
    expect_identical(simulated$z, data$z)
    expect_identical(is.na(simulated[fit$vars]), is.na(data[fit$vars]))
    expect_true(is.ordered(simulated$o))
    expect_identical(levels(simulated$o), levels(data$o))
    expect_true(all(simulated$k[!is.na(simulated$k)] >= 0))
    expect_true(all(simulated$k[!is.na(simulated$k)] %% 1 == 0))
  })
})

test_that("LTA table confidence level and both comparison fits are checked", {
  fit <- review_lta_fit()
  expect_error(get_results(fit, "initial_coefficients", level = 1), "level")
  expect_error(parameter_inference(fit, level = NA_real_), "level")
  changed <- review_lta_data()
  changed$z[1L] <- changed$z[1L] + 1
  # Initial-row transition covariates do not enter any transition likelihood.
  changed$z[2L] <- changed$z[2L] + 1
  expect_error(parameter_inference(fit, data = changed), class = "latents_bad_inference_data")
  expect_error(.lta_bootstrap_lrt(fit, fit, review_lta_data(), 0, 1, 1, 1e-8,
                                1L, quote(bootstrap_lrt())), "iter")
})

test_that("LTA comparison validates category labels and actual time values", {
  fit <- review_lta_fit()
  data <- review_lta_data()
  levels(data$c) <- c("renamed_no", "renamed_yes")
  expect_error(.lta_check_data(fit, data), class = "latents_bad_inference_data")
  data <- review_lta_data()
  data$time <- data$time * 2
  expect_error(.lta_check_data(fit, data), class = "latents_bad_inference_data")
})

test_that("negative-binomial LTA optimization respects floors and reaches stationary likelihood", {
  set.seed(77L)
  n <- 120L
  state <- rep(sample(1:2, n, replace = TRUE), each = 4L)
  data <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(1:4, n),
    y = rnorm(n * 4L, c(-2, 2)[state], 0.7),
    k = rnbinom(n * 4L, mu = c(2, 8)[state], size = 2))
  lapply(c("varying", "equal"), function(dispersion) {
    fit <- quietly(lta(data, c("y", "k"), "id", 2, time = "time", count = "k",
      count_model = "negative_binomial", count_dispersion = dispersion,
      transitions = "occasion", n_starts = 1L, max_iter = 5L, seed = 2L),
      "latents_boundary")
    expect_true(fit$converged)
    expect_true(all(fit$measurement[[1L]]$count_dispersion >= .latents_min_dispersion))
    expect_true(all(diff(fit$log_likelihood_history) >= -1e-8))
    expect_lt(max(abs(colSums(.lta_group_scores(coef(fit), fit)))), 1e-3)
  })
})

test_that("a quasi-Newton stop at the maximum (L-BFGS-B code 52) counts as converged", {
  # These datasets end L-BFGS-B with ABNORMAL_TERMINATION_IN_LNSRCH locally
  # (macOS arm64, R 4.5) at points whose scores are ~1e-6 and that a restart
  # improves by < 1e-12; other platforms hit other datasets.
  cases <- list(c(3L, 1L), c(8L, 1L), c(16L, 2L), c(25L, 1L), c(49L, 1L))
  invisible(lapply(cases, function(case) {
    set.seed(case[1L])
    n <- 120L
    state <- rep(sample(1:2, n, replace = TRUE), each = 4L)
    data <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(1:4, n),
      y = rnorm(n * 4L, c(-2, 2)[state], 0.7),
      k = rnbinom(n * 4L, mu = c(2, 8)[state], size = 2))
    fit <- quietly(lta(data, c("y", "k"), "id", 2, time = "time", count = "k",
      count_model = "negative_binomial",
      count_dispersion = c("varying", "equal")[case[2L]],
      transitions = "occasion", n_starts = 1L, seed = 2L), "latents_boundary")
    expect_true(fit$converged)
    expect_lt(max(abs(colSums(.lta_group_scores(coef(fit), fit)))), 1e-3)
  }))
})

test_that("second-order boundary detection reads only the first-order move into occasion two", {
  fit <- review_lta_fit(order = 2L)
  parameters <- .lta_unpack(rep(0, length(coef(fit))), fit)
  e <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs, parameters,
                       fit$occasion_of_row, extra = fit$extra_data)
  e$per_class[[1L]]$log_transition[[2L]][] <- -100
  expect_false(.lta_probability_boundary(e, fit$layout, fit$stayer))
})

test_that("second-order occasion inference reports every free transition coefficient", {
  set.seed(82L)
  n <- 120L
  state <- sample(1:2, 4L * n, replace = TRUE)
  data <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(1:4, n),
                     y = rnorm(4L * n, c(-2, 2)[state], 0.7))
  fit <- lta(data, "y", "id", 2, time = "time", order = 2L, transitions = "occasion",
             n_starts = 1L, max_iter = 15L, seed = 2L)
  expect_true(fit$converged)
  expect_false(fit$boundary)
  table <- parameter_inference(fit, data = data)
  expect_identical(sum(table$block == "transition2"), 8L)
  expect_true(all(is.finite(table$standard_error)))
  expect_identical(names(coef(fit)), rownames(vcov(fit)))
  expect_equal(fit$n_parameters, 15L)
})

test_that("transition BLRT validates data against the alternative as well as the null", {
  fit <- review_lta_fit()
  fit$converged <- TRUE
  alternative <- fit
  alternative$x[1L, 1L] <- alternative$x[1L, 1L] + 1
  expect_error(.lta_bootstrap_lrt(fit, alternative, review_lta_data(), 1L, 1L,
                                1L, 1e-8, 1L, quote(bootstrap_lrt())),
               class = "latents_bad_inference_data")
})

test_that("occasion-prefixed covariates remain free under second-order transitions", {
  expect_identical(.lta_transition_columns(c("(Intercept)", "occasion_3"), TRUE), 1:2)
  expect_identical(.lta_transition_columns(c("occasion_2", "occasion_3", "occasion_score"), TRUE),
                   c(1L, 3L))
})

test_that("occasion covariance charts include mixed ordinal and count coordinates", {
  data <- review_lta_data()
  data$x <- cos(seq_len(nrow(data)))
  lapply(.multilpa_structures(), function(model) {
    fit <- quietly(lta(data, c("y", "x", "c", "o", "k"), "id", 3L, time = "time",
      categorical = "c", ordinal = "o", count = "k", count_model = "negative_binomial",
      missing = "fiml", measurement = "occasion", model = model, n_starts = 1L,
      max_iter = 0L, seed = 51L), c("latents_unconverged", "latents_boundary"))
    theta <- coef(fit)
    decoded <- .lta_unpack(theta, fit)
    expect_equal(length(theta), fit$n_parameters, info = model)
    expect_equal(.lta_pack(decoded, fit$variance_model, fit$extra_data), unname(theta),
                 tolerance = 1e-9, info = model)
    e <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs, decoded,
                          fit$occasion_of_row, extra = fit$extra_data)
    expect_equal(e$log_likelihood, fit$log_likelihood, tolerance = 1e-9, info = model)
  })
})

test_that("numeric occasion covariates after generated intercepts are preserved", {
  terms <- c("occasion_2", "occasion_3", "occasion_4")
  expect_identical(.lta_transition_columns(terms, TRUE, n_occasion_intercepts = 2L), c(1L, 3L))
  expect_identical(.lta_transition_columns(terms, TRUE, 2L, 2L), 2:3)
})

test_that("homogeneous LTA inference data preserve times and category labels", {
  data <- review_lta_data()
  fit <- quietly(lta(data, c("y", "c"), "id", 2, time = "time", categorical = "c",
                      missing = "fiml", n_starts = 1L, max_iter = 0L, seed = 5L),
                 c("latents_unconverged", "latents_boundary"))
  expect_invisible(.multilpa_check_transition_data(fit, data))
  changed <- data
  changed$time <- changed$time * 2
  expect_error(.multilpa_check_transition_data(fit, changed), class = "latents_bad_inference_data")
  changed <- data
  levels(changed$c) <- c("renamed_no", "renamed_yes")
  expect_error(.multilpa_check_transition_data(fit, changed), class = "latents_bad_inference_data")
})

test_that("ordinal/count-only LTA supports all tables and likelihood reconstruction", {
  data <- review_lta_data()
  fit <- quietly(lta(data, c("o", "k"), "id", 2L, time = "time", ordinal = "o",
    count = "k", missing = "fiml", n_starts = 1L, max_iter = 0L, seed = 51L),
    c("latents_unconverged", "latents_boundary"))
  tables <- quietly(get_results(fit, "all"), "latents_no_standard_errors")
  expect_true(all(vapply(tables, is.data.frame, logical(1))))
  expect_identical(nrow(tables$profiles), 0L)
  expect_equal(length(coef(fit)), fit$n_parameters)
  e <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
    .lta_unpack(coef(fit), fit), fit$occasion_of_row, extra = fit$extra_data)
  expect_equal(e$log_likelihood, fit$log_likelihood, tolerance = 1e-11)
  expect_invisible(.lta_check_data(fit, data))
  simulated <- .lta_simulate(fit, data)
  expect_true(is.ordered(simulated$o))
  expect_identical(is.na(simulated[fit$vars]), is.na(data[fit$vars]))
})

test_that("every non-Gaussian LTA measurement exposes tables, coefficients, refit and simulation", {
  data <- review_lta_data()
  specifications <- list(
    list(vars = "o", ordinal = "o"), list(vars = "k", count = "k"),
    list(vars = "c", categorical = "c"),
    list(vars = c("c", "o", "k"), categorical = "c", ordinal = "o", count = "k"))
  lapply(c("invariant", "occasion"), function(measurement) {
    lapply(specifications, function(specification) {
      fit <- quietly(do.call(lta, c(list(data = data, id = "id", n_profiles = 2L,
        time = "time", missing = "fiml", transitions = "occasion", measurement = measurement,
        n_starts = 1L, max_iter = 0L, seed = 12L), specification)),
        c("latents_unconverged", "latents_boundary"))
      tables <- quietly(get_results(fit, "all"), "latents_no_standard_errors")
      expect_true(all(vapply(tables, is.data.frame, logical(1))))
      expect_identical(nrow(tables$profiles), 0L)
      expect_equal(length(coef(fit)), fit$n_parameters)
      parameters <- .lta_unpack(coef(fit), fit)
      e <- .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs, parameters,
                            fit$occasion_of_row, extra = fit$extra_data)
      expect_equal(e$log_likelihood, fit$log_likelihood, tolerance = 1e-11)
      expect_invisible(.lta_check_data(fit, data))
      refit <- quietly(.lta_refit(fit, data, 1L, 0L, 1e-8),
                        c("latents_unconverged", "latents_boundary"))
      expect_identical(refit$categorical, fit$categorical)
      expect_identical(refit$ordinal, fit$ordinal)
      expect_identical(refit$count, fit$count)
      expect_identical(refit$measurement_model, fit$measurement_model)
      simulated <- .lta_simulate(fit, data)
      expect_identical(is.na(simulated[fit$vars]), is.na(data[fit$vars]))
      expect_identical(lapply(simulated[fit$vars], class), lapply(data[fit$vars], class))
    })
  })
})
