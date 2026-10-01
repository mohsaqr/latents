# General latent transition model: occasion-varying and covariate-dependent
# transitions, covariate initial distribution, occasion-specific measurement,
# second order, covariance structures. References: an exact sum over every
# profile path (independent of the forward-backward code), numerical
# derivatives, external fits (depmixS4, LMest, Mplus) bundled as fixtures.

lta_references <- readRDS(test_path("fixtures", "lta-references.rds"))

# Small two-profile panel with a time-varying covariate and a group-level one.
lta_small <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(8)
  n <- 60L
  occasions <- 4L
  w <- stats::rnorm(n)
  rows <- do.call(rbind, lapply(seq_len(n), function(j) {
    z <- stats::rnorm(occasions)
    state <- integer(occasions)
    state[1L] <- 1L + (stats::runif(1L) < stats::plogis(0.3 + 0.8 * w[j]))
    for (t in seq_len(occasions)[-1L]) {  # simulation only
      move <- stats::plogis(-1 + 0.9 * z[t])
      state[t] <- if (stats::runif(1L) < move) 3L - state[t - 1L] else state[t - 1L]
    }
    data.frame(id = j, time = seq_len(occasions), z = z, w = w[j],
               y = stats::rnorm(occasions, c(-1, 1)[state], 0.7))
  }))
  rows
})

# Exact likelihood by summing over every profile path of every group.
lta_brute_force <- function(fit, data) {
  n_states <- fit$n_profiles
  occasions <- max(data$time)
  paths <- as.matrix(expand.grid(rep(list(seq_len(n_states)), occasions)))
  parameters <- .lta_parameters(fit)
  sum(vapply(split(data, data$id), function(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                    matrix(parameters$initial[, , 1L],
                                           ncol(fit$designs$initial))))
    move <- lapply(seq_len(occasions)[-1L], function(t) {
      exp(.lta_log_transition(fit$designs$transition[[t - 1L]][j, , drop = FALSE],
                              array(parameters$transition[, , , 1L],
                                    dim(parameters$transition)[1:3])))[1L, , ]
    })
    block <- function(t) fit$measurement[[if (length(fit$measurement) == 1L) 1L else t]]
    log(sum(apply(paths, 1L, function(s) {
      probability <- initial[s[1L]] *
        prod(vapply(seq_len(occasions)[-1L], function(t) move[[t - 1L]][s[t - 1L], s[t]],
                    numeric(1)))
      density <- prod(vapply(seq_len(occasions), function(t) {
        stats::dnorm(rows$y[t], block(t)$means[s[t], 1L], sqrt(block(t)$variances[s[t], 1L]))
      }, numeric(1)))
      probability * density
    })))
  }, numeric(1)))
}

test_that("the general engine reproduces the homogeneous lta() exactly", {
  activity <- c("browse", "lectures")
  old <- lta(course_engagement, activity, "student", n_profiles = 2,
             time = "sequence", n_starts = 3, seed = 1, tol = 1e-12, max_iter = 5000)
  new <- .lta_fit_general(course_engagement, activity, "student", "sequence", 2L, 1L,
                          "varying", 3L, 5000L, 1e-12, 1e-6, 1L, character(), 1e-10,
                          "observed", "homogeneous", character(), character(),
                          "invariant", "likelihood", quote(lta()))
  expect_equal(new$log_likelihood, old$log_likelihood, tolerance = 1e-9)
  expect_equal(new$n_parameters, old$n_parameters)
})

test_that("covariate transitions and initial distribution match the exact path sum", {
  fit <- lta(lta_small, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "z", initial_covariates = "w",
             n_starts = 3, seed = 1, tol = 1e-10)
  expect_true(fit$converged)
  expect_equal(fit$log_likelihood, lta_brute_force(fit, lta_small), tolerance = 1e-9)
  # 2 means + 2 variances + 2 initial + 2 origins x 2 terms.
  expect_identical(fit$n_parameters, 10L)
  coefficients <- get_results(fit, "transition_coefficients")
  expect_identical(nrow(coefficients), 4L)
  expect_true(all(is.finite(coefficients$standard_error)))
  # The simulated slope (0.9) is recovered within its interval for both moves.
  slopes <- subset(coefficients, term == "z")
  expect_true(all(slopes$conf_low < 0.9 & slopes$conf_high > 0.9))
})

test_that("occasion-varying transitions and measurement match the exact path sum", {
  fit <- lta(lta_small, "y", "id", n_profiles = 2, time = "time",
             transitions = "occasion", measurement = "occasion", n_starts = 3,
             seed = 1, tol = 1e-10)
  expect_equal(fit$log_likelihood, lta_brute_force(fit, lta_small), tolerance = 1e-9)
  # 4 occasions x (2 means + 2 variances) + 1 initial + 3 occasions x 2 origins.
  expect_identical(fit$n_parameters, 23L)
  table <- get_results(fit, "transitions")
  expect_setequal(unique(table$occasion), 2:4)
  sums <- tapply(table$probability, paste(table$occasion, table$from), sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-12)
  expect_identical(length(unique(get_results(fit, "profiles")$occasion)), 4L)
})

test_that("second-order transitions match the exact path sum and nest first order", {
  data <- lta_small
  first <- lta(data, "y", "id", n_profiles = 2, time = "time", n_starts = 3,
               seed = 1, tol = 1e-10)
  second <- lta(data, "y", "id", n_profiles = 2, time = "time", order = 2,
                n_starts = 3, seed = 1, tol = 1e-10)
  expect_gte(second$log_likelihood, first$log_likelihood - 1e-6)
  expect_equal(second$n_parameters, first$n_parameters + 4)
  n_states <- 2L
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  p <- .lta_parameters(second)
  pi0 <- exp(c(p$initial[1, 1, 1], 0))
  pi0 <- pi0 / sum(pi0)
  move1 <- function(a, b) {
    m <- stats::plogis(p$transition[1, 1, a, 1])
    if (a == b) 1 - m else m
  }
  move2 <- function(a, b, c) {
    m <- stats::plogis(p$transition2[1, 1, (a - 1L) * n_states + b, 1])
    if (b == c) 1 - m else m
  }
  block <- second$measurement[[1L]]
  exact <- sum(vapply(split(data, data$id), function(rows) {
    rows <- rows[order(rows$time), ]
    log(sum(apply(paths, 1L, function(s) {
      pi0[s[1L]] * move1(s[1L], s[2L]) * move2(s[1L], s[2L], s[3L]) *
        move2(s[2L], s[3L], s[4L]) *
        prod(stats::dnorm(rows$y, block$means[s, 1L], sqrt(block$variances[s, 1L])))
    })))
  }, numeric(1)))
  expect_equal(second$log_likelihood, exact, tolerance = 1e-9)
  expect_identical(nrow(get_results(second, "second_order_transitions")), 8L)
})

test_that("analytic scores equal numerical derivatives in every extension", {
  skip_if_not_installed("numDeriv")
  fits <- list(
    # Two group classes on 60 groups: a class can have a move that never
    # occurs. The scores are checked at a perturbed point regardless.
    covariates = quietly(lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                             n_group_classes = 2, transition_covariates = "z",
                             initial_covariates = "w", n_starts = 2, seed = 1,
                             max_iter = 20), "latents_boundary"),
    occasion = quietly(lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                           transitions = "occasion", measurement = "occasion",
                           n_starts = 2, seed = 1, max_iter = 20), "latents_boundary"),
    second = lta(lta_small, "y", "id", n_profiles = 2, time = "time", order = 2,
                 n_starts = 2, seed = 1, max_iter = 20))
  lapply(fits, \(fit) {
    set.seed(2)
    theta <- coef(fit) + stats::rnorm(length(coef(fit)), 0, 0.05)
    analytic <- colSums(.lta_group_scores(theta, fit))
    numeric_gradient <- numDeriv::grad(\(v) .lta_expectation(
      fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(v, fit),
      fit$occasion_of_row)$log_likelihood, theta)
    expect_lt(max(abs(analytic - numeric_gradient)), 1e-5 * (1 + max(abs(numeric_gradient))))
  })
})

test_that("covariate transitions agree with depmixS4 (bundled fixture)", {
  skip_on_cran()
  ref <- lta_references$depmix_covariate
  items <- paste0("y", 1:4)
  long <- do.call(rbind, lapply(1:4, \(t) {
    frame <- as.data.frame(matrix(ref$responses[, t, ], ncol = 4L))
    names(frame) <- items
    cbind(subject = seq_len(nrow(frame)), occasion = t, x = ref$x, frame)
  }))
  long[items] <- lapply(long[items], factor)
  fit <- lta(long, items, "subject", n_profiles = 2, time = "occasion",
             categorical = items, transition_covariates = "x", n_starts = 5,
             seed = 1, tol = 1e-10)
  expect_equal(fit$log_likelihood, ref$log_likelihood, tolerance = 1e-7)
  expect_identical(fit$n_parameters, ref$n_parameters)
  expect_equal(sort(abs(fit$transition_coefficients["x", 1L, , 1L])),
               sort(abs(ref$slopes)), tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("occasion-varying transitions agree with LMest (bundled fixture)", {
  skip_on_cran()
  ref <- lta_references$lmest_occasion
  M <- dim(ref$responses)[3L]
  items <- paste0("y", seq_len(M))
  long <- do.call(rbind, lapply(seq_len(dim(ref$responses)[2L]), \(t) {
    frame <- as.data.frame(matrix(ref$responses[, t, ], ncol = M))
    names(frame) <- items
    cbind(subject = seq_len(nrow(frame)), occasion = t, frame)
  }))
  long[items] <- lapply(long[items], factor)
  fit <- lta(long, items, "subject", n_profiles = ref$K, time = "occasion",
             categorical = items, transitions = "occasion", n_starts = 5, seed = 1,
             tol = 1e-10)
  expect_equal(fit$log_likelihood, ref$log_likelihood, tolerance = 1e-8)
  expect_identical(fit$n_parameters, ref$n_parameters)
})

test_that("occasion-specific measurement agrees with Mplus (bundled fixture)", {
  skip_on_cran()
  ref <- lta_references$mplus_measurement
  long <- do.call(rbind, lapply(1:2, \(t) {
    frame <- as.data.frame(ref$data[, (t - 1L) * 3L + 1:3])
    names(frame) <- paste0("u", 1:3)
    cbind(person = seq_len(nrow(frame)), occasion = t, frame)
  }))
  long[paste0("u", 1:3)] <- lapply(long[paste0("u", 1:3)], factor)
  fit <- lta(long, paste0("u", 1:3), "person", n_profiles = 2, time = "occasion",
             categorical = paste0("u", 1:3), measurement = "occasion", n_starts = 5,
             seed = 1, tol = 1e-10)
  # Mplus prints the likelihood to three decimals.
  expect_lt(abs(fit$log_likelihood - ref$log_likelihood), 5e-3)
  expect_identical(fit$n_parameters, ref$n_parameters)
})

test_that("covariance structures in lta() nest and reproduce the default", {
  activity <- c("browse", "lectures", "forum_read")
  fits <- lapply(c(default = NA, VVI = "VVI", VEI = "VEI", EEI = "EEI"), \(m) {
    suppressWarnings(lta(course_engagement, activity, "student", n_profiles = 2,
                         time = "sequence", model = if (is.na(m)) NULL else m,
                         n_starts = 3, seed = 1))
  })
  expect_equal(fits$VVI$log_likelihood, fits$default$log_likelihood)
  expect_lte(fits$EEI$log_likelihood, fits$VEI$log_likelihood + 1e-6)
  expect_lte(fits$VEI$log_likelihood, fits$VVI$log_likelihood + 1e-6)
  expect_equal(fits$VEI$n_parameters, fits$VVI$n_parameters - 2)
  expect_error(parameter_inference(fits$VEI), class = "latents_unsupported_inference")
})

test_that("tables, methods and refusals of the general model", {
  fit <- lta(lta_small, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "z", n_starts = 2, seed = 1)
  tables <- get_results(fit, "all")
  expect_true(all(vapply(tables, is.data.frame, logical(1))))
  expect_identical(nrow(tables$assignments), nrow(lta_small))
  expect_identical(names(coef(fit)), rownames(vcov(fit)))
  expect_identical(nobs(fit), 60L)
  expect_output(print(fit), "transitions on z")
  expect_s3_class(parameter_inference(fit), "data.frame")
  expect_error(lta(lta_small, "y", "id", n_profiles = 2, time = "time", order = 3),
               class = "latents_bad_argument")
  expect_error(lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                   transition_covariates = "missing_column"),
               class = "latents_bad_data")
  with_na <- lta_small
  with_na$z[3] <- NA
  expect_error(lta(with_na, "y", "id", n_profiles = 2, time = "time",
                   transition_covariates = "z"), class = "latents_bad_data")
  skip_if_not_installed("ggplot2")
  expect_true(ggplot2::is_ggplot(plot(fit)))
})

test_that("EM alone reaches a stationary point (the quasi-Newton finish masks M-step bugs)", {
  fits <- list(
    covariates = lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                     transition_covariates = "z", initial_covariates = "w",
                     n_starts = 1, seed = 1, max_iter = 5),
    second = lta(lta_small, "y", "id", n_profiles = 2, time = "time", order = 2,
                 transitions = "occasion", n_starts = 1, seed = 1, max_iter = 5))
  lapply(fits, \(fit) {
    em <- .lta_em(fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(coef(fit), fit),
                  fit$occasion_of_row, fit$variance_model, 1e-6, NULL, 1e-10,
                  20000L, 1e-14)
    expect_true(all(diff(em$history) >= -1e-8 * (1 + abs(em$history[-1L]))))
    gradient <- colSums(.lta_group_scores(.lta_pack(em$parameters, fit$variance_model), fit))
    expect_lt(max(abs(gradient)), 1e-3)
  })
})

test_that("mover-stayer: the stayer class never moves and matches the exact path sum", {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(11)
  n <- 150L
  occasions <- 4L
  stays <- stats::runif(n) < 0.35
  data <- do.call(rbind, lapply(seq_len(n), function(j) {
    s <- integer(occasions)
    s[1L] <- sample(1:2, 1L)
    for (t in 2:occasions) {  # simulation only
      s[t] <- if (!stays[j] && stats::runif(1L) < 0.3) 3L - s[t - 1L] else s[t - 1L]
    }
    data.frame(id = j, time = seq_len(occasions),
               y = stats::rnorm(occasions, c(-1, 1)[s], 0.8))
  }))
  fit <- quietly(lta(data, "y", "id", n_profiles = 2, time = "time",
                     mover_stayer = TRUE, n_starts = 3, seed = 1, tol = 1e-10),
                 "latents_boundary")
  expect_identical(names(fit$group_probabilities), c("group_class_1", "stayers"))
  stayers <- subset(get_results(fit, "transitions"), group_class == "stayers")
  expect_equal(stayers$probability, as.numeric(stayers$from == stayers$to))
  # 4 measurement + 1 class weight + 2 initial (one per class) + 2 mover moves.
  expect_identical(fit$n_parameters, 9L)
  p <- .lta_parameters(fit)
  m <- fit$measurement[[1L]]
  paths <- as.matrix(expand.grid(rep(list(1:2), occasions)))
  shares <- lapply(1:2, \(h) {
    e <- exp(c(p$initial[1, 1, h], 0))
    e / sum(e)
  })
  move <- function(a, b) {
    mv <- stats::plogis(p$transition[1, 1, a, 1])
    if (a == b) 1 - mv else mv
  }
  exact <- sum(vapply(split(data, data$id), \(rows) {
    log(sum(apply(paths, 1L, \(s) {
      density <- prod(stats::dnorm(rows$y, m$means[s, 1L], sqrt(m$variances[s, 1L])))
      fit$group_probabilities[1L] * shares[[1L]][s[1L]] *
        prod(vapply(2:occasions, \(t) move(s[t - 1L], s[t]), numeric(1))) * density +
        fit$group_probabilities[2L] * shares[[2L]][s[1L]] * all(s == s[1L]) * density
    })))
  }, numeric(1)))
  expect_equal(fit$log_likelihood, exact, tolerance = 1e-9)
  # Stayers and movers that never leave a profile are hard to separate in a
  # small sample, so the fit may sit on a boundary; standard errors follow
  # that flag either way.
  errors <- quietly(get_results(fit, "transition_coefficients"),
                    "latents_no_standard_errors")$standard_error
  if (isTRUE(fit$boundary)) {
    expect_true(all(is.na(errors)))
  } else {
    expect_true(all(is.finite(errors)))
  }
})

test_that("missing indicators (FIML) in the extended model match the exact path sum", {
  data <- lta_small
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(2)
  data$y[sample(nrow(data), 30L)] <- NA
  fit <- lta(data, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "z", initial_covariates = "w", missing = "fiml",
             n_starts = 3, seed = 1, tol = 1e-10)
  p <- .lta_parameters(fit)
  m <- fit$measurement[[1L]]
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  exact <- sum(vapply(split(data, data$id), \(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                    matrix(p$initial[, , 1L], ncol(fit$designs$initial))))
    moves <- lapply(2:4, \(t) exp(.lta_log_transition(
      fit$designs$transition[[t - 1L]][j, , drop = FALSE],
      array(p$transition[, , , 1L], dim(p$transition)[1:3])))[1L, , ])
    observed <- !is.na(rows$y)
    log(sum(apply(paths, 1L, \(s) {
      initial[s[1L]] * prod(vapply(2:4, \(t) moves[[t - 1L]][s[t - 1L], s[t]], numeric(1))) *
        prod(stats::dnorm(rows$y[observed], m$means[s[observed], 1L],
                          sqrt(m$variances[s[observed], 1L])))
    })))
  }, numeric(1)))
  expect_equal(fit$log_likelihood, exact, tolerance = 1e-9)
  skip_if_not_installed("numDeriv")
  theta <- coef(fit) + 0.05
  analytic <- colSums(.lta_group_scores(theta, fit))
  numeric_gradient <- numDeriv::grad(\(v) .lta_expectation(
    fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(v, fit),
    fit$occasion_of_row)$log_likelihood, theta)
  expect_lt(max(abs(analytic - numeric_gradient)), 1e-5 * (1 + max(abs(numeric_gradient))))
})

test_that("covariance structures combine with the extensions and nest", {
  activity <- c("browse", "lectures")
  fits <- lapply(c("VEI", "VVI", "VVV"), \(m) {
    quietly(lta(course_engagement, activity, "student", n_profiles = 2,
                time = "sequence", transition_covariates = "previous_grade",
                model = m, n_starts = 2, seed = 1), "latents_boundary")
  })
  expect_lte(fits[[1L]]$log_likelihood, fits[[2L]]$log_likelihood + 1e-6)
  expect_lte(fits[[2L]]$log_likelihood, fits[[3L]]$log_likelihood + 1e-6)
  expect_equal(fits[[2L]]$n_parameters - fits[[1L]]$n_parameters, 1)
  expect_error(vcov(fits[[1L]]), class = "latents_unsupported_inference")
  expect_error(vcov(fits[[3L]]), class = "latents_unsupported_inference")
})
