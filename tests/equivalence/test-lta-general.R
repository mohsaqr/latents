# General latent transition model: occasion-varying and covariate-dependent
# transitions, covariate initial distribution, occasion-specific measurement,
# second order, covariance structures. References: an exact sum over every
# profile path (independent of the forward-backward code), numerical
# derivatives, external fits (depmixS4, LMest, Mplus) bundled as fixtures.
# The path sums build each group's initial and transition probabilities from
# the fitted coefficients with the package's internal link functions
# (.lta_log_initial, .lta_log_transition); the sum over paths, the densities
# and the mixing are written here.

lta_references <- readRDS(test_path("fixtures", "lta-references.rds"))

lta_path_reference <- "independent: exact sum over every profile path, base R"

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
  parameters <- latents:::.lta_parameters(fit)
  sum(vapply(split(data, data$id), function(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(latents:::.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                              matrix(parameters$initial[, , 1L],
                                                     ncol(fit$designs$initial))))
    move <- lapply(seq_len(occasions)[-1L], function(t) {
      exp(latents:::.lta_log_transition(fit$designs$transition[[t - 1L]][j, , drop = FALSE],
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

# The package's analytic scores (the internal .lta_group_scores behind the
# standard errors of lta()) against numDeriv::grad() of the package's
# likelihood at `theta`. The bound is the original test's: an absolute
# difference of at most 1e-5 (1 + the largest |numerical gradient|), loose
# because the Richardson differences of a likelihood of this size carry an
# error far above machine precision.
lta_check_scores <- function(fit, theta, setting) {
  analytic <- colSums(latents:::.lta_group_scores(theta, fit))
  numeric_gradient <- numDeriv::grad(\(v) latents:::.lta_expectation(
    fit$x, fit$codes, fit$layout, fit$designs, latents:::.lta_unpack(v, fit),
    fit$occasion_of_row)$log_likelihood, theta)
  expect_reference(analytic, numeric_gradient, reference = "numDeriv", covers = "lta()",
                   quantity = sprintf("score at a perturbed point, %s", setting),
                   tolerance = 1e-5 * (1 + max(abs(numeric_gradient))), measure = "absolute",
                   note = "Package side: the internal analytic scores (.lta_group_scores) behind vcov() of lta(); reference: numDeriv::grad() of the package's likelihood (.lta_expectation).")
}

test_that("analytic scores equal numerical derivatives in every extension", {
  skip_on_cran()
  skip_without_reference("numDeriv")
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
  # Each point is drawn from the same seed, as in the original test.
  perturbed <- function(fit) {
    set.seed(2)
    coef(fit) + stats::rnorm(length(coef(fit)), 0, 0.05)
  }
  lta_check_scores(fits$covariates, perturbed(fits$covariates),
                   "two group classes, covariate transitions and initial distribution")
  lta_check_scores(fits$occasion, perturbed(fits$occasion),
                   "occasion-varying transitions and measurement")
  lta_check_scores(fits$second, perturbed(fits$second), "second-order transitions")
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
  mplus_note <- paste("Read from tests/equivalence/fixtures/lta-references.rds,",
                      "element mplus_measurement:", ref$source)
  # Mplus prints the likelihood to three decimals.
  expect_reference(as.numeric(logLik(fit)), ref$log_likelihood, reference = "Mplus",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, occasion-specific binary measurement",
                   tolerance = 5e-3, measure = "absolute", note = mplus_note)
  # A count of free parameters: equal or not.
  expect_reference(attr(logLik(fit), "df"), ref$n_parameters, reference = "Mplus",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "number of parameters, occasion-specific binary measurement",
                   tolerance = 0, measure = "absolute", note = mplus_note)
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
  p <- latents:::.lta_parameters(fit)
  m <- fit$measurement[[1L]]
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  exact <- sum(vapply(split(data, data$id), \(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(latents:::.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                              matrix(p$initial[, , 1L], ncol(fit$designs$initial))))
    moves <- lapply(2:4, \(t) exp(latents:::.lta_log_transition(
      fit$designs$transition[[t - 1L]][j, , drop = FALSE],
      array(p$transition[, , , 1L], dim(p$transition)[1:3])))[1L, , ])
    observed <- !is.na(rows$y)
    log(sum(apply(paths, 1L, \(s) {
      initial[s[1L]] * prod(vapply(2:4, \(t) moves[[t - 1L]][s[t - 1L], s[t]], numeric(1))) *
        prod(stats::dnorm(rows$y[observed], m$means[s[observed], 1L],
                          sqrt(m$variances[s[observed], 1L])))
    })))
  }, numeric(1)))
  expect_reference(as.numeric(logLik(fit)), exact, reference = lta_path_reference,
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, missing indicators (FIML) with covariates",
                   tolerance = 1e-9, measure = "relative")
  skip_without_reference("numDeriv")
  lta_check_scores(fit, coef(fit) + 0.05, "missing indicators (FIML) with covariates")
})


test_that("covariate transitions and initial distribution match the exact path sum", {
  skip_on_cran()
  fit <- lta(lta_small, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "z", initial_covariates = "w",
             n_starts = 3, seed = 1, tol = 1e-10)
  expect_true(fit$converged)
  expect_reference(as.numeric(logLik(fit)), lta_brute_force(fit, lta_small),
                   reference = lta_path_reference, covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, covariate transitions and initial distribution",
                   tolerance = 1e-9, measure = "relative")
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
  skip_on_cran()
  fit <- lta(lta_small, "y", "id", n_profiles = 2, time = "time",
             transitions = "occasion", measurement = "occasion", n_starts = 3,
             seed = 1, tol = 1e-10)
  expect_reference(as.numeric(logLik(fit)), lta_brute_force(fit, lta_small),
                   reference = lta_path_reference, covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, occasion-varying transitions and measurement",
                   tolerance = 1e-9, measure = "relative")
  # 4 occasions x (2 means + 2 variances) + 1 initial + 3 occasions x 2 origins.
  expect_identical(fit$n_parameters, 23L)
  table <- get_results(fit, "transitions")
  expect_setequal(unique(table$occasion), 2:4)
  sums <- tapply(table$probability, paste(table$occasion, table$from), sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-12)
  expect_identical(length(unique(get_results(fit, "profiles")$occasion)), 4L)
})

test_that("second-order transitions match the exact path sum and nest first order", {
  skip_on_cran()
  data <- lta_small
  first <- lta(data, "y", "id", n_profiles = 2, time = "time", n_starts = 3,
               seed = 1, tol = 1e-10)
  second <- lta(data, "y", "id", n_profiles = 2, time = "time", order = 2,
                n_starts = 3, seed = 1, tol = 1e-10)
  expect_gte(second$log_likelihood, first$log_likelihood - 1e-6)
  expect_equal(second$n_parameters, first$n_parameters + 4)
  n_states <- 2L
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  p <- latents:::.lta_parameters(second)
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
  expect_reference(as.numeric(logLik(second)), exact, reference = lta_path_reference,
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, second-order transitions",
                   tolerance = 1e-9, measure = "relative")
  expect_identical(nrow(get_results(second, "second_order_transitions")), 8L)
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
  depmix_note <- paste("Read from tests/equivalence/fixtures/lta-references.rds,",
                       "element depmix_covariate:", ref$source)
  expect_reference(as.numeric(logLik(fit)), ref$log_likelihood, reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, covariate transitions",
                   tolerance = 1e-7, measure = "relative", note = depmix_note)
  # A count of free parameters: equal or not.
  expect_reference(attr(logLik(fit), "df"), ref$n_parameters, reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "number of parameters, covariate transitions",
                   tolerance = 0, measure = "absolute", note = depmix_note)
  # The covariate slope of each move. The two programs label the profiles
  # independently, so the absolute slopes are compared in increasing order.
  # 1e-3 relative, the original tolerance: depmixS4 stops at its own
  # convergence criterion, and slopes converge more slowly than the likelihood.
  expect_reference(sort(abs(unname(fit$transition_coefficients["x", 1L, , 1L]))),
                   sort(abs(unname(ref$slopes))), reference = "depmixS4", covers = "lta()",
                   quantity = "absolute covariate slopes of the two moves, in increasing order",
                   tolerance = 1e-3, measure = "relative", key = c("smaller", "larger"),
                   note = depmix_note)
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
  lmest_note <- paste("Read from tests/equivalence/fixtures/lta-references.rds,",
                      "element lmest_occasion:", ref$source)
  expect_reference(as.numeric(logLik(fit)), ref$log_likelihood, reference = "LMest",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, occasion-varying transitions",
                   tolerance = 1e-8, measure = "relative", note = lmest_note)
  # A count of free parameters: equal or not.
  expect_reference(attr(logLik(fit), "df"), ref$n_parameters, reference = "LMest",
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "number of parameters, occasion-varying transitions",
                   tolerance = 0, measure = "absolute", note = lmest_note)
})

test_that("mover-stayer: the stayer class never moves and matches the exact path sum", {
  skip_on_cran()
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
                 c("latents_boundary", "latents_small_classes"))
  expect_identical(names(fit$group_probabilities), c("group_class_1", "stayers"))
  stayers <- subset(get_results(fit, "transitions"), group_class == "stayers")
  expect_equal(stayers$probability, as.numeric(stayers$from == stayers$to))
  # 4 measurement + 1 class weight + 2 initial (one per class) + 2 mover moves.
  expect_identical(fit$n_parameters, 9L)
  p <- latents:::.lta_parameters(fit)
  m <- fit$measurement[[1L]]
  paths <- as.matrix(expand.grid(rep(list(1:2), occasions)))
  shares <- lapply(1:2, \(h) {
    logits <- c(p$initial[1, 1, h], 0)
    e <- exp(logits - max(logits))
    e / sum(e)
  })
  expect_true(all(is.finite(unlist(shares))))
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
  expect_reference(as.numeric(logLik(fit)), exact, reference = lta_path_reference,
                   covers = c("lta()", "logLik.multilpa_lta()"),
                   quantity = "log likelihood, mover-stayer",
                   tolerance = 1e-9, measure = "relative")
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
