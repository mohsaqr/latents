# Round-2 review of the group-class families (validation/review-2026-10-26/
# round-2/additive.md): groups of sampling weight zero are absent from every
# count, entropy, classification table, starting value, refusal and warning
# (ADD2-1 to ADD2-4); collapsed classes are flagged (ADD2-5); a boundary fit's
# documented shortfall (ADD2-6); bootstrap_lrt() needs two replicates (ADD2-7).

# Raw ratings: classes in `means` (H x 2), between variances `between`
# (H x 2), within variances `within` (H x 2), groups of sizes `sizes`.
r2_draw <- function(seed, sizes, means, between, within, prefix = "G") {
  withr::with_seed(seed, {
    n_groups <- length(sizes)
    classes <- sample.int(nrow(means), n_groups, replace = TRUE)
    intercepts <- means[classes, , drop = FALSE] +
      matrix(stats::rnorm(n_groups * 2L), n_groups, 2L) *
      sqrt(between[classes, , drop = FALSE])
    index <- rep(seq_len(n_groups), sizes)
    ratings <- intercepts[index, , drop = FALSE] +
      matrix(stats::rnorm(length(index) * 2L), length(index), 2L) *
      sqrt(within[classes[index], , drop = FALSE])
    data.frame(g = paste0(prefix, index), y1 = ratings[, 1L], y2 = ratings[, 2L],
               stringsAsFactors = FALSE)
  })
}

r2_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# 80 weighted groups, and the same with 40 far-away groups of weight zero.
r2_padded <- function() {
  data <- r2_draw(21L, withr::with_seed(1L, sample(2:8, 80L, replace = TRUE)),
                  means = matrix(c(-1.5, 1.5, 1, -1), 2L),
                  between = matrix(c(0.6, 0.8, 0.5, 0.9), 2L),
                  within = matrix(1, 2L, 2L))
  weights <- withr::with_seed(9L, stats::runif(80L, 0.3, 2.5))
  data$w <- weights[as.integer(sub("G", "", data$g))]
  extra <- r2_draw(77L, rep(c(1L, 12L), 20L), means = matrix(c(4, 6, -3, 5), 2L),
                   between = matrix(0.5, 2L, 2L), within = matrix(1, 2L, 2L),
                   prefix = "Z")
  extra$w <- 0
  list(data = data, padded = rbind(data, extra))
}

r2_fit <- function(data, ...) {
  r2_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                    family = "additive", weights = "w", n_starts = 3L, seed = 1L, ...))
}

test_that("zero-weight groups leave the information ratio and the weak-class warning unchanged (ADD2-1)", {
  sets <- r2_padded()
  fit <- r2_fit(sets$data)
  padded <- r2_fit(sets$padded)
  expect_equal(padded$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  effective <- function(x) {
    .additive_effective_groups(x, .additive_inference(x, "robust", warn_weak = FALSE))
  }
  expect_equal(effective(padded)$information_ratio, effective(fit)$information_ratio,
               tolerance = 1e-8)
  expect_gt(min(effective(fit)$information_ratio), 0.6)
  warning_text <- function(x) {
    inference <- .additive_inference(x, "robust", warn_weak = FALSE)
    tryCatch({
      .additive_warn_weak(x, inference)
      NA_character_
    }, latents_weak_class = conditionMessage)
  }
  expect_false(is.na(warning_text(fit)))
  # Rounded to the printed precision, the two messages are the same text.
  expect_identical(warning_text(padded), warning_text(fit))
})

test_that("zero-weight groups are in no count, entropy or classification table (ADD2-2)", {
  sets <- r2_padded()
  fit <- r2_fit(sets$data)
  padded <- r2_fit(sets$padded)
  fit_row <- get_results(fit, "fit")
  padded_row <- get_results(padded, "fit")
  expect_equal(padded_row$entropy, fit_row$entropy, tolerance = 1e-8)
  expect_identical(padded_row$n_groups, 80L)
  expect_identical(padded_row$n_obs, fit_row$n_obs)
  classes <- r2_quiet(get_results(fit, "group_classes"))
  padded_classes <- r2_quiet(get_results(padded, "group_classes"))
  expect_identical(padded_classes$n_assigned, classes$n_assigned)
  expect_equal(padded_classes$mean_posterior, classes$mean_posterior, tolerance = 1e-8)
  expect_equal(get_results(padded, "classification"),
               get_results(fit, "classification"), tolerance = 1e-8)
  expect_identical(sum(classes$n_assigned), 80L)
  # Per-group tables keep every group: a zero-weight group still has a posterior.
  expect_identical(nrow(get_results(padded, "groups")), 120L)
  expect_output(print(padded), "2 group classes, 80 groups")
})

test_that("identification is checked on the groups of positive weight (ADD2-3)", {
  base <- r2_draw(5L, rep(5L, 40L), means = matrix(c(-1, 1, 0.5, -0.5), 2L),
                  between = matrix(0.5, 2L, 2L), within = matrix(1, 2L, 2L))
  base$w <- 1
  extra <- r2_draw(3L, rep(4L, 5L), means = matrix(0, 2L, 2L),
                   between = matrix(0.5, 2L, 2L), within = matrix(1, 2L, 2L),
                   prefix = "Z")
  extra$w <- 0
  fit_call <- function(data, classes = 2L) {
    multilpa(data, c("y1", "y2"), "g", n_group_classes = classes,
             family = "additive", weights = "w", n_starts = 1L, seed = 1L)
  }
  # (a) An indicator constant in every positive-weight group.
  constant <- base
  constant$y2 <- 1
  expect_error(fit_call(rbind(constant, extra)), class = "latents_bad_data")
  # (b) Every positive-weight group a singleton.
  singletons <- r2_draw(4L, rep(1L, 60L), means = matrix(0, 1L, 2L),
                        between = matrix(0.5, 1L, 2L), within = matrix(1, 1L, 2L))
  singletons$w <- 1
  expect_error(fit_call(rbind(singletons, extra), classes = 1L),
               class = "latents_unidentified")
  # (c) More group classes than positive-weight groups.
  few <- base
  few$w <- as.numeric(few$g %in% c("G1", "G2", "G3"))
  expect_error(fit_call(few, classes = 4L), class = "latents_unidentified")
  # (d) Zero within scatter on y2 in every positive-weight group.
  flat <- base
  flat$y2 <- stats::ave(flat$y1, flat$g) + 0.1
  stats <- .additive_prepare(rbind(flat, extra), c("y1", "y2"), "g", "w")
  expect_identical(stats$zero_within_scatter, c(y1 = FALSE, y2 = TRUE))
  # Without weights the checks are those of every group, as before.
  expect_error(fit_call(transform(constant, w = 1)), class = "latents_bad_data")
})

test_that("starting values ignore zero-weight groups (ADD2-4)", {
  sets <- r2_padded()
  structure <- .additive_structure("additive", "varying")
  stats <- .additive_prepare(sets$data, c("y1", "y2"), "g", "w")
  padded <- .additive_prepare(sets$padded, c("y1", "y2"), "g", "w")
  starts <- function(s) {
    lapply(1:3, function(index) withr::with_seed(11L, .additive_initialize(
      s, 2L, structure, 1e-6, index, kind = if (index == 3L) "partition" else "kmeans")))
  }
  expect_equal(starts(padded), starts(stats), tolerance = 1e-12)
  # One start: the fit reached from it does not depend on the padding.
  one <- function(data) {
    r2_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 3L,
                      family = "additive", weights = "w", n_starts = 1L, seed = 1L))
  }
  expect_equal(one(sets$padded)$log_likelihood, one(sets$data)$log_likelihood,
               tolerance = 1e-10)
})

test_that("collapsed group classes are flagged, warned and refused inference (ADD2-5)", {
  parameters <- list(
    means = matrix(c(1, 1, 2, 2), 2L, dimnames = list(c("a", "b"), c("y1", "y2"))),
    between = matrix(c(0.5, 0.5, 0, 0), 2L),
    within = matrix(c(1, 1, 2, 2), 2L))
  expect_identical(.additive_collapsed_classes(parameters), "a = b")
  # Scale and origin free: the same parameters in other units still collapse.
  moved <- parameters
  moved$means <- 1e6 + 1e3 * moved$means
  moved$between <- 1e6 * moved$between
  moved$within <- 1e6 * moved$within
  expect_identical(.additive_collapsed_classes(moved), "a = b")
  apart <- parameters
  apart$within[2L, 1L] <- 1 + 1e-6
  expect_identical(.additive_collapsed_classes(apart), character())
  expect_identical(.additive_collapsed_classes(lapply(parameters, function(block) {
    block[1L, , drop = FALSE]
  })), character())
  # A dispersion fit to data whose classes differ only in means converges
  # onto two identical classes.
  data <- withr::with_seed(3L, {
    sizes <- sample(2:8, 50L, replace = TRUE)
    means <- matrix(c(-1, 1, 0.5, -0.5), 2L)
    classes <- sample.int(2L, 50L, TRUE, c(0.5, 0.5))
    intercepts <- means[classes, ] + matrix(stats::rnorm(100L), 50L, 2L) *
      sqrt(matrix(c(0.3, 0.6, 0.4, 0.2), 2L)[classes, ])
    index <- rep(seq_len(50L), sizes)
    ratings <- intercepts[index, ] + matrix(stats::rnorm(2L * length(index)), ncol = 2L)
    data.frame(g = paste0("G", index), y1 = ratings[, 1L], y2 = ratings[, 2L])
  })
  conditions <- character()
  fit <- withCallingHandlers(
    multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L, family = "dispersion",
             n_starts = 6L, seed = 2L),
    warning = function(w) {
      conditions <<- c(conditions, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  expect_true(fit$collapsed_classes)
  expect_true("latents_weak_class" %in% conditions)
  expect_error(.additive_inference(fit, "observed"), class = "latents_boundary_fit")
  # EM leaves the pair a little apart on the flat ridge, by an amount that
  # depends on which start wins; the likelihood test flags it from any seed
  # (the parameter test alone caught one seed in twelve).
  flagged <- vapply(c(1L, 3L, 7L), function(seed) {
    collapsed <- quietly(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                                  family = "dispersion", n_starts = 6L, seed = seed),
                         classes = "latents_weak_class")
    isTRUE(collapsed$collapsed_classes)
  }, logical(1))
  expect_true(all(flagged))
})

test_that("a boundary fit stops within the documented 0.002 standard errors (ADD2-6)", {
  data <- r2_draw(14L, rep(6L, 40L), means = matrix(c(2, -1), 1L),
                  between = matrix(0, 1L, 2L), within = matrix(1, 1L, 2L))
  fit <- r2_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1L,
                           family = "additive", n_starts = 1L, seed = 1L))
  stats <- fit$sufficient_statistics
  n <- 6
  closed <- vapply(c("y1", "y2"), function(v) {
    a <- stats$averages[, v]
    s <- stats$scatter[, v]
    mu <- mean(a)
    within <- sum(s) / (length(a) * (n - 1))
    between <- mean((a - mu)^2) - within / n
    if (between < 0) {
      between <- 0
      within <- (sum(s) + n * sum((a - mu)^2)) / (length(a) * n)
    }
    c(mu, within, between)
  }, numeric(3))
  free <- closed[3L, ] > 0
  log_likelihood <- function(v) {
    between <- c(0, 0)
    between[free] <- v[5L]
    .additive_expectation(stats, list(means = matrix(v[1:2], 1L), within = v[3:4],
                                      between = matrix(between, 1L),
                                      weights = 1))$log_likelihood
  }
  skip_if_not(sum(free) == 1L)
  optimum <- c(closed[1L, ], closed[2L, ], closed[3L, free])
  standard_error <- sqrt(diag(solve(-numDeriv::hessian(log_likelihood, optimum))))
  estimate <- c(fit$means, fit$within_variances, fit$between_variances[free])
  expect_true(any(fit$between_zero))
  expect_lt(max(abs(estimate - optimum) / standard_error), 0.002)
})

test_that("bootstrap_lrt() for group-class fits needs at least two replicates (ADD2-7)", {
  data <- r2_draw(8L, rep(4L, 30L), means = matrix(c(-1, 1, 1, -1), 2L),
                  between = matrix(0.4, 2L, 2L), within = matrix(1, 2L, 2L))
  null_model <- r2_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1L,
                                  family = "additive", n_starts = 1L, seed = 1L))
  alternative <- r2_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                                   family = "additive", n_starts = 2L, seed = 1L))
  expect_error(bootstrap_lrt(null_model, alternative, iter = 1L, seed = 1L),
               "at least 2")
})
