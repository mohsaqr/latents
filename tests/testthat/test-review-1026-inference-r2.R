# Round-2 review of the inference area, 2026-10-26: INF2-1 to INF2-6 and
# round-1 INF-3. Each test fails on the code before its fix.

r2_seeded <- function(seed, code) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  force(code)
}

# Two profiles in groups of six; `offsets` adds a unit-level shift that person
# centring removes.
r2_profiles <- function(seed, n_groups = 40L, per_group = 6L, offsets = FALSE) {
  r2_seeded(seed, {
    n <- n_groups * per_group
    profile <- sample(1:2, n, TRUE, c(0.6, 0.4))
    shift <- if (offsets) rep(stats::rnorm(n_groups), each = per_group) else 0
    data.frame(g = rep(seq_len(n_groups), each = per_group),
               a = 5 + shift + stats::rnorm(n, c(-1.2, 1.8)[profile], c(1, 0.8)[profile]),
               b = -3 + shift + stats::rnorm(n, c(1, -1.5)[profile], c(0.9, 1.1)[profile]))
  })
}

# Mean of the four ratios (reported SE of a profile mean over the empirical
# SD of that mean across replications), profiles ordered by their mean on `a`.
r2_calibration <- function(reps, offsets, vcov_type = NULL, ...,
                           generate = r2_profiles, first_seed = 5000L) {
  rows <- lapply(seq_len(reps), function(r) {
    data <- generate(first_seed + r, offsets = offsets)
    fit <- multilpa(data, c("a", "b"), "g", 2, 1, n_starts = 2, seed = r, tol = 1e-9, ...)
    table <- if (is.null(vcov_type)) parameter_inference(fit) else
      parameter_inference(fit, vcov_type = vcov_type)
    order <- order(fit$means[, "a"])
    errors <- stats::setNames(table$standard_error, .multilpa_parameter_names(table))
    wanted <- sprintf("measurement.mean.profile_%d.%s", rep(order, 2), rep(c("a", "b"), each = 2))
    c(as.vector(fit$means[order, ]), errors[wanted])
  })
  values <- do.call(rbind, rows)
  mean(colMeans(values[, 5:8]) / apply(values[, 1:4], 2, stats::sd))
}

test_that("grand-centred SEs are calibrated against the spread of the reported means (INF2-1)", {
  skip_on_cran()
  # Before the fix the ratios averaged 0.81 (grand centring) and 0.64 (grand
  # scaling) in 200 replications; after, 0.98 and 1.05.
  centred <- r2_calibration(60L, FALSE, centering = "grand")
  expect_gt(centred, 0.85)
  expect_lt(centred, 1.15)
  scaled <- r2_calibration(60L, FALSE, scaling = "grand")
  expect_gt(scaled, 0.85)
  expect_lt(scaled, 1.15)
})

test_that("a robust grand-centred covariance gives the identically-zero contrast no variance (INF2-1)", {
  data <- r2_profiles(1001L)
  fit <- multilpa(data, c("a", "b"), "g", 2, 1, centering = "grand", n_starts = 2,
                  seed = 1, tol = 1e-10)
  # The estimates obey sum_k pi_k mu_ka = 0 exactly; the old covariance gave
  # it a standard error of about 0.1, the mean errors' own size.
  expect_lt(abs(sum(fit$profile_probabilities[1, ] * fit$means[, "a"])), 1e-12)
  table <- parameter_inference(fit, vcov_type = "robust")
  covariance <- attr(table, "covariance")
  weights <- stats::setNames(numeric(nrow(covariance)), rownames(covariance))
  weights[sprintf("measurement.mean.profile_%d.a", 1:2)] <- fit$profile_probabilities[1, ]
  weights[sprintf("profile.probability.profile_%d.group_class_1", 1:2)] <- fit$means[, "a"]
  contrast <- sqrt(max(0, drop(t(weights) %*% covariance %*% weights)))
  expect_lt(contrast, 1e-4 * min(table$standard_error[table$parameter == "mean"]))
  # An untransformed fit's errors are left exactly as they were.
  plain <- multilpa(data, c("a", "b"), "g", 2, 1, n_starts = 2, seed = 1, tol = 1e-10)
  plain_table <- parameter_inference(plain)
  expect_equal(attr(plain_table, "covariance_unconstrained"),
               solve(attr(plain_table, "hessian")), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("person-centred fits default to the sandwich, which is no longer refused (INF2-2)", {
  data <- r2_profiles(31L, offsets = TRUE)
  fit <- multilpa(data, c("a", "b"), "g", 2, 1, centering = "person", n_starts = 2,
                  seed = 1, tol = 1e-10)
  table <- parameter_inference(fit)
  expect_identical(attr(table, "vcov_type"), "robust")
  expect_true(all(is.finite(table$standard_error[table$parameter == "mean"])))
  expect_error(parameter_inference(fit, vcov_type = "observed"),
               class = "latents_unsupported_inference")
  expect_error(parameter_inference(fit, vcov_type = "opg"),
               class = "latents_unsupported_inference")
  # The sandwich is A^-1 B A^-1 times G / (G - 1) whatever the rank of B.
  theta <- unname(coef(fit, "unconstrained"))
  scores <- .multilpa_group_scores(theta, fit$indicator_data, fit, NULL)
  expect_lt(qr(scores)$rank, ncol(scores))
  bread <- solve(attr(table, "hessian"))
  expected <- bread %*% crossprod(scores) %*% bread * 40 / 39
  expect_equal(unname(attr(table, "covariance_unconstrained")), unname(expected),
               tolerance = 1e-6)
})

test_that("person-centred sandwich SEs are calibrated (INF2-2)", {
  skip_on_cran()
  # Observed-information errors averaged 0.83 of the spread in 300
  # replications; the sandwich, which was refused, 0.96.
  ratio <- r2_calibration(60L, TRUE, vcov_type = "robust", centering = "person")
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)
})

test_that("measurement tables keep their errors when a response probability is on its bound (INF2-3)", {
  data <- r2_seeded(11L, {
    n <- 300L
    profile <- sample(1:2, n, TRUE, c(0.55, 0.45))
    frame <- data.frame(g = rep(seq_len(60L), each = 5L),
                        a = stats::rnorm(n, c(-1.5, 1.5)[profile]),
                        b = stats::rnorm(n, c(1, -1)[profile]))
    frame$item <- factor(vapply(profile, function(k) sample(c("x", "y", "z"), 1L,
      prob = list(c(0.55, 0.45, 0), c(0.2, 0.3, 0.5))[[k]]), ""))
    frame
  })
  fit <- multilpa(data, c("a", "b", "item"), "g", 2, 1, categorical = "item",
                  n_starts = 4, seed = 3, tol = 1e-12, max_iter = 20000)
  expect_true(any(unlist(fit$response_probabilities) <= 1e-6))
  profiles <- get_results(fit, "profiles")
  expect_true(all(is.finite(profiles$mean_standard_error)))
  expect_true(all(is.finite(profiles$variance_standard_error)))
  fixed <- parameter_inference(fit, boundary = "fix")
  means <- fixed[fixed$parameter == "mean", ]
  expect_equal(profiles$mean_standard_error,
               means$standard_error[match(paste(profiles$profile, profiles$indicator),
                                          paste(sub("profile_", "", means$outcome), means$term))])
})

test_that("an optimum spans at most the tolerance (INF2-5)", {
  chain <- -100 - (0:29) * 0.9e-4
  optimum <- .multilpa_optimum_index(chain, 1e-4)
  spans <- tapply(chain, optimum, function(values) max(values) - min(values))
  expect_true(all(spans <= 1e-4))
  expect_gt(max(optimum), 1L)
  expect_identical(.multilpa_optimum_index(c(-100, -100.00005, -101, NA, -100.00002,
                                             -101.00004), 1e-4),
                   c(1L, 1L, 2L, NA, 1L, 2L))
  expect_identical(.multilpa_optimum_index(-5, 1), 1L)
})

test_that("bootstrap resamples do not depend on what the refits draw (INF2-6)", {
  data <- data.frame(g = rep(1:10, each = 2), v = 1:20)
  draws <- function(consume) {
    seen <- list()
    r2_seeded(1L, .latents_cluster_bootstrap(
      data, data$g, 10L, "g", 3L,
      refit = function(resampled) {
        seen[[length(seen) + 1L]] <<- resampled$v
        stats::runif(consume)
        list(converged = TRUE)
      },
      extract = function(fit) 0, reference_names = "x"))
    seen
  }
  expect_identical(draws(0), draws(5))
})

test_that("a CSV round trip of the fitting data is accepted by every verb (INF-3)", {
  data <- r2_seeded(4L, data.frame(g = rep(1:30, each = 6), a = stats::rnorm(180, 10),
                                   b = stats::rnorm(180, 3) / 7))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data, path, row.names = FALSE)
  csv <- utils::read.csv(path)
  expect_true(any(csv$b != data$b))
  changed <- csv
  changed$a[5] <- changed$a[5] + 1e-9
  for (centering in c("none", "grand")) {
    fit <- multilpa(data, c("a", "b"), "g", 2, 1, centering = centering, n_starts = 2, seed = 1)
    expect_equal(parameter_inference(fit, data = csv)$standard_error,
                 parameter_inference(fit)$standard_error)
    expect_error(parameter_inference(fit, data = changed),
                 class = "latents_bad_inference_data")
    expect_no_error(suppressWarnings(parameter_inference(
      fit, data = csv, method = "bootstrap", iter = 2, n_starts = 1, seed = 1)))
  }
  person <- multilpa(data, c("a", "b"), "g", 2, 1, centering = "person", n_starts = 2, seed = 1)
  expect_equal(parameter_inference(person, data = csv)$standard_error,
               parameter_inference(person)$standard_error)
  expect_error(parameter_inference(person, data = changed),
               class = "latents_bad_inference_data")
})

# A membership covariate `z` moves units between the two profiles.
r2_covariate_profiles <- function(seed, n_groups = 40L, per_group = 6L, offsets = FALSE) {
  r2_seeded(seed, {
    n <- n_groups * per_group
    z <- stats::rnorm(n)
    profile <- ifelse(stats::runif(n) < stats::plogis(0.4 + 0.8 * z), 1L, 2L)
    shift <- if (offsets) rep(stats::rnorm(n_groups), each = per_group) else 0
    data.frame(g = rep(seq_len(n_groups), each = per_group), z = z,
               a = 5 + shift + stats::rnorm(n, c(-1.2, 1.8)[profile], c(1, 0.8)[profile]),
               b = -3 + shift + stats::rnorm(n, c(1, -1.5)[profile], c(0.9, 1.1)[profile]))
  })
}

test_that("covariate fits: grand-transformed SEs are calibrated (INF2-1)", {
  skip_on_cran()
  # 200 replications: grand centring 0.86 before, 1.02 after; grand scaling
  # 0.61 before, 0.98 after.
  centred <- r2_calibration(60L, FALSE, centering = "grand", profile_covariates = "z",
                            generate = r2_covariate_profiles, first_seed = 7000L)
  expect_gt(centred, 0.85)
  expect_lt(centred, 1.15)
  scaled <- r2_calibration(60L, FALSE, scaling = "grand", profile_covariates = "z",
                           generate = r2_covariate_profiles, first_seed = 7000L)
  expect_gt(scaled, 0.85)
  expect_lt(scaled, 1.15)
})

test_that("covariate fits: person-centred fits use the calibrated sandwich (INF2-2)", {
  skip_on_cran()
  fit <- multilpa(r2_covariate_profiles(7001L, offsets = TRUE), c("a", "b"), "g", 2, 1,
                  profile_covariates = "z", centering = "person", n_starts = 2, seed = 1,
                  tol = 1e-9)
  expect_identical(attr(parameter_inference(fit), "vcov_type"), "robust")
  expect_error(parameter_inference(fit, vcov_type = "observed"),
               class = "latents_unsupported_inference")
  # The sandwich was refused as rank deficient; 200 replications gave 1.09
  # (observed-information errors, the old default, 0.94).
  ratio <- r2_calibration(60L, TRUE, vcov_type = "robust", centering = "person",
                          profile_covariates = "z", generate = r2_covariate_profiles,
                          first_seed = 7000L)
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)
})
