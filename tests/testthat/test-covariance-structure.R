# The constrained diagonal covariance family. What is tested here is the
# statistics, not the plumbing: the estimator against the formula it solves
# and the nesting lattice the six models form. The comparisons against mclust
# live in equivalence/test-mclust.R, which does not ship.

activity <- c("browse", "lectures", "forum_read")

structure_fit <- function(..., n_profiles = 3L) {
  # The constraints are imposed by every M-step, so they hold at any iterate:
  # one start and a short run are enough to test them.
  quietly(multilpa(engagement_small, activity, "student",
                   n_profiles = n_profiles, n_group_classes = 1, n_starts = 1,
                   seed = 1, max_iter = 50, ...))
}

test_that("volume, shape and orientation resolve to all fourteen codes", {
  resolve <- multilpa:::.multilpa_resolve_structure
  spherical <- c(EII = "equal", VII = "varying")
  for (code in names(spherical)) {
    expect_identical(resolve("varying", "diagonal", volume = spherical[[code]],
                             shape = "spherical"), code, info = code)
  }
  # The other twelve are volume x shape x orientation, read off the letters.
  letters3 <- c("EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE",
                "EEV", "VEV", "EVV", "VVV")
  word <- function(letter) if (letter == "E") "equal" else "varying"
  for (code in letters3) {
    pieces <- strsplit(code, "")[[1L]]
    expect_identical(
      resolve("varying", "diagonal", volume = word(pieces[1L]),
              shape = word(pieces[2L]),
              orientation = if (pieces[3L] == "I") "axis" else word(pieces[3L])),
      code, info = code)
  }
  # A spherical shape is a multiple of the identity, so it has no orientation.
  expect_error(resolve("varying", "diagonal", shape = "spherical",
                       orientation = "equal"), class = "multilpa_bad_argument")
  # And an axis-parallel orientation is not a full covariance.
  expect_error(resolve("varying", "full", orientation = "axis"),
               class = "multilpa_bad_argument")
})

test_that("legacy arguments still resolve as they always did", {
  resolve <- multilpa:::.multilpa_resolve_structure
  expect_identical(resolve("varying", "diagonal"), "VVI")
  expect_identical(resolve("equal", "diagonal"), "EEI")
  expect_identical(resolve("varying", "full"), "VVV")
  expect_identical(resolve("equal", "full"), "EEE")
  expect_identical(resolve("varying", "diagonal", volume = "equal"), "EVI")
  expect_identical(resolve("varying", "diagonal", shape = "equal"), "VEI")
  expect_identical(resolve("equal", "diagonal", shape = "varying"), "EVI")
  expect_identical(resolve("varying", "diagonal", shape = "spherical"), "VII")
  expect_identical(resolve("equal", "diagonal", shape = "spherical"), "EII")
  # An unrestricted covariance is the varying-orientation corner of the family,
  # so naming a volume there picks the model that constrains only the volume.
  expect_identical(resolve("varying", "full", volume = "equal"), "EVV")
})

test_that("each M-step solves the formula it claims to", {
  solve <- multilpa:::.multilpa_structure_variances
  # Scatter chosen so the profiles differ in both volume and shape.
  diagonals <- rbind(c(4, 1, 1), c(1, 1, 9))
  weights <- c(10, 20)
  d <- 3L

  # VVI and EEI are the two the package already had.
  expect_equal(solve(diagonals, weights, "VVI", 1e-10), diagonals / weights)
  expect_equal(solve(diagonals, weights, "EEI", 1e-10),
               matrix(colSums(diagonals) / 30, 2L, d, byrow = TRUE))
  # Spherical: one variance over every indicator.
  expect_equal(solve(diagonals, weights, "EII", 1e-10)[1L, 1L],
               sum(diagonals) / (30 * d))
  expect_equal(solve(diagonals, weights, "VII", 1e-10)[, 1L],
               rowSums(diagonals) / (weights * d))

  # EVI: Sigma_k = lambda * A_k with |A_k| = 1, so A_k is the profile's own
  # scatter rescaled to determinant one and lambda averages what was divided
  # out. Written out here rather than called, so the test is the formula.
  determinants <- apply(diagonals, 1L, function(row) prod(row)^(1 / d))
  shapes <- diagonals / determinants
  expect_equal(solve(diagonals, weights, "EVI", 1e-10),
               shapes * (sum(determinants) / 30))
  # Its defining constraint: every shape has determinant one, so every
  # profile's covariance has the same determinant.
  evi <- solve(diagonals, weights, "EVI", 1e-10)
  expect_equal(prod(evi[1L, ]), prod(evi[2L, ]))

  # VEI: shared shape, so the profiles' variances are proportional.
  vei <- solve(diagonals, weights, "VEI", 1e-10)
  expect_equal(vei[1L, ] / vei[2L, ],
               rep(vei[1L, 1L] / vei[2L, 1L], d))
  expect_equal(prod(vei[1L, ] / exp(mean(log(vei[1L, ])))), 1)
})

test_that("relaxing a constraint cannot lower the maximized likelihood", {
  specification <- list(
    EII = list(volume = "equal", shape = "spherical"),
    VII = list(volume = "varying", shape = "spherical"),
    EEI = list(variance_model = "equal"),
    VEI = list(volume = "varying", shape = "equal"),
    EVI = list(volume = "equal", shape = "varying"),
    VVI = list(variance_model = "varying"))
  fit_code <- function(code, start = NULL) {
    quietly(do.call(multilpa, c(list(
      engagement_small, activity, "student", n_profiles = 3L,
      n_group_classes = 1, n_starts = 1,
      seed = 1, max_iter = 100, start = start), specification[[code]])))
  }
  fits <- lapply(stats::setNames(nm = names(specification)), fit_code)
  expect_identical(vapply(fits, `[[`, character(1), "covariance_structure"),
                   c(EII = "EII", VII = "VII", EEI = "EEI", VEI = "VEI",
                     EVI = "EVI", VVI = "VVI"))
  # Every edge of the nesting lattice. The wider model is started from the
  # narrower model's solution, which satisfies the wider constraints, so EM
  # can only raise the likelihood: the inequality holds by construction and
  # does not depend on random starts finding the global maximum.
  nested <- list(c("EII", "VII"), c("EII", "EEI"), c("VII", "VEI"),
                 c("EEI", "VEI"), c("EEI", "EVI"), c("VEI", "VVI"),
                 c("EVI", "VVI"))
  invisible(lapply(nested, function(edge) {
    narrow <- fits[[edge[1L]]]
    wide <- fit_code(edge[2L], start = starting_values(narrow))
    expect_lte(narrow$log_likelihood, wide$log_likelihood + 1e-6,
               label = paste(edge, collapse = " in "))
    expect_lt(narrow$n_parameters, wide$n_parameters)
  }))
})

test_that("a constrained structure refuses inference rather than misreporting it", {
  fit <- structure_fit(volume = "equal", shape = "varying", n_profiles = 2L)
  expect_error(parameter_inference(fit),
               class = "multilpa_unsupported_inference")
  expect_error(vcov(fit), class = "multilpa_unsupported_inference")
  # The two structures the coordinates do cover still work.
  free <- structure_fit(variance_model = "varying", n_profiles = 2L)
  expect_s3_class(parameter_inference(free), "data.frame")
})

test_that("a constrained structure refuses a held variance block", {
  free <- structure_fit(variance_model = "varying", n_profiles = 2L)
  start <- starting_values(free)
  expect_error(
    multilpa(engagement_small, activity, "student", n_profiles = 2,
             n_group_classes = 1, n_starts = 1, seed = 1, start = start,
             fixed = "variances", volume = "equal", shape = "varying"),
    class = "multilpa_bad_argument")
  # Holding only the means leaves the spread free, so it is allowed.
  expect_s3_class(
    multilpa(engagement_small, activity, "student", n_profiles = 2,
             n_group_classes = 1, n_starts = 1, seed = 1, start = start,
             fixed = "means", volume = "equal", shape = "varying"),
    "multilpa")
})

test_that("the structure is reported, and the default one is unchanged", {
  default <- structure_fit(n_profiles = 2L)
  expect_identical(default$covariance_structure, "VVI")
  expect_identical(get_results(default, "model")$covariance_structure, "VVI")
  expect_output(print(default), "VVI", fixed = TRUE)
  # Naming the pair that the old arguments already meant must not change the fit.
  named <- structure_fit(volume = "varying", shape = "varying", n_profiles = 2L)
  expect_equal(named$log_likelihood, default$log_likelihood)
  expect_equal(named$variances, default$variances)
})

test_that("the equal-shape iteration is bounded and says when it stops short", {
  solve <- multilpa:::.multilpa_structure_variances
  diagonals <- rbind(c(4, 1, 1), c(1, 1, 9))
  # One iteration cannot reach the fixed point from this start, and the
  # estimator says so instead of returning the last value as if it were the
  # maximiser.
  expect_warning(solve(diagonals, c(10, 20), "VEI", 1e-10, max_iter = 1L),
                 class = "multilpa_no_converge")
  expect_silent(solve(diagonals, c(10, 20), "VEI", 1e-10))
})

test_that("the structure grid crosses models with class counts", {
  candidates <- quietly(enumerate_classes(
    engagement_small, activity, "student", n_profiles = 2:3,
    n_group_classes = 1, model = c("EEI", "EVI", "EEE"),
    n_starts = 2, seed = 1, max_iter = 50))
  grid <- get_results(candidates, "candidates")
  expect_true("model" %in% names(grid))
  expect_identical(nrow(grid), 6L)
  expect_setequal(unique(grid$model), c("EEI", "EVI", "EEE"))
  # The class counts alone no longer name one candidate, and saying so beats
  # returning whichever came first.
  expect_error(candidate_fit(candidates, n_profiles = 3, n_group_classes = 1),
               class = "multilpa_unknown_candidate")
  fitted <- candidate_fit(candidates, n_profiles = 3, n_group_classes = 1,
                          model = "EVI")
  expect_identical(fitted$covariance_structure, "EVI")
  # And the print method names the column, so the grid is readable.
  expect_output(print(candidates), "model", fixed = TRUE)
})

test_that("the grid refuses a model code this package does not fit", {
  expect_error(
    enumerate_classes(engagement_small, activity, "student", n_profiles = 2,
                      n_group_classes = 1, model = "XYZ", n_starts = 1),
    class = "multilpa_bad_argument")
})

test_that("an ellipsoidal structure keeps a full covariance array", {
  fit <- structure_fit(volume = "equal", shape = "varying",
                       orientation = "varying", n_profiles = 2L)
  expect_identical(fit$covariance_structure, "EVV")
  # The orientation cannot live in a diagonal block, so the fit carries the
  # array whatever `covariance_model` defaulted to.
  expect_identical(fit$covariance_model, "full")
  expect_equal(dim(fit$covariances), c(length(activity), length(activity), 2L))
  covariances <- get_results(fit, "covariances")
  off_diagonal <- covariances$covariance[
    covariances$indicator != covariances$indicator_2]
  expect_true(any(abs(off_diagonal) > 1e-8))
})

test_that("the assignments table carries what modal assignment discards", {
  fit <- structure_fit(n_profiles = 2L)
  assignments <- get_results(fit, "assignments")
  expect_true("uncertainty" %in% names(assignments))
  # One minus the posterior of the profile that was assigned.
  posteriors <- get_results(fit, "posteriors", format = "wide")
  columns <- grep("^posterior_profile_", names(posteriors), value = TRUE)
  expect_equal(assignments$uncertainty,
               1 - apply(posteriors[columns], 1L, max))
  expect_true(all(assignments$uncertainty >= 0 & assignments$uncertainty < 1))
})

test_that("a warm-started M-step lands where a cold one does", {
  set.seed(3)
  n <- 400L; d <- 3L; k <- 3L
  X <- matrix(stats::rnorm(n * d), n, d)
  X[1:150, ] <- X[1:150, ] %*% diag(c(3, 1, 0.5)) + 2
  z <- matrix(stats::runif(n * k), n, k)
  z <- z / rowSums(z)
  weights <- colSums(z)
  means <- t(crossprod(z, X) / weights)
  scatter <- lapply(seq_len(k), function(profile) {
    residuals <- sweep(X, 2L, means[, profile], "-")
    crossprod(residuals, residuals * z[, profile])
  })
  invisible(lapply(c("VEE", "EVE", "VVE"), function(code) {
    cold <- multilpa:::.multilpa_structure_covariances(scatter, weights, code, 1e-12)
    # Seeding with the answer must be a fixed point, and seeding with something
    # else must still arrive at it: a warm start may only save work.
    warm <- multilpa:::.multilpa_structure_covariances(scatter, weights, code,
                                                       1e-12, start = cold)
    elsewhere <- multilpa:::.multilpa_structure_covariances(
      scatter, weights, code, 1e-12, start = cold[, , k:1, drop = FALSE])
    expect_equal(warm, cold, tolerance = 1e-6, info = code)
    expect_equal(elsewhere, cold, tolerance = 1e-6, info = code)
  }))
})

test_that("capped M-steps are counted and reported once, not once each", {
  reset <- multilpa:::.multilpa_reset_structure_log
  note <- multilpa:::.multilpa_warn_structure
  report <- multilpa:::.multilpa_report_structure_log
  reset()
  expect_silent(report())
  note("VVE", 1000L)
  note("VVE", 1000L)
  note("VVE", 1000L)
  # One warning, carrying the count: an M-step runs hundreds of times per fit,
  # so warning from inside one would say the same thing hundreds of times.
  expect_warning(report(), "3 M-step", class = "multilpa_no_converge")
  reset()
  expect_silent(report())
})

test_that("every structure's defining constraint holds in the fitted parameters", {
  # Counting a structure's parameters and reporting its name proves nothing
  # about what was fitted. Each model is defined by a constraint on the
  # covariances, and that is what is checked here: a fit reported as EEI whose
  # profiles have different variances is the wrong model under the right label.
  indicators <- c("browse", "lectures", "forum_read")
  volume_of <- function(block, d) prod(eigen(block, symmetric = TRUE,
                                             only.values = TRUE)$values)^(1 / d)
  shape_of <- function(block, d) {
    values <- sort(eigen(block, symmetric = TRUE, only.values = TRUE)$values,
                   decreasing = TRUE)
    values / prod(values)^(1 / d)
  }
  invisible(lapply(multilpa:::.multilpa_structures(), function(code) {
    arguments <- multilpa:::.multilpa_structure_arguments(code)
    # The two shared-orientation models legitimately report M-steps that stop
    # at their cap on a fixture this small; that is the package saying so, not
    # a failure, and every other warning still reaches testthat.
    fit <- quietly(do.call(multilpa, c(
      list(engagement_small, indicators, "student", n_profiles = 3,
           n_group_classes = 1, n_starts = 1, seed = 1, max_iter = 50),
      arguments)), c(.multilpa_expected_warnings, "multilpa_no_converge"))
    expect_identical(fit$covariance_structure, code)
    d <- length(indicators)
    blocks <- lapply(seq_len(3L), function(profile) {
      if (is.null(fit$covariances)) diag(fit$variances[profile, ], d) else
        matrix(fit$covariances[, , profile], d, d)
    })
    volumes <- vapply(blocks, volume_of, numeric(1), d = d)
    shapes <- lapply(blocks, shape_of, d = d)
    if (startsWith(code, "E")) {
      expect_equal(volumes, rep(volumes[1L], 3L), tolerance = 1e-6, info = code)
    }
    if (substr(code, 2L, 2L) == "E") {
      expect_equal(shapes[[2L]], shapes[[1L]], tolerance = 1e-6, info = code)
      expect_equal(shapes[[3L]], shapes[[1L]], tolerance = 1e-6, info = code)
    }
    if (code %in% c("EII", "VII")) {
      # Spherical: every indicator has the same variance within a profile.
      invisible(lapply(seq_len(3L), function(profile) {
        expect_equal(unname(fit$variances[profile, ]),
                     rep(unname(fit$variances[profile, 1L]), d),
                     tolerance = 1e-8, info = code)
      }))
    }
    if (endsWith(code, "I")) {
      # Axis-parallel: no off-diagonal covariance at all.
      invisible(lapply(blocks, function(block) {
        expect_equal(max(abs(block - diag(diag(block), d))), 0, info = code)
      }))
    }
    if (code %in% c("EEI", "EEE")) {
      expect_equal(blocks[[2L]], blocks[[1L]], tolerance = 1e-8, info = code)
    }
  }))
})

test_that("a structure named by its pieces fits what the old arguments fit", {
  indicators <- c("browse", "lectures", "forum_read")
  common <- list(engagement_small, indicators, "student", n_profiles = 3,
                 n_group_classes = 1, n_starts = 1, seed = 1, max_iter = 50)
  # The four structures `variance_model` and `covariance_model` already named
  # are maximized by a different branch of the M-step from the constrained
  # ones, so the two routes to them have to be checked against each other.
  legacy <- list(
    EEI = list(variance_model = "equal"),
    VVI = list(variance_model = "varying"),
    EEE = list(covariance_model = "full", variance_model = "equal"),
    VVV = list(covariance_model = "full", variance_model = "varying"))
  invisible(lapply(names(legacy), function(code) {
    old <- quietly(do.call(multilpa, c(common, legacy[[code]])))
    new <- quietly(do.call(multilpa, c(
      common, multilpa:::.multilpa_structure_arguments(code))))
    expect_identical(old$covariance_structure, code)
    expect_identical(new$covariance_structure, code)
    expect_equal(new$log_likelihood, old$log_likelihood, info = code)
    expect_equal(new$variances, old$variances, info = code)
    expect_identical(new$n_parameters, old$n_parameters)
  }))
})
