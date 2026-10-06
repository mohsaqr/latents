# General latent transition model: occasion-varying and covariate-dependent
# transitions, covariate initial distribution, occasion-specific measurement,
# second order, covariance structures. The comparisons with an exact sum over
# every profile path, numerical derivatives and external fits (depmixS4,
# LMest, Mplus) are in tests/equivalence/.

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

test_that("the general engine reproduces the homogeneous lta() exactly", {
  skip_on_cran()
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
  skip_on_cran()
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

test_that("covariance structures combine with the extensions and nest", {
  skip_on_cran()
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
