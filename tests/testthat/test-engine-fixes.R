# Defects found while recording the engine golden baseline (phase E0 of
# validation/ENGINE_DESIGN.md), each fixed before the consolidation so the
# later phases can be held to identical results.

# Regression tests from past audits; CI runs them on every platform.
skip_on_cran()

fixes_noise_frame <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(1)
  data.frame(a = c(stats::rnorm(100), stats::rnorm(100, 4), stats::runif(5, -20, 20)),
             b = c(stats::rnorm(100), stats::rnorm(100, 3), stats::runif(5, -20, 20)))
}

fixes_noise_fit <- function() {
  lpa(fixes_noise_frame(), c("a", "b"), 2, noise = TRUE, n_starts = 2, seed = 1)
}

test_that("a noise fit summarizes, and refuses the tables it cannot form", {
  fit <- fixes_noise_fit()
  expect_true(any(fit$subject_profiles == 0L))
  expect_s3_class(summary(fit), "summary_multilpa")
  tables <- get_results(fit, "all")
  expect_false(any(c("classification_errors", "bch_weights") %in% names(tables)))
  expect_error(get_results(fit, "classification_errors"),
               class = "latents_unsupported_noise")
  expect_error(get_results(fit, "bch_weights"), class = "latents_unsupported_noise")
  # Invariant: the tables a noise fit does report stay probabilities.
  averages <- get_results(fit, "average_posteriors", level = "individuals")
  sums <- tapply(averages$average_posterior, averages$assigned_class, sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-12)
})

test_that("relabelling a noise fit keeps its noise units at profile 0", {
  fit <- fixes_noise_fit()
  swapped <- .multilpa_permute_profiles(fit, c(2L, 1L))
  noise <- fit$subject_profiles == 0L
  expect_identical(swapped$subject_profiles[noise], fit$subject_profiles[noise])
  expect_identical(swapped$subject_profiles[!noise], 3L - fit$subject_profiles[!noise])
  agreement <- sensitivity(fit, seeds = 1:2)$agreement
  expect_false(anyNA(agreement))
  expect_true(all(agreement > 0.99))
})

test_that("transition parameter_inference() uses the caller's step in its table", {
  skip_on_cran()
  fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
             "student", n_profiles = 2, time = "sequence",
             transitions = "occasion", n_starts = 1, seed = 1)
  expect_s3_class(fit, "multilpa_lta")
  coarse <- parameter_inference(fit, step = 1e-2)
  fine <- parameter_inference(fit, step = 1e-4)
  covariance <- attr(coarse, "covariance_unconstrained")
  expect_equal(coarse$standard_error[coarse$block == "transition"],
               unname(sqrt(diag(covariance)))[startsWith(colnames(covariance),
                                                         "transition.")],
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(coarse$standard_error, fine$standard_error,
                                tolerance = 1e-12)))
})

test_that("a general transition fit refuses a single transition network", {
  skip_if_not_installed("tna")
  fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
             "student", n_profiles = 2, time = "sequence",
             transitions = "occasion", n_starts = 1, seed = 1)
  expect_error(get_tna(fit), class = "latents_unsupported_tna")
  expect_error(get_group_tna(fit), class = "latents_unsupported_tna")
})

test_that("a negative definite additive information is refused with a class", {
  case <- data.frame(group = rep(1:8, each = 4L),
                     y = rep(seq(-2, 3, length.out = 8L), each = 4L) +
                       c(-1, -0.3, 0.3, 1))
  fit <- multilpa(case, "y", "group", n_group_classes = 1, family = "additive",
                  tol = 1e-12, seed = 1, n_starts = 1)
  local_mocked_bindings(.additive_information = function(...) -diag(3))
  expect_error(.additive_inference(fit, "observed"),
               class = "latents_singular_information")
})

test_that("a group covariate varying within groups is refused as bad data", {
  data <- student_esm
  data$varies <- seq_len(nrow(data))
  # The covariate-free stage may sit on a variance bound; that warning is
  # not what is tested here.
  expect_error(withCallingHandlers(
    multilpa(data, c("happy", "worried"), "student", n_profiles = 2,
             n_group_classes = 2, group_covariates = "varies", n_starts = 1,
             seed = 1),
    latents_boundary = \(w) invokeRestart("muffleWarning")),
    class = "latents_bad_data")
})

test_that("a growth mixture is not trapped where the level-only starts lead", {
  skip_on_cran()
  # Classes that differ in slope, with persons scattered in level: every start
  # of the model without random effects splits persons by level, and all of
  # them led the growth fit to a maximum 11 log-likelihood units below the
  # one EM reaches from the truth, with slopes 1.2 and 1.7 for 0 and 3. The
  # start from the persons' own trajectories finds it.
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(4)
  waves <- 0:3
  data <- do.call(rbind, lapply(seq_len(120), function(i) {
    k <- 1L + (i > 60)
    effects <- c(stats::rnorm(1, 0, 2), stats::rnorm(1, 0, 0.4))
    mean <- if (k == 1L) 50 + 3 * waves else 55 + 0 * waves
    data.frame(student = i, wave = waves, trajectory = k,
               score = mean + effects[1L] + effects[2L] * waves +
                 stats::rnorm(length(waves), 0, 1.5))
  }))
  fit <- mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                            class_level = "group", random = "wave",
                            random_covariance = "equal", n_starts = 3, seed = 1)
  truth <- fit$params
  truth$beta[] <- c(50, 3, 55, 0)
  truth$sigma2 <- c(2.25, 2.25)
  truth$random_covariance <- list(diag(c(4, 0.16)))
  truth$gamma[] <- 0
  from_truth <- .growth_em(fit$spec, fit$stats, truth, 2000L, 1e-10)
  expect_gte(fit$log_likelihood, from_truth$expectation$log_likelihood - 1e-4)
  expect_identical(utils::tail(fit$starts$kind, 1L), "trajectories")
  slopes <- get_results(fit, "coefficients")
  expect_equal(sort(subset(slopes, term == "wave")$estimate), c(0, 3), tolerance = 0.1)
})
