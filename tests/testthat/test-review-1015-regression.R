# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

# ---- Review 2026-10-15: mixture regressions and growth mixtures ----

test_that("a factor predictor can be shared by every class through `common`", {
  shared <- mixture_regression(score ~ hours + strategy, study_hours, n_classes = 1,
                               common = "strategy", n_starts = 0, seed = 1,
                               vcov_type = "none")
  reference <- stats::lm(score ~ hours + strategy, study_hours)
  n <- nrow(study_hours)
  maximum <- -n / 2 * (log(2 * pi * mean(stats::residuals(reference)^2)) + 1)
  expect_equal(shared$log_likelihood, maximum, tolerance = 1e-8)
  expect_true("coefficient.common.strategysurface" %in% names(coef(shared)))
  # A term is matched by its variables, in whatever order they are written.
  crossed <- mixture_regression(score ~ hours * sleep, study_hours, n_classes = 2,
                                common = ~ sleep:hours, n_starts = 1, seed = 1,
                                vcov_type = "none")
  expect_true("coefficient.common.hours:sleep" %in% names(coef(crossed)))
  expect_error(mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                  common = "strategy", n_starts = 1),
               class = "latents_bad_argument")
})

test_that("the random-effects plot of a degenerate class builds the same way each time", {
  skip_if_not_installed("ggplot2")
  fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 40),
                            n_classes = 2, id = "student", class_level = "group",
                            random = "intercept", n_starts = 1, seed = 1,
                            vcov_type = "none")
  # As the fit marks a covariance at its floor.
  fit$degenerate_classes <- "class 2"
  drawn <- plot(fit, what = "random")
  set.seed(11)
  state <- .Random.seed
  built <- ggplot2::ggplot_build(drawn)
  expect_identical(.Random.seed, state)
  expect_identical(built$data, ggplot2::ggplot_build(drawn)$data)
  expect_true(any(vapply(built$data, function(layer) "label" %in% names(layer),
                         logical(1))))
})

test_that("enumerated growth mixtures record a call that names their data", {
  small <- subset(growth_scores, student <= 40)
  grid <- enumerate_regressions(score ~ wave, small, n_classes = 1:2, id = "student",
                                class_level = "group", random = "intercept",
                                n_starts = 1, seed = 1)
  model <- get_results(grid, "model", n_classes = 2)
  expect_identical(model$call[[1L]], quote(mixture_regression))
  expect_identical(model$call$data, quote(small))
})

test_that("a fit does not keep its fitting frame in a formula environment", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1,
                            seed = 1, vcov_type = "none")
  expect_identical(environment(fit$spec$membership), baseenv())
  growth <- mixture_regression(score ~ wave, subset(growth_scores, student <= 40),
                               n_classes = 1, id = "student", class_level = "group",
                               random = "wave", n_starts = 0, seed = 1,
                               vcov_type = "none")
  expect_identical(environment(growth$spec$random), baseenv())
  # The fitting frame (data, every start's result) is not reachable.
  expect_false(exists("results", envir = environment(fit$spec$membership),
                      inherits = FALSE))
})

test_that("confint() honours and validates `parm`", {
  growth <- mixture_regression(score ~ wave, subset(growth_scores, student <= 40),
                               n_classes = 1, id = "student", class_level = "group",
                               random = "intercept", n_starts = 0, seed = 1)
  picked <- confint(growth, parm = "coefficient.class_1.wave")
  expect_identical(picked$parameter, "coefficient.class_1.wave")
  expect_identical(confint(growth, parm = 2)$parameter, names(coef(growth))[2])
  expect_error(confint(growth, parm = "nope"), class = "latents_bad_argument")
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 1, n_starts = 0,
                            seed = 1)
  expect_identical(rownames(confint(fit, parm = 1)), names(coef(fit))[1])
  expect_error(confint(fit, parm = "nope"), class = "latents_bad_argument")
})

test_that("a recovery table keeps a missing known class instead of failing", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, id = "student",
                            class_level = "group", n_starts = 1, seed = 1,
                            vcov_type = "none")
  known <- study_hours
  known$truth <- as.character(known$student_type)
  known$truth[known$student == 1] <- NA
  recovered <- get_results(fit, "recovery", data = known, truth = "truth")
  expect_identical(sum(recovered$n), 150L)
  expect_true(anyNA(recovered$truth))
  # Missing for part of a group only: the known class varies within it.
  known$truth[known$student == 2][1L] <- NA
  expect_error(get_results(fit, "recovery", data = known, truth = "truth"),
               class = "latents_bad_data")
})

test_that("class trajectories are drawn when the classes differ in level only", {
  fit <- mixture_regression(score ~ 1, subset(growth_scores, student <= 40),
                            n_classes = 2, id = "student", class_level = "group",
                            random = "wave", n_starts = 1, seed = 1)
  trajectories <- get_results(fit, "trajectories")
  expect_identical(nrow(trajectories), 120L)
  # No fixed time slope: each class's mean trajectory is flat.
  levels_by_class <- tapply(trajectories$estimate, trajectories$class,
                            function(v) diff(range(v)))
  expect_true(all(levels_by_class < 1e-8))
})

test_that("an enumeration whose every cell fails marks no best row, silently", {
  tiny <- study_hours[seq_len(8), ]
  expect_no_warning(grid <- enumerate_regressions(score ~ hours, tiny, n_classes = 8:9,
                                                  n_starts = 1, seed = 1))
  expect_false(any(as.data.frame(grid)$best_bic))
})

test_that("a single-level ordinal mixture with too few covariate patterns is refused", {
  graded <- study_hours
  graded$grade <- cut(graded$score, c(-Inf, 50, 65, 80, Inf), ordered_result = TRUE)
  # One binary factor: 2 patterns x 3 free category probabilities = 6 < 9.
  expect_error(mixture_regression(grade ~ strategy, graded, family = "ordinal",
                                  n_classes = 2, n_starts = 1, seed = 1),
               class = "latents_not_identified")
  # A continuous predictor, or several outcomes per class assignment, identifies it.
  # Here class 2 has no member in the top grade, so its top threshold diverges
  # and the separation warning is genuine. A diverging threshold is wherever
  # EM stops: about 30 on the logit scale at the default `tol`, on either side
  # of the warning's cut-off of 30 (29.7 after the unit-free EM stop of review
  # 2026-10-26 round 2, K2-6), and 38.5 at `tol = 1e-9`, which is used here.
  expect_warning(continuous <- mixture_regression(grade ~ hours, graded,
                                                  family = "ordinal", n_classes = 2,
                                                  n_starts = 1, seed = 1,
                                                  vcov_type = "none", tol = 1e-9),
                 class = "latents_separation")
  expect_s3_class(continuous, "latents_mixture_regression")
  expect_s3_class(mixture_regression(grade ~ strategy, graded, family = "ordinal",
                                     n_classes = 2, id = "student",
                                     class_level = "group", n_starts = 1, seed = 1,
                                     vcov_type = "none"),
                  "latents_mixture_regression")
})

test_that("the classification plot names the classes, not only colours them", {
  skip_if_not_installed("ggplot2")
  fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 40),
                            n_classes = 2, id = "student", class_level = "group",
                            random = "intercept", n_starts = 1, seed = 1,
                            vcov_type = "none")
  built <- ggplot2::ggplot_build(plot(fit, what = "classification"))
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  expect_true(all(c("Class 1", "Class 2") %in% labels))
})

test_that("a growth start that cannot be built raises latents_failed_starts", {
  build <- .growth_start
  calls <- 0L
  local_mocked_bindings(.growth_start = function(spec, stats, result, first = FALSE) {
    calls <<- calls + 1L
    if (calls == 2L) stop("no starting covariance for this start")
    build(spec, stats, result, first)
  })
  expect_warning(
    fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 40),
                              n_classes = 2, id = "student", class_level = "group",
                              random = "intercept", n_starts = 2, seed = 1,
                              vcov_type = "none"),
    class = "latents_failed_starts")
  expect_identical(sum(as.data.frame(get_results(fit, "starts"))$stage == "failed"), 1L)
})
