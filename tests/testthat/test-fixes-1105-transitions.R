# Counterexamples for the transition review; fitting checks run in development.
skip_on_cran()

fix1105_transition_data <- function(occasions = 4L, persons = 30L) {
  .latents_local_seed(51116L)
  result <- data.frame(id = rep(seq_len(persons), each = occasions),
                       time = rep(seq_len(occasions), persons))
  state <- sample(0:1, nrow(result), replace = TRUE)
  result$a <- rnorm(nrow(result), 3 * state, .7)
  result$b <- rnorm(nrow(result), -3 * state, .7)
  result
}

fix1105_transition_capture <- function(expr) {
  conditions <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    conditions[[length(conditions) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, conditions = conditions)
}

test_that("short-panel latent transition ridges are refused before fitting", {
  two <- fix1105_transition_data(2L)
  expect_error(lta(two, c("a", "b"), "id", 2, "time", n_group_classes = 2,
                   n_starts = 1, seed = 1), class = "latents_unidentified")
  expect_error(lta(two, c("a", "b"), "id", 2, "time", mover_stayer = TRUE,
                   n_starts = 1, seed = 1), class = "latents_unidentified")
  three <- fix1105_transition_data(3L)
  expect_error(lta(three, c("a", "b"), "id", 2, "time", n_group_classes = 2,
                   transitions = "occasion", n_starts = 1, seed = 1),
               class = "latents_unidentified")
  expect_error(lta(three, c("a", "b"), "id", 2, "time", mover_stayer = TRUE,
                   order = 2, n_starts = 1, seed = 1), class = "latents_unidentified")
  expect_invisible(.lta_check_short_panel(3L, 2L, "homogeneous"))
  expect_invisible(.lta_check_short_panel(4L, 2L, "occasion", 2L, TRUE))
  # Covariate-dependent identification is deliberately outside this guard.
  expect_invisible(.lta_check_short_panel(2L, 2L, "occasion", has_covariates = TRUE))
})

test_that("positive but negligible transition support is distinct from an empty origin", {
  pairs <- array(c(.05, .95, .2, .8), c(1L, 2L, 2L))
  expectation <- list(moments = list(list(pairs = list(pairs), pairs2 = NULL)))
  support <- .lta_origin_support(expectation, FALSE, 2L)
  expect_equal(unname(support$first[, 1]), c(.25, 1.75))
  expect_identical(.lta_empty_origins(expectation, FALSE, 2L)$first,
                   matrix(FALSE, 2L, 1L))
  weak <- .lta_weak_origins(support, "homogeneous")
  expect_identical(weak$first, matrix(c(TRUE, FALSE), 2L, 1L))
  expect_warning(.lta_warn_weak_origins(weak$first),
                 class = "latents_weak_transition_support")
  pairs[, 1L, ] <- 0
  expectation$moments[[1]]$pairs <- list(pairs)
  expect_identical(.lta_empty_origins(expectation, FALSE, 2L)$first,
                   matrix(c(TRUE, FALSE), 2L, 1L))
  expect_false(any(.lta_weak_origins(.lta_origin_support(expectation, FALSE, 2L),
                                   "homogeneous")$first))
})

test_that("first-order support in a second-order chain reads only its first move", {
  first <- array(c(0, 1, 0, 1), c(1L, 2L, 2L))
  later <- array(c(1, 0, 1, 0), c(1L, 2L, 2L))
  second <- array(rep(.25, 8), c(1L, 4L, 2L))
  expectation <- list(moments = list(list(pairs = list(first, later), pairs2 = list(second))))
  support <- .lta_origin_support(expectation, FALSE, 2L)
  expect_equal(unname(support$first[, 1]), c(0, 2))
  expect_equal(unname(support$second[, 1]), rep(.5, 4))
  expect_identical(.lta_empty_origins(expectation, FALSE, 2L)$first,
                   matrix(c(TRUE, FALSE), 2L, 1L))
})

test_that("bootstrap extraction failures stay inside the counted replicate loop", {
  withr::local_seed(51117)
  index <- 0L
  result <- .latents_cluster_bootstrap(data.frame(id = 1:3, a = c(1, 2, 3)),
    1:3, 3L, "id", 3L,
    refit = function(data) { index <<- index + 1L; list(converged = TRUE, index = index) },
    extract = function(fit) {
      if (fit$index == 1L) stop("missing occasion signature")
      if (fit$index == 2L) return(c(1, Inf))
      c(2, 3)
    }, reference_names = c("a", "b"))
  expect_true(all(is.na(result$estimates[1:2, ])))
  expect_equal(unname(result$estimates[3, ]), c(2, 3))
  expect_match(result$messages[1], "missing occasion signature")
  expect_match(result$messages[2], "nonfinite or have a different shape")
  expect_true(is.na(result$messages[3]))
})

test_that("coefficient schemas are checked before transition bootstrap alignment", {
  data <- fix1105_transition_data()
  reference <- fix1105_transition_capture(lta(data, c("a", "b"), "id", 2,
    "time", transitions = "occasion", n_starts = 1, seed = 1, max_iter = 0))$value
  shorter <- fix1105_transition_capture(lta(subset(data, time <= 3), c("a", "b"),
    "id", 2, "time", transitions = "occasion", n_starts = 1, seed = 1,
    max_iter = 0))$value
  shorter$converged <- TRUE # controlled refit response, not a numerical fit assertion
  local_mocked_bindings(.lta_refit = function(...) shorter,
                         .lta_align_estimates = function(...) stop("alignment must not run"))
  result <- .lta_bootstrap_replicates(reference, data, 2L, 1L, 1L, 1e-8,
                                     names(.lta_bootstrap_estimates(reference)))
  expect_true(all(is.na(result$estimates)))
  expect_true(all(grepl("occasion or design term is absent", result$messages)))
  expect_false(any(grepl("alignment must not run", result$messages)))
})

test_that("evaluate-only general fits do not claim an unsuccessful optimization", {
  data <- fix1105_transition_data()
  result <- fix1105_transition_capture(lta(data, c("a", "b"), "id", 2, "time",
    transitions = "occasion", n_starts = 1, seed = 1, max_iter = 0))
  expect_false(any(vapply(result$conditions, inherits, logical(1), "latents_unconverged")))
  expect_identical(result$value$iterations, 0L)
})

test_that("trailing dots survive logit identification tables and relabelling", {
  names <- c("transition.group_class_1.profile_1->profile_2.titleZed.",
              "initial.group_class_1.profile_1.titleZed.")
  expect_identical(.lta_coefficient_term(names), rep("titleZed.", 2))
  expect_identical(tail(.lta_name_pieces(names)[[1]], 1), "")
  moved <- .lta_relabel_profiles(stats::setNames(c(.3, -.2), names), c(2L, 1L))
  expect_true(all(endsWith(names(moved), ".titleZed.")))
})

test_that("an extreme covariate row alone is not a probability boundary", {
  initial <- matrix(log(c(.5, .5)), 3, 2, byrow = TRUE)
  moves <- array(log(.5), c(3L, 2L, 2L))
  moves[1L, , ] <- log(matrix(c(1e-10, 1 - 1e-10, .5, .5), 2, byrow = TRUE))
  expectation <- list(per_class = list(list(log_initial = initial,
    log_transition = list(moves), log_transition2 = list())))
  layout <- list(slot = matrix(1:6, 3, 2), within = matrix(TRUE, 3, 2))
  expect_false(.lta_probability_boundary(expectation, layout))
  moves[, 1L, 1L] <- log(1e-10)
  expectation$per_class[[1]]$log_transition <- list(moves)
  expect_true(.lta_probability_boundary(expectation, layout))
})

test_that("second-order first-move tables label occasion two and expose support", {
  data <- fix1105_transition_data()
  fit <- fix1105_transition_capture(lta(data, c("a", "b"), "id", 2, "time",
    order = 2, n_starts = 1, seed = 1, max_iter = 0))$value
  table <- get_results(fit, "transitions")
  expect_identical(unique(table$occasion), 2L)
  expect_true(all(c("origin_expected_count", "supported") %in% names(table)))
  second <- get_results(fit, "second_order_transitions")
  expect_true(all(c("origin_expected_count", "supported") %in% names(second)))
})

test_that("both transition engines honor the explicit all-row individual BIC", {
  data <- fix1105_transition_data()
  data$a[1] <- data$b[1] <- NA_real_
  data$w <- 1
  invisible(lapply(c(FALSE, TRUE), function(drop) {
    fit <- function(weighted) fix1105_transition_capture(lta(data, c("a", "b"),
      "id", 2, "time", missing = "fiml", drop_empty = drop,
      weights = if (weighted) "w" else NULL, max_iter = 0, n_starts = 1, seed = 1))$value
    homogeneous <- fit(FALSE)
    general <- fit(TRUE)
    expected_n <- nrow(data) - as.integer(drop)
    invisible(lapply(list(homogeneous, general), function(x) {
      expect_equal((x$bic_individual + 2 * x$log_likelihood) / x$n_parameters,
                    log(expected_n), tolerance = 1e-12)
      expect_identical(x$n_observations, nrow(data))
      expect_equal(x$n_informative, nrow(data) - 1L)
    }))
  }))
})

test_that("covariance charts alone govern spherical-to-diagonal transition nesting", {
  data <- fix1105_transition_data()
  fit <- function(model) fix1105_transition_capture(lta(data, c("a", "b"),
    "id", 2, "time", model = model, n_starts = 1, seed = 1, max_iter = 0))$value
  spherical <- fit("EII")
  diagonal <- fit("EEI")
  full <- fit("EEE")
  expect_true(.lta_nested(spherical, diagonal))
  expect_true(.lta_nested(spherical, full))
  expect_false(.lta_nested(diagonal, spherical))
})

test_that("unit-free transition logit checks allow finite effects and detect separation", {
  design <- cbind(1, c(-6.5, seq(-2, 2, length.out = 100)))
  expect_false(.lta_extreme_logits(matrix(c(-1.5, 2), 2), design))
  expect_true(.lta_extreme_logits(matrix(c(-1.5, 200), 2), design))
  shifted <- cbind(1, design[, 2] * 1000 + 3000)
  expect_false(.lta_extreme_logits(matrix(c(-7.5, .002), 2), shifted))
})

test_that("support on older weighted reporting moments uses the sampling masses", {
  pairs <- array(c(.05, .1, .15, .2, .25, .3, .35, .4), c(2L, 2L, 2L))
  expectation <- list(moments = list(list(pairs = list(pairs), pairs2 = NULL)))
  support <- .lta_origin_support(expectation, FALSE, 2L, c(0, 2))
  expected <- 2 * apply(pairs[2L, , , drop = FALSE], 2L, sum)
  expect_equal(unname(support$first[, 1L]), unname(expected))
  expect_equal(unname(support$first_by_occasion[[1]][, 1L]), unname(expected))
})

test_that("chain coordinates cannot exceed the latent path distribution dimension", {
  expect_error(.lta_check_short_panel(4L, 3L, "homogeneous", order = 2L,
                                      n_profiles = 2L), class = "latents_unidentified")
  expect_invisible(.lta_check_short_panel(4L, 2L, "homogeneous", order = 2L,
                                          n_profiles = 2L))
})

test_that("the structural one-group-class probability is not a boundary", {
  data <- fix1105_transition_data()
  fit <- fix1105_transition_capture(lta(data, c("a", "b"), "id", 2, "time",
    n_group_classes = 1, n_starts = 1, seed = 1, max_iter = 0))$value
  expect_equal(unname(fit$group_probabilities), 1)
  expect_false(fit$probability_boundary)
})
