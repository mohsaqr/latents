# Wald standard errors for the ten constrained covariance structures in the
# homogeneous lta() (class `multilpa_transitions`), whose table reports
# natural probabilities, variances and covariances: the information is taken
# in the structure's chart and carried to the table by the delta method.
#
# Two oracles, with criteria fixed before either was run:
# * the independent reference in helper-lta-structure-reference.R (its own
#   forward recursion, density, charts and numerical Hessian), retained with
#   its datasets in fixtures/lta-structure-references.rds (built by
#   equivalence/lta-structure-se/make-references.R): standard errors within
#   1e-5 relative, log likelihood within 1e-8 relative;
# * the extended engine (class `multilpa_lta`, analytic chart scores, its own
#   likelihood code), fitted to the same data and structure, whose covariance
#   of its logits and chart coordinates is carried to the same natural
#   parameters by the delta method: within 1e-4 relative.

skip_on_cran()
skip_if_not_installed("numDeriv")

lta_structure_references <- readRDS(test_path("fixtures", "lta-structure-references.rds"))

lta_structure_refit <- function(case, structure) {
  suppressWarnings(lta(case$data, case$vars, "person", 2L, time = "wave",
                       n_group_classes = case$n_group_classes, model = structure,
                       missing = case$missing, n_starts = 4L, seed = 1L,
                       tol = 1e-12, max_iter = 20000L))
}

# The extended engine's natural parameters, in parameter_inference()'s order,
# as a function of its own packed coordinates.
lta_extended_natural <- function(extended, structure) {
  k <- extended$n_profiles
  d <- length(extended$continuous)
  h <- extended$n_group_classes
  function(value) {
    p <- .lta_unpack(value, extended)
    block <- p$measurement[[1L]]
    sigma <- lapply(seq_len(k), function(profile) {
      if (is.null(block$covariances)) diag(block$variances[profile, ], d) else
        matrix(block$covariances[, , profile], d, d)
    })
    rows <- if (identical(structure, "EII")) 1L else seq_len(k)
    lower <- which(lower.tri(matrix(0, d, d), diag = TRUE))
    spread <- unlist(lapply(rows, function(profile) {
      if (endsWith(structure, "I")) diag(sigma[[profile]]) else sigma[[profile]][lower]
    }))
    initial <- unlist(lapply(seq_len(h), function(class) {
      exp(.lta_log_initial(extended$designs$initial[1L, , drop = FALSE],
                           matrix(p$initial[, , class], ncol(extended$designs$initial))))
    }))
    transition <- unlist(lapply(seq_len(h), function(class) {
      moves <- exp(.lta_log_transition(
        extended$designs$transition[[1L]][1L, , drop = FALSE],
        array(p$transition[, , , class], dim(p$transition)[1:3])))
      as.vector(t(matrix(moves, k, k)))
    }))
    c(as.vector(t(sweep(block$means, 2L, extended$centers, "+"))), spread,
      p$group_probabilities, initial, transition)
  }
}

lta_extended_se <- function(case, structure) {
  extended <- suppressWarnings(.lta_fit_general(
    case$data, case$vars, "person", "wave", 2L, case$n_group_classes, "varying", 4L,
    20000L, 1e-12, 1e-6, 1L, character(), 1e-10, "observed", "homogeneous",
    character(), character(), "invariant", "likelihood", quote(lta()),
    missing = case$missing,
    covariance_model = if (endsWith(structure, "I")) "diagonal" else "full",
    structure = structure))
  natural <- lta_extended_natural(extended, structure)
  theta <- .lta_theta(extended)
  map <- numDeriv::jacobian(natural, theta)
  list(fit = extended, estimate = natural(theta),
       standard_error = sqrt(pmax(diag(map %*% vcov(extended) %*% t(map)), 0)))
}

test_that("constrained-structure standard errors of lta() match the independent reference", {
  worst <- unlist(lapply(names(lta_structure_references$cases), function(name) {
    case <- lta_structure_references$cases[[name]]
    vapply(case$results, function(reference) {
      label <- sprintf("%s %s", name, reference$structure)
      fit <- lta_structure_refit(case, reference$structure)
      expect_true(fit$converged, label = label)
      expect_identical(fit$n_parameters, reference$n_parameters, label = label)
      expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-8,
                   label = label)
      # latents' estimate is a stationary point of the independent likelihood.
      expect_lt(max(abs(reference$gradient)), 1e-3, label = label)
      ours <- parameter_inference(fit)
      expect_equal(ours$estimate, reference$estimate, tolerance = 1e-6, label = label)
      informative <- reference$standard_error > 1e-8
      relative <- abs(ours$standard_error - reference$standard_error)[informative] /
        reference$standard_error[informative]
      expect_lt(max(relative), 1e-5, label = label)
      # A degenerate coordinate (a single group class's share) has no error.
      expect_true(all(ours$standard_error[!informative] < 1e-8), label = label)
      max(relative)
    }, numeric(1))
  }))
  message("largest relative SE difference by case and structure:\n",
          paste(sprintf("%s %.2e", names(worst), worst), collapse = "\n"))
})

test_that("the retained reference is reproduced by the reference code", {
  case <- lta_structure_references$cases$one_class
  retained <- case$results$VEI
  fit <- lta_structure_refit(case, "VEI")
  live <- reference_lta_structure_se(fit, case$data,
                                     jacobian = function(f, x) numDeriv::jacobian(f, x))
  expect_equal(live$log_likelihood, retained$log_likelihood, tolerance = 1e-10)
  expect_equal(live$standard_error, retained$standard_error, tolerance = 1e-6)
})

test_that("constrained-structure standard errors of lta() agree with the extended engine", {
  cases <- lta_structure_references$cases
  checks <- c(lapply(names(cases$one_class$results),
                     function(structure) list("one_class", structure)),
              list(list("two_classes", "VEE"), list("fiml_ragged", "EVV"),
                   list("three_indicators", "EVE")))
  worst <- vapply(checks, function(check) {
    case <- cases[[check[[1L]]]]
    structure <- check[[2L]]
    label <- paste(check, collapse = " ")
    fit <- lta_structure_refit(case, structure)
    extended <- lta_extended_se(case, structure)
    expect_equal(extended$fit$log_likelihood, fit$log_likelihood, tolerance = 1e-9,
                 label = label)
    expect_identical(extended$fit$n_parameters, fit$n_parameters, label = label)
    ours <- parameter_inference(fit)
    expect_equal(ours$estimate, extended$estimate, tolerance = 1e-5, label = label)
    informative <- extended$standard_error > 1e-8
    relative <- abs(ours$standard_error - extended$standard_error)[informative] /
      extended$standard_error[informative]
    expect_lt(max(relative), 1e-4, label = label)
    max(relative)
  }, numeric(1))
  message("largest relative SE difference from the extended engine:\n",
          paste(sprintf("%s %.2e", vapply(checks, paste, character(1), collapse = " "),
                        worst), collapse = "\n"))
})
