# Wald inference for the ten constrained covariance structures in the
# homogeneous lta() (`multilpa_transitions`). The contracts are checked here:
# the analytic chart scores are the gradient of the package's own transition
# likelihood, the coordinates count the structure's parameters, the natural
# table carries the delta method, and an unknown structure is refused by
# class. The comparisons with an independent likelihood, the extended engine
# and depmixS4 are in tests/equivalence/test-lta-structure-inference.R and
# equivalence/test-lta-structure-depmixs4.R.

skip_on_cran()

lta_structure_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(31)
  n <- 90L
  waves <- 4L
  frame <- do.call(rbind, lapply(seq_len(n), function(person) {
    path <- Reduce(function(previous, wave) {
      if (stats::runif(1L) < 0.8) previous else 3L - previous
    }, seq_len(waves - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
    data.frame(person = person, wave = seq_len(waves), state = path)
  }))
  s <- frame$state
  frame$a <- stats::rnorm(nrow(frame), c(0, 2.2)[s], c(1, 1.3)[s])
  frame$b <- stats::rnorm(nrow(frame), c(0, 1.8)[s], c(0.8, 1.1)[s]) + 0.3 * frame$a
  frame$q <- ifelse(stats::runif(nrow(frame)) < c(0.75, 0.3)[s], "yes", "no")
  frame$state <- NULL
  frame
})

lta_structure_fit <- function(structure, categorical = character()) {
  quietly(lta(lta_structure_panel, c("a", "b", categorical), "person", 2L,
              time = "wave", model = structure, categorical = categorical,
              n_starts = 2L, seed = 1L, tol = 1e-10, max_iter = 20000L))
}

# The package's likelihood and analytic score on the inference coordinates.
lta_structure_pieces <- function(fit) {
  x <- matrix(as.numeric(fit$indicator_data), nrow(fit$indicator_data))
  centers <- colMeans(x, na.rm = TRUE)
  centred <- sweep(x, 2L, centers, "-")
  layout <- .multilpa_sequence_layout(fit$group_index, fit$time_values,
                                      fit$n_groups, fit$occasions)
  view <- .multilpa_transition_view(fit, centers)
  blocks <- .multilpa_transition_shared_blocks(view, "unconstrained")
  widths <- c(measurement = length(blocks$measurement), group = length(blocks$group))
  list(theta = unname(.multilpa_transition_encode(fit, "unconstrained", centers)),
       loglik = function(value) {
         .multilpa_transition_expectation(
           centred, layout, .multilpa_transition_decode(
             value, view, widths[["measurement"]], widths[["group"]]),
           fit$categorical_data)$log_likelihood
       },
       score = function(value) {
         colSums(.multilpa_transition_group_scores(value, centred, layout, view, widths,
                                                   fit$categorical_data))
       })
}

lta_structure_check_score <- function(fit, label) {
  pieces <- lta_structure_pieces(fit)
  expect_identical(length(pieces$theta), as.integer(fit$n_parameters), label = label)
  expect_equal(pieces$loglik(pieces$theta), fit$log_likelihood, tolerance = 1e-10,
               label = label)
  set.seed(9)
  away <- pieces$theta + 0.05 * stats::rnorm(length(pieces$theta))
  numerical <- vapply(seq_along(away), function(index) {
    step <- replace(numeric(length(away)), index, 1e-5)
    (pieces$loglik(away + step) - pieces$loglik(away - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(pieces$score(away) - numerical)), 1e-4, label = label)
}

test_that("chart scores of the homogeneous lta() are the likelihood's gradient", {
  invisible(lapply(c("EII", "VII", "EVI", "VEI", "EEV", "VEV", "EVE", "VVE", "VEE",
                     "EVV"), function(structure) {
    fit <- lta_structure_fit(structure)
    lta_structure_check_score(fit, structure)
    inference <- parameter_inference(fit)
    expect_identical(.multilpa_parameter_names(inference), names(coef(fit)),
                     label = structure)
    expect_identical(rownames(vcov(fit)), names(coef(fit)), label = structure)
    free <- vcov(fit, scale = "unconstrained")
    expect_gt(min(eigen(free, symmetric = TRUE, only.values = TRUE)$values), 0,
              label = structure)
    reported <- inference$parameter != "probability"
    expect_true(all(is.finite(inference$standard_error)), label = structure)
    expect_true(all(inference$standard_error[reported] > 0), label = structure)
  }))
})

test_that("a constrained structure with a categorical indicator scores and infers", {
  # The chart sits between the means and the response logits; its width must
  # not swallow the categorical block.
  invisible(lapply(c("VEI", "EVE"), function(structure) {
    fit <- lta_structure_fit(structure, categorical = "q")
    lta_structure_check_score(fit, paste(structure, "mixed"))
    inference <- parameter_inference(fit)
    expect_true(any(inference$parameter == "response"))
    expect_true(all(is.finite(inference$standard_error)))
  }))
})

test_that("the natural spread table of a constrained structure is the delta method", {
  fit <- lta_structure_fit("VEE")
  covariance <- vcov(fit)
  free <- vcov(fit, scale = "unconstrained")
  centers <- colMeans(fit$indicator_data)
  view <- .multilpa_transition_view(fit, centers)
  jacobian <- .multilpa_transition_jacobian(fit, view)
  expect_equal(unname(covariance), unname(jacobian %*% free %*% t(jacobian)),
               tolerance = 1e-10)
  # The robust and OPG covariances exist too (90 groups, 11 parameters).
  expect_true(all(is.finite(parameter_inference(fit, vcov_type = "robust")$standard_error)))
  expect_true(all(is.finite(parameter_inference(fit, vcov_type = "opg")$standard_error)))
})

test_that("EII reports its one shared spread block once", {
  fit <- lta_structure_fit("EII")
  inference <- parameter_inference(fit)
  spread <- inference[inference$parameter == "variance", ]
  expect_identical(unique(spread$outcome), "shared")
  expect_identical(nrow(spread), 2L)
  expect_equal(spread$estimate[[1L]], spread$estimate[[2L]])
})

test_that("a covariance structure outside the fourteen is refused by class", {
  fit <- lta_structure_fit("VEI")
  fit$covariance_structure <- "XYZ"
  expect_error(parameter_inference(fit), class = "latents_unsupported_inference")
  expect_error(vcov(fit), class = "latents_unsupported_inference")
})
