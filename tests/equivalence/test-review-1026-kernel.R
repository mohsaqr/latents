# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1026-kernel.R,
# so that they run in CI's equivalence job rather than under R CMD check.

kernel_review_softmax <- function(scores) exp(scores) / rowSums(exp(scores))

# A d x d x k covariance array as one row per cell, for a keyed comparison.
kernel_review_cells <- function(sigma, structure) {
  cells <- expand.grid(row = seq_len(dim(sigma)[1L]), column = seq_len(dim(sigma)[2L]),
                       profile = seq_len(dim(sigma)[3L]))
  data.frame(structure = structure, cells, covariance = as.vector(sigma))
}

test_that("structure posteriors match brute force, and survive a -5e15 offset per row", {
  set.seed(5)
  n <- 30L
  k <- 3L
  n_group_classes <- 2L
  n_groups <- 6L
  group_index <- rep(seq_len(n_groups), each = 5L)
  w <- cbind(1, stats::rnorm(n))
  v <- cbind(1, stats::rnorm(n_groups))
  w_group <- cbind(1, stats::rnorm(n_groups))
  gamma <- cbind(0, matrix(stats::rnorm(2L * (k - 1L)), 2L))
  delta <- cbind(0, matrix(stats::rnorm(2L * (n_group_classes - 1L)), 2L))
  class_logits <- cbind(0, matrix(stats::rnorm((n_group_classes + 2L) * (k - 1L)),
                                  n_group_classes + 2L))
  observation <- list(nesting = "observation", w = w)
  group <- list(nesting = "group", w = w_group, group_index = group_index)
  two_level <- list(nesting = "two-level", w = w, v = v, n = n,
                    n_group_classes = n_group_classes, group_index = group_index,
                    n_groups = n_groups)
  brute <- function(relative) {
    density <- exp(relative)
    prior <- kernel_review_softmax(w %*% gamma)
    group_prior <- kernel_review_softmax(w_group %*% gamma)
    group_density <- exp(rowsum(relative, group_index))
    eta <- kernel_review_softmax(v %*% delta)
    conditional <- lapply(seq_len(n_group_classes), function(h) {
      kernel_review_softmax(latents:::.mixture_two_level_design(two_level, h) %*% class_logits)
    })
    by_class <- vapply(seq_len(n_group_classes), function(h) {
      as.vector(exp(rowsum(log(rowSums(density * conditional[[h]])), group_index)))
    }, numeric(n_groups))
    rho <- eta * by_class / rowSums(eta * by_class)
    list(tau = density * prior / rowSums(density * prior),
         log_likelihood = sum(log(rowSums(density * prior))),
         group_tau = group_density * group_prior / rowSums(group_density * group_prior),
         rho = rho,
         two_level_tau = Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
           rho[group_index, h] * density * conditional[[h]] /
             rowSums(density * conditional[[h]])
         })))
  }
  reference <- "independent: brute-force mixture posteriors by Bayes' rule, base R"
  covers <- "mixture_regression()"
  note <- "Package side: the internal class-structure E-step (.mixture_structure_expectation) of mixture_regression(); the reference uses the package's two-level design matrix (.mixture_two_level_design) as input."
  log_density <- matrix(stats::rnorm(n * k, -3), n, k)
  truth <- brute(log_density)
  first <- latents:::.mixture_structure_expectation(observation, list(gamma = gamma), log_density)
  expect_reference(first$tau, truth$tau, reference = reference, covers = covers,
                   quantity = "observation-level posteriors", tolerance = 1e-13,
                   measure = "relative", note = note)
  expect_reference(first$log_likelihood, truth$log_likelihood, reference = reference, covers = covers,
                   quantity = "observation-level log likelihood", tolerance = 1e-13,
                   measure = "relative", note = note)
  second <- latents:::.mixture_structure_expectation(group, list(gamma = gamma), log_density)
  expect_reference(second$group_tau, truth$group_tau, reference = reference, covers = covers,
                   quantity = "group-level posteriors", tolerance = 1e-13,
                   measure = "relative", note = note)
  third <- latents:::.mixture_structure_expectation(
    two_level, list(delta = delta, class_logits = class_logits), log_density)
  expect_reference(third$rho, truth$rho, reference = reference, covers = covers,
                   quantity = "two-level group-class posteriors", tolerance = 1e-13,
                   measure = "relative", note = note)
  expect_reference(third$tau, truth$two_level_tau, reference = reference, covers = covers,
                   quantity = "two-level row-class posteriors", tolerance = 1e-13,
                   measure = "relative", note = note)
  # A large per-row constant changes no posterior. Before 8ac3576 the
  # observation posteriors were off by about 0.44 here. The brute force sees
  # the densities relative to each row's maximum.
  extreme <- log_density + (-5e15 * stats::runif(n, 0.5, 1.5))
  shifted_truth <- brute(extreme - apply(extreme, 1L, max))
  shifted <- latents:::.mixture_structure_expectation(observation, list(gamma = gamma), extreme)
  expect_reference(shifted$tau, shifted_truth$tau, reference = reference, covers = covers,
                   quantity = "observation-level posteriors, -5e15 offset per row", tolerance = 1e-13,
                   measure = "relative", note = note)
  shifted_group <- latents:::.mixture_structure_expectation(group, list(gamma = gamma), extreme)
  expect_reference(shifted_group$group_tau, shifted_truth$group_tau, reference = reference, covers = covers,
                   quantity = "group-level posteriors, -5e15 offset per row", tolerance = 1e-13,
                   measure = "relative", note = note)
  shifted_two <- latents:::.mixture_structure_expectation(
    two_level, list(delta = delta, class_logits = class_logits), extreme)
  expect_reference(shifted_two$rho, shifted_truth$rho, reference = reference, covers = covers,
                   quantity = "two-level group-class posteriors, -5e15 offset per row", tolerance = 1e-13,
                   measure = "relative", note = note)
  expect_reference(shifted_two$tau, shifted_truth$two_level_tau, reference = reference, covers = covers,
                   quantity = "two-level row-class posteriors, -5e15 offset per row", tolerance = 1e-13,
                   measure = "relative", note = note)
  expect_equal(rowSums(shifted_two$tau), rep(1, n), tolerance = 1e-14)
})

test_that("the forward-backward recursions equal path enumeration", {
  set.seed(9)
  n_groups <- 4L
  k <- 3L
  n_occasions <- 4L
  lengths <- c(4L, 2L, 3L, 1L)
  layout <- list(slot = matrix(1, n_groups, n_occasions), n_occasions = n_occasions,
                 within = outer(lengths, seq_len(n_occasions), ">="))
  emission <- lapply(seq_len(n_occasions), function(t) {
    matrix(stats::rnorm(n_groups * k, -1), n_groups, k)
  })
  initial <- matrix(c(.5, .3, .2), n_groups, k, byrow = TRUE)
  initial[2L, ] <- c(.1, .1, .8)
  transition <- matrix(c(.7, .2, .1, .1, .8, .1, .3, .3, .4), k, k, byrow = TRUE)
  paths <- as.matrix(expand.grid(rep(list(seq_len(k)), n_occasions)))
  brute <- lapply(seq_len(n_groups), function(g) {
    prefixes <- unique(paths[, seq_len(lengths[g]), drop = FALSE])
    log_path <- apply(prefixes, 1L, function(state) {
      log(initial[g, state[1L]]) + emission[[1L]][g, state[1L]] +
        sum(vapply(seq_along(state)[-1L], function(t) {
          log(transition[state[t - 1L], state[t]]) + emission[[t]][g, state[t]]
        }, numeric(1)))
    })
    total <- log(sum(exp(log_path)))
    list(log_likelihood = total, posterior = vapply(seq_len(lengths[g]), function(t) {
      vapply(seq_len(k), function(state) {
        sum(exp(log_path[prefixes[, t] == state] - total))
      }, numeric(1))
    }, numeric(k)))
  })
  scaled <- latents:::.latents_forward_backward(emission, layout, log(initial), log(transition))
  logged <- latents:::.latents_forward_backward(emission, layout, log(initial),
    lapply(seq_len(n_occasions - 1L), function(t) {
      array(rep(log(transition), each = n_groups), c(n_groups, k, k))
    }))
  reference <- "independent: enumeration of every state path, base R"
  note <- "Package side: the internal forward-backward recursions (.latents_forward_backward) of lta()."
  expected <- vapply(brute, `[[`, numeric(1), "log_likelihood")
  groups <- sprintf("group %d", seq_len(n_groups))
  expect_reference(scaled$log_scaled, expected, reference = reference, covers = "lta()",
                   quantity = "sequence log likelihood, one transition matrix", tolerance = 1e-13,
                   measure = "relative", key = groups, note = note)
  expect_reference(logged$log_scaled, expected, reference = reference, covers = "lta()",
                   quantity = "sequence log likelihood, transition array per occasion", tolerance = 1e-13,
                   measure = "relative", key = groups, note = note)
  # One row per group, observed occasion and state.
  cells <- do.call(rbind, lapply(seq_len(n_groups), function(g) {
    data.frame(group = g, occasion = rep(seq_len(lengths[g]), each = k),
               state = rep(seq_len(k), lengths[g]))
  }))
  smoothed <- function(pass) {
    data.frame(cells, posterior = mapply(function(g, t, state) {
      exp(pass$alpha[[t]][g, state] + pass$beta[[t]][g, state] - expected[g])
    }, cells$group, cells$occasion, cells$state))
  }
  enumerated <- data.frame(cells, posterior = mapply(function(g, t, state) {
    brute[[g]]$posterior[state, t]
  }, cells$group, cells$occasion, cells$state))
  keys <- c("group", "occasion", "state")
  expect_reference(smoothed(scaled), enumerated, reference = reference, covers = "lta()",
                   quantity = "smoothed state posteriors, one transition matrix", tolerance = 1e-12,
                   measure = "relative", key = keys, note = note)
  expect_reference(smoothed(logged), enumerated, reference = reference, covers = "lta()",
                   quantity = "smoothed state posteriors, transition array per occasion", tolerance = 1e-12,
                   measure = "relative", key = keys, note = note)
})

test_that("the analytic FIML score is the numerical gradient of the likelihood", {
  skip_without_reference("numDeriv")
  skip_without_reference("MASS")
  set.seed(4)
  y <- rbind(MASS::mvrnorm(60L, c(-1.5, -1), matrix(c(1, .5, .5, 1), 2L)),
             MASS::mvrnorm(60L, c(1.5, 1), matrix(c(1, -.3, -.3, 1.5), 2L)))
  y <- y[sample(nrow(y)), ]
  y[sample(length(y), 30L)] <- NA
  frame <- data.frame(g = rep(seq_len(20L), each = 6L), a = y[, 1L], b = y[, 2L])
  centre <- colMeans(y, na.rm = TRUE)
  x <- sweep(y, 2L, centre, "-")
  scores <- do.call(rbind, lapply(c("varying", "equal"), function(variance_model) {
    fit <- suppressMessages(multilpa(frame, c("a", "b"), "g", 2L, 2L,
      n_starts = 2L, seed = 3L, missing = "fiml", covariance_model = "full",
      variance_model = variance_model))
    centred <- fit
    centred$means <- sweep(fit$means, 2L, centre, "-")
    theta <- latents:::.multilpa_coefficients(centred, "unconstrained") + 0.03
    log_likelihood <- function(point) {
      latents:::.multilpa_expectation(x, fit$group_index, latents:::.multilpa_decode(point, fit))$log_likelihood
    }
    data.frame(variance_model = variance_model, coordinate = seq_along(theta),
               score = unname(-latents:::.multilpa_score(theta, x, fit)),
               group_score = unname(colSums(latents:::.multilpa_group_scores(theta, x, fit))),
               numerical = numDeriv::grad(log_likelihood, theta,
                                          method.args = list(eps = 1e-4, d = 1e-3, r = 6)))
  }))
  keys <- c("variance_model", "coordinate")
  numerical <- data.frame(scores[keys], score = scores$numerical)
  # Scaled: score entries can sit near zero, where numDeriv's error is
  # absolute; a relative measure would judge its rounding noise. 1e-6, the
  # original test's tolerance.
  expect_reference(scores[c(keys, "score")], numerical, reference = "independent: numDeriv derivative of the package's own likelihood (self-consistency)",
                   covers = c("multilpa()", "vcov.multilpa()"),
                   quantity = "FIML score at a point off the optimum", tolerance = 1e-6,
                   measure = "scaled", key = keys,
                   note = "Package side: the internal analytic score (.multilpa_score) behind the standard errors of multilpa(); reference: numDeriv::grad() of the package's FIML log likelihood.")
  expect_reference(data.frame(scores[keys], score = scores$group_score), numerical, reference = "independent: numDeriv derivative of the package's own likelihood (self-consistency)",
                   covers = c("multilpa()", "vcov.multilpa()"),
                   quantity = "sum of the group FIML scores at a point off the optimum", tolerance = 1e-6,
                   measure = "scaled", key = keys,
                   note = "Package side: the internal group scores (.multilpa_group_scores) behind the robust standard errors of multilpa(); reference: numDeriv::grad() of the package's FIML log likelihood.")
})

test_that("each covariance M-step reproduces mclust's, or attains a lower objective", {
  skip_without_reference("mclust")
  skip_without_reference("MASS")
  set.seed(11)
  d <- 3L
  k <- 3L
  x <- rbind(MASS::mvrnorm(100L, c(0, 0, 0), matrix(c(1, .5, .2, .5, 2, .3, .2, .3, 1.5), 3L)),
             MASS::mvrnorm(100L, c(3, 1, -1), matrix(c(2, -.4, 0, -.4, 1, .1, 0, .1, .7), 3L)),
             MASS::mvrnorm(100L, c(-2, 3, 2), diag(c(.5, 1.5, 3))))
  z <- matrix(stats::runif(300L * k), 300L, k) +
    3 * outer(rep(seq_len(k), each = 100L), seq_len(k), "==")
  z <- z / rowSums(z)
  weights <- colSums(z)
  means <- t(crossprod(x, z)) / weights
  scatter <- lapply(seq_len(k), function(profile) {
    residuals <- sweep(x, 2L, means[profile, ], "-")
    crossprod(residuals, residuals * z[, profile])
  })
  objective <- function(sigma) sum(vapply(seq_len(k), function(profile) {
    weights[profile] * as.numeric(determinant(sigma[, , profile])$modulus) +
      sum(solve(sigma[, , profile]) * scatter[[profile]])
  }, numeric(1)))
  structures <- latents:::.multilpa_structures()
  ours <- lapply(structures, \(structure) latents:::.multilpa_structure_covariances(scatter, weights, structure, 1e-8))
  reference <- lapply(structures, \(structure) {
    do.call(mclust::mstep, list(data = x, modelName = structure, z = z),
            envir = asNamespace("mclust"))$parameters$variance$sigma
  })
  # Every structure attains mclust's objective or a lower one.
  invisible(Map(\(mine, theirs, structure) {
    expect_lte(objective(mine), objective(theirs) + 1e-9 * abs(objective(theirs)), label = structure)
  }, ours, reference, structures))
  # Closed forms, and fixed points solved to tolerance: the covariances are
  # mclust's. VEE, EVE and VVE are compared by their objective only (above).
  # 1e-8 relative, the original test's tolerance (the VEI and VEV fixed points).
  compared <- !structures %in% c("VEE", "EVE", "VVE")
  expect_reference(do.call(rbind, Map(kernel_review_cells, ours[compared], structures[compared])),
                   do.call(rbind, Map(kernel_review_cells, reference[compared], structures[compared])),
                   reference = "mclust", covers = c("lpa()", "multilpa()"),
                   quantity = "covariance M-step, every closed-form structure", tolerance = 1e-8,
                   measure = "relative", key = c("structure", "profile", "row", "column"),
                   note = "Package side: the internal covariance M-step (.multilpa_structure_covariances) of lpa()/multilpa(); reference: mclust::mstep().")
  expect_reference(vapply(structures, \(structure) latents:::.multilpa_structure_parameters(structure, k, d), integer(1)),
                   vapply(structures, \(structure) as.integer(mclust::nVarParams(structure, d, k)), integer(1)),
                   reference = "mclust", covers = c("lpa()", "multilpa()", "logLik.multilpa()"),
                   quantity = "number of covariance parameters, three profiles, three indicators",
                   tolerance = 0,
                   note = "Package side: the internal parameter count (.multilpa_structure_parameters) behind the degrees of freedom of logLik(); reference: mclust::nVarParams().")
})

test_that("a floored axis-parallel M-step stays the constrained maximizer", {
  objective <- function(variances, scatter, weights) {
    sum(weights * rowSums(log(variances))) + sum(scatter / variances)
  }
  # Profile 2 has an effective membership below one: its scatter diagonal is
  # below the floor while its variance is not. The floor does not bind at the
  # maximum, so the closed form must equal the exact unconstrained one.
  weights <- c(50, 0.5)
  scatter <- rbind(c(40, 60), c(8e-7, 1))
  floor_value <- 1e-6
  evi <- latents:::.multilpa_structure_variances(scatter, weights, "EVI", floor_value)
  log_det <- rowSums(log(scatter)) / 2
  exact <- (scatter / exp(log_det)) * sum(exp(log_det)) / sum(weights)
  expect_gt(min(exact), floor_value)
  expect_reference(evi, exact, reference = "independent: closed-form EVI maximiser, base R",
                   covers = c("lpa()", "multilpa()"),
                   quantity = "EVI variances, floor not binding", tolerance = 1e-12,
                   measure = "relative",
                   note = "Package side: the internal axis-parallel M-step (.multilpa_structure_variances) of lpa()/multilpa().")
  expect_equal(objective(evi, scatter, weights),
               objective(latents:::.multilpa_bounded_variances(scatter, weights, "EVI", floor_value),
                         scatter, weights), tolerance = 1e-13)
  vei <- latents:::.multilpa_structure_variances(scatter, weights, "VEI", floor_value)
  expect_equal(vei, latents:::.multilpa_bounded_variances(scatter, weights, "VEI", floor_value),
               tolerance = 1e-7)
  # A binding floor: the KKT solve is no worse than constrained Nelder-Mead
  # in the volume/shape chart, and the result stays inside the family.
  weights <- c(40, 30, 0.6)
  scatter <- rbind(c(40, 20, 60), c(30, 60, 3), c(1e-7, 0.3, 0.9))
  floor_value <- 1e-3
  ours <- latents:::.multilpa_structure_variances(scatter, weights, "EVI", floor_value)
  design <- latents:::.multilpa_eigenvalue_design(3L, 3L, "EVI")
  chart <- function(theta) matrix(exp(design %*% theta), 3L, 3L, byrow = TRUE)
  found <- stats::constrOptim(c(log(2), rep(0, ncol(design) - 1L)),
    function(theta) objective(chart(theta), scatter, weights), NULL,
    ui = design, ci = rep(log(floor_value), nrow(design)) - 1e-12,
    method = "Nelder-Mead", control = list(maxit = 20000L, reltol = 1e-14))
  expect_gte(min(ours), floor_value * (1 - 1e-12))
  expect_lte(objective(ours, scatter, weights), found$value + 1e-9)
  expect_lt(diff(range(rowMeans(log(ours)))), 1e-12)
})

test_that("the noise hypervolume is mclust's, in logs", {
  skip_without_reference("mclust")
  set.seed(7)
  x <- matrix(stats::rnorm(150L), 50L, 3L) %*% matrix(c(1, .3, 0, 0, 1, .5, 0, 0, 2), 3L)
  note <- "Package side: the internal log hypervolume (.multilpa_log_hypervolume) of the noise component of multilpa()."
  expect_reference(latents:::.multilpa_log_hypervolume(x), log(mclust::hypvol(x)),
                   reference = "mclust", covers = "multilpa()",
                   quantity = "log hypervolume of the data, three variables", tolerance = 1e-12,
                   measure = "relative", note = note)
  expect_reference(latents:::.multilpa_log_hypervolume(x[, 1L, drop = FALSE]),
                   log(diff(range(x[, 1L]))),
                   reference = "independent: log range of one variable, base R", covers = "multilpa()",
                   quantity = "log hypervolume of one variable", tolerance = 1e-14,
                   measure = "relative", note = note)
})
