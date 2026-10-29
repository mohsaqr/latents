# Covariances for chosen pairs of continuous indicators
# (`lpa(..., local_dependence = "y1:y2")`) against Latent GOLD 6.1, the Mplus 9
# Demo, mclust and an independent implementation written from the formula
# (features/pair-covariance/reference.R). The kit, its outputs and the design
# are described in features/pair-covariance/README.md. The comparisons read
# the retained outputs offline; latents is refitted here.
#
# Tolerances (measured values in the README).
# - Latent GOLD's `write=` file holds full precision; it converges to a
#   gradient tolerance of 1e-12, latents (pc_fit()) to a relative change of
#   1e-12. The measured log-likelihood gap is 1.7e-9; the tolerance is 1e-6.
# - Mplus's RESULTS file prints eight significant digits (0.30437437E+04),
#   so its log likelihood and BIC are known to half a unit of the eighth
#   digit, 5e-5 at |LL| ~ 3000 and 5e-4 above 10000; the tolerance is that
#   plus 2e-6.
# - A log-likelihood gap delta moves a parameter by about sqrt(2 delta)
#   standard errors in the worst direction; estimates are compared on that
#   scale, |difference| <= 2e-3 standard errors (measured: 6.0e-5).
# - Standard errors: relative 2e-3 (measured: 1.0e-5), the curvature at points
#   apart by the optimizer gap and each program's numerical derivatives.
#   With missing values Latent GOLD's standard errors differ from the
#   observed-data Hessian by up to 6% (a convention, README); there latents
#   is compared with the reference's Hessian instead.
# - Posteriors: Latent GOLD prints six decimals (5e-7), Mplus twelve; they
#   move with the optimizer gap: 2e-5 (measured: 1.1e-5; 2.0e-5 when latents
#   stopped at 1e-10).

sys.source(file.path("features", "pair-covariance", "helpers.R"), envir = environment())
sys.source(file.path("features", "pair-covariance", "reference.R"), envir = environment())

pc_tolerance <- list(log_likelihood = 1e-6, scaled_estimate = 2e-3, standard_error = 2e-3,
                     posterior = 2e-5)

#' Half a unit of the eighth significant digit Mplus writes for `value`
pc_mplus_resolution <- function(value) 0.5e-8 * 10^(floor(log10(abs(value))) + 1)

pc_all_cases <- pc_cases()

#' Fit every case (and every reference-study dataset) on two R workers: a
#' socket cluster whose workers load this source tree.
pc_parallel <- function(items, worker, ...) {
  cluster <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  setup <- function(root, helpers, reference) {
    suppressMessages(pkgload::load_all(root, quiet = TRUE))
    sys.source(helpers, envir = globalenv())
    sys.source(reference, envir = globalenv())
    NULL
  }
  environment(setup) <- globalenv()
  environment(worker) <- globalenv()
  parallel::clusterCall(cluster, setup, normalizePath(".."),
                        normalizePath(pc_root("helpers.R")),
                        normalizePath(pc_root("reference.R")))
  parallel::clusterCall(cluster, setwd, normalizePath("."))
  parallel::clusterApplyLB(cluster, items, worker, ...)
}

pc_fits <- stats::setNames(pc_parallel(pc_all_cases, function(case) pc_fit(case)),
                           names(pc_all_cases))

# Every reference run: the case, the program, and the output's name. A floor
# case has two Mplus runs: the floor as an inequality (MODEL CONSTRAINT) and
# its companion `_fixed`, the binding bound held as an equality.
pc_runs <- do.call(rbind, lapply(pc_all_cases, function(case) {
  rbind(
    data.frame(case = case$name, program = case$spec$programs, output = case$name),
    if (!is.null(case$spec$floor))
      data.frame(case = case$name, program = "mplus_fixed",
                 output = paste0(case$name, "_fixed")))
}))

pc_available <- function(output, program) {
  file.exists(pc_root("returned", paste0(output, if (program == "latentgold") ".lst" else ".out")))
}

#' One reference against latents: every compared quantity, as one row each
pc_compare <- function(run) {
  case <- pc_all_cases[[run$case]]
  fit <- pc_fits[[run$case]]
  reference <- if (run$program == "latentgold") pc_latentgold(run$output, case) else
    pc_mplus(run$output, case)
  aligned <- ref_relabel(pc_posteriors(fit, reference$id), reference$posteriors)
  floor <- !is.null(case$spec$floor)
  joined <- pc_join(pc_latents_table(fit, inference = !floor), reference$table, aligned$order)
  resolution <- if (run$program == "latentgold") 0 else
    pc_mplus_resolution(reference$log_likelihood)
  list(case = case, fit = fit, reference = reference, aligned = aligned, joined = joined,
       floor = floor, gap = fit$log_likelihood - reference$log_likelihood,
       ll_tolerance = pc_tolerance$log_likelihood + resolution,
       bic_tolerance = 2 * pc_tolerance$log_likelihood +
         if (run$program == "latentgold") 5e-5 else pc_mplus_resolution(reference$bic))
}

pc_rows <- list()

invisible(lapply(seq_len(nrow(pc_runs)), function(index) {
  run <- pc_runs[index, ]
  test_that(sprintf("%s: latents agrees with %s", run$output, run$program), {
    expect_true(pc_available(run$output, run$program), info = "no retained output")
    skip_if_not(pc_available(run$output, run$program))
    result <- pc_compare(run)
    reference <- result$reference
    fit <- result$fit
    joined <- result$joined
    if (run$program != "latentgold") expect_true(reference$terminated, info = run$output)
    expect_identical(reference$id, seq_len(nrow(result$case$data)))
    # latents must reach the best maximum the program found. Where the program
    # stopped lower, that is recorded and its estimates are not compared
    # (they describe another solution).
    expect_gte(result$gap, -result$ll_tolerance)
    lower <- result$gap > result$ll_tolerance
    missing_values <- isTRUE(result$case$spec$missing)
    compare_se <- !lower && !result$floor && !(missing_values && run$program == "latentgold")
    # The tight scale of an estimate is its standard error: the reference's,
    # else latents' (Mplus prints none for the proportions), else (a floor
    # fit, which has none) 1e-3 of the estimate.
    scale <- ifelse(!is.na(joined$standard_error_reference), joined$standard_error_reference,
                    ifelse(!is.na(joined$standard_error_latents), joined$standard_error_latents,
                           1e-3 * pmax(1, abs(joined$estimate_reference))))
    compared <- !is.na(joined$estimate_reference)
    pc_rows[[length(pc_rows) + 1L]] <<- data.frame(
      case = run$case, program = run$program, log_likelihood_gap = result$gap,
      ll_tolerance = result$ll_tolerance, reference_lower = lower,
      posterior_gap = result$aligned$gap, n_parameters = fit$n_parameters,
      reference_n_parameters = reference$n_parameters,
      bic_difference = fit$bic - reference$bic,
      max_scaled_estimate = if (lower) NA_real_ else
        max(abs(joined$estimate_latents - joined$estimate_reference)[compared] / scale[compared]),
      max_relative_se = if (!compare_se) NA_real_ else max(abs(
        joined$standard_error_latents / joined$standard_error_reference - 1), na.rm = TRUE),
      max_relative_se_latentgold_missing = if (missing_values && run$program == "latentgold")
        max(abs(joined$standard_error_latents / joined$standard_error_reference - 1)) else
          NA_real_)
    # The companion holds one bound as an equality, which latents counts as a
    # parameter at its bound.
    expect_equal(fit$n_parameters,
                 reference$n_parameters + if (run$program == "mplus_fixed") 1 else 0)
    if (lower) return(invisible())
    expect_lte(abs(result$gap), result$ll_tolerance)
    if (run$program != "mplus_fixed") {
      expect_lte(abs(fit$bic - reference$bic), result$bic_tolerance)
    }
    expect_lte(result$aligned$gap, pc_tolerance$posterior)
    expect_true(all(abs(joined$estimate_latents - joined$estimate_reference)[compared] <=
                      pc_tolerance$scaled_estimate * scale[compared]),
                info = paste(run$output, run$program))
    if (compare_se) {
      with_se <- !is.na(joined$standard_error_reference)
      expect_true(all(abs(joined$standard_error_latents[with_se] /
                            joined$standard_error_reference[with_se] - 1) <=
                        pc_tolerance$standard_error), info = paste(run$output, run$program))
    }
  })
}))

test_that("Mplus refuses the indicator in units of 1e-4 as having no variance", {
  # A limit of Mplus, not a disagreement: the scales cases are compared with
  # Latent GOLD only (README).
  invisible(lapply(c("scales_01", "scales_02"), function(name) {
    lines <- pc_read_lines(pc_root("returned", paste0(name, ".out")))
    expect_true(any(grepl("One or more variables have a variance of zero", lines)))
    expect_false(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
  }))
})

test_that("with missing values latents' standard errors are the observed-data Hessian's", {
  # Latent GOLD's differ (README); the reference differentiates the FIML log
  # likelihood itself.
  invisible(lapply(Filter(\(case) isTRUE(case$spec$missing), pc_all_cases), function(case) {
    fit <- pc_fits[[case$name]]
    x <- as.matrix(case$data[case$spec$vars])
    par <- list(proportions = as.vector(fit$profile_probabilities), means = unname(fit$means),
                covariances = unname(fit$covariances))
    errors <- ref_standard_errors(par, x, case$spec$vars, case$spec$pairs, FALSE)
    table <- pc_latents_table(fit)
    raw <- with(table, ifelse(quantity == "proportion", sprintf("probability:%d", class),
                       ifelse(quantity == "covariance",
                              sprintf("covariance:%d:%s", class, term),
                              sprintf("%s:%d:%s", quantity, class, term))))
    expect_true(all(raw %in% names(errors$standard_error)))
    expect_lte(max(abs(table$standard_error / errors$standard_error[raw] - 1)), 1e-4)
  }))
})

test_that("the kit covers every configuration against both programs", {
  table <- do.call(rbind, pc_rows)
  expect_gte(length(pc_all_cases), 60L)
  expect_gte(sum(table$program == "latentgold"), 55L)
  expect_gte(sum(table$program == "mplus"), 60L)
  expect_identical(sum(table$program == "mplus_fixed"), 4L)
  dir.create(pc_root("results"), showWarnings = FALSE)
  utils::write.csv(table, pc_root("results", "external.csv"), row.names = FALSE)
})


# -- the floor: a variance or a block's eigenvalue held at min_variance ---------------

test_that("floor fits satisfy the KKT conditions and no feasible point is higher", {
  floor_cases <- Filter(\(case) !is.null(case$spec$floor), pc_all_cases)
  expect_length(floor_cases, 4L)
  invisible(lapply(floor_cases, function(case) {
    fit <- pc_fits[[case$name]]
    floor <- case$spec$min_variance
    x <- as.matrix(case$data[case$spec$vars])
    par <- list(proportions = as.vector(fit$profile_probabilities), means = unname(fit$means),
                covariances = unname(fit$covariances))
    expect_true(fit$boundary, info = case$name)
    expect_error(parameter_inference(fit), class = "latents_boundary_fit")
    # The likelihood the reference computes at latents' estimate is latents' own.
    expect_lte(abs(ref_loglik(par, x) - fit$log_likelihood), 1e-9 * abs(fit$log_likelihood))
    blocks <- ref_blocks(case$spec$vars, case$spec$pairs)
    kkt <- do.call(rbind, lapply(seq_len(case$spec$k), function(j) {
      do.call(rbind, lapply(blocks, function(block) {
        cbind(class = j, size = length(block), ref_block_kkt(par, x, block, j, floor))
      }))
    }))
    lowest <- pmin(kkt$lambda_1, kkt$lambda_2, na.rm = TRUE)
    # Every eigenvalue is feasible, and the floor binds somewhere.
    expect_true(all(lowest >= floor * (1 - 1e-10)), info = case$name)
    at_floor <- abs(lowest / floor - 1) <= 1e-8
    expect_true(any(at_floor), info = case$name)
    # Stationary in every free direction: the gradient is zero there (scaled
    # by the size of a step of one standard deviation of the coordinate).
    n <- nrow(x)
    free_1 <- ifelse(kkt$size == 1L, !at_floor, abs(kkt$lambda_1 / floor - 1) > 1e-8)
    free_2 <- kkt$size == 2L & abs(kkt$lambda_2 / floor - 1) > 1e-8
    expect_true(all(abs(kkt$gradient_1[free_1] * kkt$lambda_1[free_1]) <= 1e-3),
                info = case$name)
    expect_true(all(abs(kkt$gradient_2[free_2] * kkt$lambda_2[free_2]) <= 1e-3),
                info = case$name)
    expect_true(all(abs(kkt$gradient_angle[kkt$size == 2L]) <= 1e-3), info = case$name)
    # At the floor the likelihood still wants a smaller eigenvalue: the
    # multiplier -d loglik / d lambda is nonnegative.
    floor_gradient <- ifelse(kkt$size == 1L | abs(kkt$lambda_2 / floor - 1) <= 1e-8,
                             ifelse(kkt$size == 1L, kkt$gradient_1, kkt$gradient_2), NA)
    expect_true(all(floor_gradient[at_floor] < 0), info = case$name)
    # Constrained optimality: an independent maximizer in a chart that cannot
    # cross the floor finds nothing higher.
    reference <- ref_fit(x, case$spec$vars, case$spec$pairs, case$spec$k,
                         equal = identical(case$spec$variance, "equal"), floor = floor,
                         n_random = 6L, seed = case$seed)
    expect_lte(reference$log_likelihood - fit$log_likelihood, 1e-6)
    expect_lte(abs(reference$log_likelihood - fit$log_likelihood), 1e-5)
  }))
})

# -- two indicators: the pair is the full covariance ------------------------------------

test_that("with two indicators the pair is VVV or EEE, as in mclust", {
  skip_if_not_installed("mclust")
  # Mclust() calls mclustBIC() unqualified (as in test-mclust.R).
  mclustBIC <- mclust::mclustBIC
  two <- Filter(\(case) isTRUE(case$spec$compare_full), pc_all_cases)
  expect_length(two, 4L)
  invisible(lapply(two, function(case) {
    fit <- pc_fits[[case$name]]
    code <- if (identical(case$spec$variance, "equal")) "EEE" else "VVV"
    full <- lpa(case$data, case$spec$vars, case$spec$k, model = code, n_starts = 20,
                seed = 1, tol = 1e-10, max_iter = 5000)
    expect_lte(abs(full$log_likelihood - fit$log_likelihood), 1e-7)
    expect_identical(full$n_parameters, fit$n_parameters)
    expect_lte(max(abs(full$subject_posteriors[, ref_relabel(
      full$subject_posteriors, fit$subject_posteriors)$order] - fit$subject_posteriors)), 1e-5)
    x <- as.matrix(case$data[case$spec$vars])
    clustered <- mclust::Mclust(x, G = case$spec$k, modelNames = code, verbose = FALSE)
    expect_lte(clustered$loglik - fit$log_likelihood, 1e-6)
    expect_identical(as.integer(clustered$df), fit$n_parameters)
    # mclust's EM from latents' own posteriors stays at latents' estimate.
    em <- getExportedValue("mclust", paste0("me", code))(
      data = x, z = fit$subject_posteriors, control = mclust::emControl(tol = c(1e-12, 1e-12)))
    expect_lte(abs(em$loglik - fit$log_likelihood), 1e-7)
    expect_lte(max(abs(t(em$parameters$mean) - unname(fit$means))), 1e-5)
    expect_lte(max(abs(em$parameters$variance$sigma - unname(fit$covariances))), 1e-5)
  }))
})

# -- the reference study: many cheap datasets against the formula ---------------------

pc_study_one <- function(index) {
  case <- ref_study_case(index)
  x <- as.matrix(case$data)
  fit <- lpa(case$data, case$vars, case$k, local_dependence = case$pairs,
             variance_model = if (case$equal) "equal" else "varying", n_starts = 10,
             seed = 1, tol = 1e-10, max_iter = 5000,
             missing = if (case$missing) "fiml" else "error")
  par <- list(proportions = as.vector(fit$profile_probabilities), means = unname(fit$means),
              covariances = unname(fit$covariances))
  pattern <- ref_pattern(case$vars, case$pairs)
  sigmas <- lapply(seq_len(case$k), \(j) matrix(fit$covariances[, , j], length(case$vars)))
  reference <- ref_fit(x, case$vars, case$pairs, case$k, equal = case$equal,
                       floor = fit$min_variance, seed = case$seed)
  row <- data.frame(
    index = index, n = case$n, k = case$k, d = length(case$vars),
    n_pairs = length(case$pairs), equal = case$equal, missing = case$missing,
    log_likelihood = fit$log_likelihood, reference_at_latents = ref_loglik(par, x),
    reference_maximum = reference$log_likelihood,
    off_pattern_zero = all(vapply(sigmas, \(s) all(s[!pattern] == 0), logical(1))),
    symmetric = all(vapply(sigmas, \(s) identical(s, t(s)), logical(1))),
    max_asymmetry = max(vapply(sigmas, \(s) max(abs(s - t(s))) / max(abs(s)), numeric(1))),
    min_eigenvalue = min(vapply(sigmas, \(s) min(eigen(s, symmetric = TRUE,
                                                        only.values = TRUE)$values), numeric(1))),
    min_variance = fit$min_variance, boundary = isTRUE(fit$boundary),
    equal_shared = !case$equal ||
      all(vapply(sigmas, \(s) identical(s, sigmas[[1L]]), logical(1))),
    n_parameters = fit$n_parameters,
    n_parameters_formula = (case$k - 1L) + case$k * length(case$vars) +
      (if (case$equal) 1L else case$k) * (length(case$vars) + length(case$pairs)),
    decrement = NA_real_, max_relative_se = NA_real_, max_relative_se_richardson = NA_real_,
    max_scaled_gradient = NA_real_)
  if (row$boundary) return(row)
  inference <- parameter_inference(fit)
  errors <- ref_standard_errors(par, x, case$vars, case$pairs, case$equal, richardson = TRUE)
  owner <- \(j) if (case$equal) "shared" else sprintf("profile_%s", j)
  # Each raw coordinate's row in parameter_inference().
  lookup <- function(name) {
    parts <- strsplit(name, ":", fixed = TRUE)[[1L]]
    hit <- switch(parts[1L],
      probability = inference$outcome == sprintf("profile_%s", parts[2L]) &
        inference$parameter == "probability" & inference$term %in% "group_class_1",
      mean = inference$outcome == sprintf("profile_%s", parts[2L]) &
        inference$parameter == "mean" & inference$term %in% parts[3L],
      variance = inference$outcome == owner(if (case$equal) 1L else parts[2L]) &
        inference$parameter == "covariance" &
        inference$term %in% sprintf("%s:%s", parts[3L], parts[3L]),
      covariance = inference$outcome == owner(if (case$equal) 1L else parts[2L]) &
        inference$parameter == "covariance" &
        inference$term %in% sprintf("%s:%s", parts[4L], parts[3L]))
    stopifnot("one inference row per coordinate" = sum(hit) == 1L)
    inference$standard_error[hit]
  }
  latents_se <- vapply(names(errors$standard_error), lookup, numeric(1))
  row$max_relative_se <- max(abs(latents_se / errors$standard_error - 1))
  richardson <- c(errors$standard_error[1L], errors$richardson_se)
  row$max_relative_se_richardson <- max(abs(latents_se / richardson - 1))
  row$decrement <- errors$decrement
  row$max_scaled_gradient <- max(abs(errors$gradient * errors$standard_error[-1L]))
  row
}

pc_study <- do.call(rbind, pc_parallel(seq_len(120L), pc_study_one))

test_that("the reference likelihood at latents' estimates is latents' likelihood", {
  expect_identical(nrow(pc_study), 120L)
  relative <- abs(pc_study$reference_at_latents - pc_study$log_likelihood) /
    abs(pc_study$log_likelihood)
  expect_true(all(relative <= 1e-9), info = paste(format(max(relative)), collapse = " "))
  expect_identical(pc_study$n_parameters, pc_study$n_parameters_formula)
})

test_that("latents reaches the independent maximum", {
  # Both maximize the same likelihood from their own starts: latents must not
  # stop below the reference by more than the convergence gap (EM to 1e-10
  # relative change). A reference maximum below latents' is the reference's
  # own local maximum.
  gap <- pc_study$reference_maximum - pc_study$log_likelihood
  expect_true(all(gap <= 1e-6), info = paste(pc_study$index[gap > 1e-6], collapse = " "))
  expect_gte(mean(abs(gap) <= 1e-6), 0.9)
})

test_that("the covariance pattern holds exactly", {
  expect_true(all(pc_study$off_pattern_zero))
  expect_true(all(pc_study$max_asymmetry <= 4 * .Machine$double.eps))
  # Exactly symmetric: the triangles of the M-step's crossprod() used to differ
  # in the last bit in about half the fits (FINDINGS.md, finding 1, fixed).
  expect_true(all(pc_study$symmetric),
              info = sprintf("%d of %d fits not exactly symmetric",
                             sum(!pc_study$symmetric), nrow(pc_study)))
  expect_true(all(pc_study$equal_shared))
  expect_true(all(pc_study$min_eigenvalue >= pc_study$min_variance))
})

test_that("standard errors and the score agree with numerical derivatives", {
  interior <- subset(pc_study, !boundary)
  expect_gte(nrow(interior), 100L)
  # The reference's central-difference Hessian (step 1e-4) is accurate to
  # about 1e-6 relative; latents differentiates its analytic score centrally.
  # 1e-4 covers both. (numDeriv's Richardson extrapolation, recorded beside
  # it, is off by up to 13% on a perfectly separated class: README.)
  expect_true(all(interior$max_relative_se <= 1e-4),
              info = paste(interior$index[interior$max_relative_se > 1e-4], collapse = " "))
  # A zero score: the Newton decrement g' H^-1 g (twice the log-likelihood
  # still available) is at the scale of the EM stopping rule.
  expect_true(all(interior$decrement <= 1e-6),
              info = paste(interior$index[interior$decrement > 1e-6], collapse = " "))
  dir.create(pc_root("results"), showWarnings = FALSE)
  utils::write.csv(pc_study, pc_root("results", "reference-study.csv"), row.names = FALSE)
})

# -- the active-set floor solver (a40d5c7) against an independent optimum ------------

test_that("the constrained diagonal floor update is the constrained optimum (KKT)", {
  skip_if_not_installed("nnls")
  problems <- lapply(seq_len(200L), ref_floor_problem)
  rows <- do.call(rbind, lapply(problems, function(problem) {
    solved <- .multilpa_bounded_variances(problem$scatter, problem$weights,
                                          problem$structure, problem$floor)
    reference <- ref_floor_solve(problem$scatter, problem$weights, problem$structure,
                                 problem$floor)
    kkt <- ref_floor_kkt(solved, problem$scatter, problem$weights, problem$structure,
                         problem$floor)
    data.frame(structure = problem$structure,
               gap = ref_floor_objective(solved, problem$scatter, problem$weights) -
                 ref_floor_objective(reference$sigma, problem$scatter, problem$weights),
               feasible = kkt$feasible, stationarity = kkt$stationarity,
               min_multiplier = kkt$min_multiplier, n_active = kkt$n_active,
               family = kkt$family)
  }))
  # The floor binds in most problems, so the active-set path is exercised.
  expect_gte(sum(rows$n_active > 0L), 150L)
  # Feasible, on the family, stationary with nonnegative multipliers.
  expect_true(all(rows$feasible >= -1e-12))
  expect_true(all(rows$family <= 1e-10))
  expect_true(all(rows$stationarity <= 1e-9), info = format(max(rows$stationarity)))
  # nnls returns nonnegative multipliers by construction; stationarity above
  # is what certifies that a nonnegative set exists.
  expect_true(all(rows$min_multiplier >= 0))
  # No feasible point of the independent barrier solve is lower; the barrier
  # stops a little inside the floor, so it may be higher by its own gap.
  expect_true(all(rows$gap <= 1e-10), info = format(max(rows$gap)))
  dir.create(pc_root("results"), showWarnings = FALSE)
  utils::write.csv(rows, pc_root("results", "floor-solver.csv"), row.names = FALSE)
})
