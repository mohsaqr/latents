# The ML step-three distal outcome model (three_step(method = "ml")) against
# Latent GOLD 6.1 (`step3 modal ml` on latents' own posteriors, and end to
# end from its own step one), the Mplus 9 Demo (DU3STEP and DE3STEP) and an
# independent implementation (features/ml-three-step/reference.R), on the
# simulab kit in features/ml-three-step/kit/ and a sweep of 120 more
# datasets. The kit, its runs and the conventions are described in
# README.md there; the two disagreements, both fixed, are reproduced by
# repro-local-maximum.R (finding 1) and repro-zero-se.R (finding 2) and are
# described in FINDINGS.md.

ml3_dir <- if (dir.exists(file.path("features", "ml-three-step"))) {
  file.path("features", "ml-three-step")
} else file.path("equivalence", "features", "ml-three-step")
sys.source(file.path(ml3_dir, "helpers.R"), envir = environment())
sys.source(file.path(ml3_dir, "reference.R"), envir = environment())

skip_if_not(file.exists(file.path(ml3_kit(), "cases.rds")), "the ML three-step kit is not built")
skip_if_not_installed("numDeriv")

ml3_cases <- readRDS(file.path(ml3_kit(), "cases.rds"))
ml3_kit_cases <- Filter(\(case) !identical(case$config, "sweep"), ml3_cases)
ml3_sweep_cases <- Filter(\(case) identical(case$config, "sweep"), ml3_cases)

# Tolerances. Latent GOLD runs with tolerance 1e-12, and latents and the
# reference polish to a vanishing score, so on identical posteriors the
# estimates agree far inside 1e-5 standard errors and relative 1e-5.
ml3_tol <- list(estimate_se = 1e-5, relative = 1e-5, log_likelihood = 1e-6,
                reference_estimate_se = 1e-6, reference_relative = 1e-6)

# Every comparison is recorded (maxima.csv, comparisons.csv beside the kit).
ml3_log <- new.env()
ml3_log$rows <- list()
ml3_check <- function(case, variance, against, quantity, difference, tolerance) {
  difference <- max(abs(difference))
  ml3_log$rows[[length(ml3_log$rows) + 1L]] <- data.frame(
    case = case, variance = variance, against = against, quantity = quantity,
    difference = difference, tolerance = tolerance)
  invisible(NULL)
}

# One expectation per compared quantity of a test block, naming every case
# outside its tolerance, so a finding on many cases is one failure (testthat
# stops a run after ten) and the rest of the file still runs.
ml3_log$asserted <- 0L
ml3_assert <- function() {
  index <- seq_along(ml3_log$rows)
  fresh <- do.call(rbind, ml3_log$rows[index > ml3_log$asserted])
  ml3_log$asserted <- length(ml3_log$rows)
  if (is.null(fresh)) return(invisible(NULL))
  fresh$ok <- is.finite(fresh$difference) & fresh$difference <= fresh$tolerance
  keys <- unique(fresh[c("against", "quantity")])
  invisible(lapply(seq_len(nrow(keys)), \(i) {
    rows <- fresh[fresh$against == keys$against[i] & fresh$quantity == keys$quantity[i], ]
    bad <- rows[!rows$ok, ]
    expect_true(nrow(bad) == 0L, info = sprintf(
      "vs %s, %s: %d of %d outside tolerance: %s", keys$against[i], keys$quantity[i],
      nrow(bad), nrow(rows),
      paste(sprintf("%s %s %.3g (tol %.3g)", bad$case, bad$variance, bad$difference,
                    bad$tolerance), collapse = "; ")))
  }))
}

ml3_fit_all <- function(cases) {
  fits <- parallel::mclapply(cases, \(case) suppressWarnings(ml3_fit(case)), mc.cores = 2L,
                             mc.preschedule = FALSE)
  lost <- vapply(fits, \(fit) !inherits(fit, "multilpa"), logical(1))
  if (any(lost)) stop(sprintf("fits lost: %s", paste(names(fits)[lost], collapse = ", ")),
                      call. = FALSE)
  fits
}

#' three_step(method = "ml") with its warnings and error kept
ml3_run <- function(fit, case, ...) {
  warnings <- character()
  result <- withCallingHandlers(
    tryCatch(three_step(fit, case$data, case$outcome, level = case$level, method = "ml", ...),
             error = \(e) e),
    warning = \(w) {
      warnings <<- c(warnings, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  list(result = if (inherits(result, "error")) NULL else as.data.frame(result),
       error = if (inherits(result, "error")) result else NULL, warnings = unique(warnings))
}

ml3_inner <- function(fit, case, variance) {
  inputs <- ml3_level_inputs(fit, case)
  suppressWarnings(.multilpa_step_ml_fit(.multilpa_level_assignments(fit, case$level),
                                         inputs$values, variance))
}

ml3_fits <- ml3_fit_all(ml3_cases)
ml3_variances <- c("varying", "equal")
ml3_lg_cases <- Filter(\(case) isTRUE(case$lg), ml3_kit_cases)

test_that("the kit's posteriors are latents' own fit, refitted", {
  invisible(lapply(ml3_cases, \(case) {
    inputs <- ml3_level_inputs(ml3_fits[[case$name]], case)
    expect_lt(max(abs(inputs$posteriors - case$posteriors)), 1e-8, label = case$name)
    expect_identical(as.integer(inputs$modal), as.integer(case$modal), label = case$name)
    # The reference assigns by the largest posterior; so does latents.
    expect_identical(as.integer(inputs$modal), ref_modal(inputs$posteriors), label = case$name)
  }))
  expect_identical(length(ml3_kit_cases), 73L)
  expect_identical(length(ml3_sweep_cases), 120L)
})

test_that("step three matches Latent GOLD on latents' posteriors", {
  invisible(lapply(ml3_lg_cases, \(case) {
    fit <- ml3_fits[[case$name]]
    k <- case$spec$k
    lapply(ml3_variances, \(variance) {
      v <- substr(variance, 1L, 1L)
      inner <- ml3_inner(fit, case, variance)
      default <- ml3_lg_result(case$name, paste0(v, "r"), k, variance)
      random <- ml3_lg_result(paste0("x_", case$name), paste0(v, "r"), k, variance)
      expect_false(is.null(default), label = sprintf("%s %s Latent GOLD output", case$name, v))
      # latents must reach the best maximum Latent GOLD found, default start
      # or 100 random start sets.
      best <- max(default$log_likelihood, random$log_likelihood)
      ml3_check(case$name, variance, "Latent GOLD best", "log likelihood shortfall",
                max(0, best - inner$log_likelihood), ml3_tol$log_likelihood)
      same <- Filter(\(lg) !is.null(lg) && abs(lg$log_likelihood - inner$log_likelihood) <
                       ml3_tol$log_likelihood, list(default = default, random = random))
      if (!length(same)) return(NULL)
      lg <- same[[1L]]
      against <- paste("Latent GOLD", names(same)[1L])
      ours <- ml3_run(fit, case, variance = variance, vcov_type = "independent")$result
      pairs <- ml3_run(fit, case, variance = variance, vcov_type = "independent",
                       contrast = "pairs")$result
      ml3_check(case$name, variance, against, "log likelihood",
                inner$log_likelihood - lg$log_likelihood, ml3_tol$log_likelihood)
      ml3_check(case$name, variance, against, "mean (SE units)",
                (ours$estimate - lg$means) / ours$standard_error, ml3_tol$estimate_se)
      # Latent GOLD's robust variance carries the same N / (N - 1) factor.
      ml3_check(case$name, variance, against, "mean SE (relative)",
                ours$standard_error / lg$mean_se - 1, ml3_tol$relative)
      ml3_check(case$name, variance, against, "variance (relative)",
                ours$variance / rep_len(lg$variances, k) - 1, ml3_tol$relative)
      # Relative to the statistic, or to 1 when it is near zero (e_null_01).
      ml3_check(case$name, variance, against, "pair Wald (relative)",
                (pairs$statistic^2 - lg$pair_wald) / pmax(1, lg$pair_wald), ml3_tol$relative)
      first <- pairs$reference_class == 1L
      ml3_check(case$name, variance, against, "contrast with class 1 (SE units)",
                (pairs$estimate[first] - lg$contrasts_first) / pairs$standard_error[first],
                ml3_tol$estimate_se)
      ml3_check(case$name, variance, against, "contrast SE (relative)",
                pairs$standard_error[first] / lg$contrasts_first_se - 1, ml3_tol$relative)
      NULL
    })
  }))
  ml3_assert()
})

test_that("Latent GOLD's model-based and design-based variances match the reference and latents", {
  invisible(lapply(ml3_lg_cases, \(case) {
    fit <- ml3_fits[[case$name]]
    inputs <- ml3_level_inputs(fit, case)
    k <- case$spec$k
    lapply(ml3_variances, \(variance) {
      v <- substr(variance, 1L, 1L)
      inner <- ml3_inner(fit, case, variance)
      standard <- ml3_lg_result(case$name, paste0(v, "s"), k, variance)
      if (abs(standard$log_likelihood - inner$log_likelihood) > ml3_tol$log_likelihood) return(NULL)
      reference <- ref_step3(inputs$posteriors, inputs$values, variance)
      # The standard (inverse observed information) errors are not latents'
      # output; they check the reference's information matrix.
      ml3_check(case$name, variance, "Latent GOLD standard", "reference model SE (relative)",
                reference$classes$standard_error_model / standard$mean_se - 1, ml3_tol$relative)
      ml3_check(case$name, variance, "Latent GOLD standard", "reference model pair Wald (relative)",
                (reference$pairs$wald_model - standard$pair_wald) / pmax(1, standard$pair_wald),
                ml3_tol$relative)
      robust <- ml3_lg_result(case$name, paste0(v, "r"), k, variance)
      ml3_check(case$name, variance, "Latent GOLD robust", "reference variance SE (relative)",
                reference$variance_se_robust / robust$variance_se - 1, ml3_tol$relative)
      if (!isTRUE(case$psu)) return(NULL)
      psu <- ml3_lg_result(case$name, paste0(v, "p"), k, variance)
      ours <- ml3_run(fit, case, variance = variance, vcov_type = "cluster")$result
      pairs <- ml3_run(fit, case, variance = variance, vcov_type = "cluster",
                       contrast = "pairs")$result
      ml3_check(case$name, variance, "Latent GOLD psuid", "cluster mean SE (relative)",
                ours$standard_error / psu$mean_se - 1, ml3_tol$relative)
      ml3_check(case$name, variance, "Latent GOLD psuid", "cluster pair Wald (relative)",
                (pairs$statistic^2 - psu$pair_wald) / pmax(1, psu$pair_wald), ml3_tol$relative)
      NULL
    })
  }))
  ml3_assert()
})

#' latents against the reference on one case and variance model
ml3_against_reference <- function(case, fit, variance, against = "reference") {
  inputs <- ml3_level_inputs(fit, case)
  inner <- ml3_inner(fit, case, variance)
  reference <- ref_step3(inputs$posteriors, inputs$values, variance, ci_level = 0.9,
                         clusters = inputs$clusters)
  # The reference keeps the best of several starts; latents must reach it.
  ml3_check(case$name, variance, against, "log likelihood shortfall",
            max(0, reference$log_likelihood - inner$log_likelihood), ml3_tol$log_likelihood)
  if (abs(reference$log_likelihood - inner$log_likelihood) > ml3_tol$log_likelihood) {
    return(invisible(NULL))
  }
  run <- ml3_run(fit, case, variance = variance, vcov_type = "cluster", ci_level = 0.9)
  ours <- run$result
  if (is.null(ours)) {
    fail(sprintf("%s %s: three_step() failed: %s", case$name, variance,
                 conditionMessage(run$error)))
    return(invisible(NULL))
  }
  pairs <- ml3_run(fit, case, variance = variance, vcov_type = "cluster", ci_level = 0.9,
                   contrast = "pairs")$result
  # A class carried by one unit (sweep_041 under equal variances, finding 2):
  # the reference's sandwich there is zero to rounding, and latents reports
  # no error at all. That is checked as such; its errors, interval and pairs
  # are left out of the comparisons below, and the mean is measured on the
  # model-based error's scale.
  se <- reference$classes$standard_error
  model_se <- reference$classes$standard_error_model
  single <- is.na(ours$standard_error)
  if (any(single)) {
    ml3_check(case$name, variance, against, "single-unit class: reference sandwich / model SE",
              se[single] / model_se[single], 1e-6)
    ml3_check(case$name, variance, against, "single-unit class: latents interval is NA",
              as.numeric(!is.na(c(ours$conf_low[single], ours$conf_high[single]))), 0)
  }
  se <- ifelse(single, model_se, se)
  kept <- !single
  tol <- ml3_tol$reference_relative
  ml3_check(case$name, variance, against, "log likelihood",
            inner$log_likelihood - reference$log_likelihood, ml3_tol$log_likelihood)
  ml3_check(case$name, variance, against, "mean (SE units)",
            (ours$estimate - reference$classes$estimate) / se, ml3_tol$reference_estimate_se)
  ml3_check(case$name, variance, against, "mean SE (relative)",
            ((ours$standard_error - reference$classes$standard_error) / se)[kept], tol)
  ml3_check(case$name, variance, against, "90% CI (SE units)",
            c((ours$conf_low - reference$classes$conf_low) / se,
              (ours$conf_high - reference$classes$conf_high) / se)[c(kept, kept)],
            ml3_tol$reference_estimate_se)
  ml3_check(case$name, variance, against, "variance (relative)",
            ours$variance / reference$classes$variance - 1, tol)
  ml3_check(case$name, variance, against, "effective n (relative)",
            ours$effective_n / reference$classes$effective_n - 1, tol)
  expect_identical(pairs$class, reference$pairs$class)
  expect_identical(pairs$reference_class, reference$pairs$reference_class)
  # A pair with a single-unit class has no error in latents.
  paired <- !is.na(pairs$standard_error)
  expect_identical(paired, !(pairs$class %in% which(single) |
                               pairs$reference_class %in% which(single)))
  pairs <- pairs[paired, , drop = FALSE]
  reference$pairs <- reference$pairs[paired, , drop = FALSE]
  # latents adjusts over the pairs that have a test (p.adjust() drops NA).
  if (!all(paired)) {
    reference$pairs$p_adjusted <- stats::p.adjust(reference$pairs$p_value, "BH")
  }
  pair_se <- reference$pairs$standard_error
  ml3_check(case$name, variance, against, "pair difference (SE units)",
            (pairs$estimate - reference$pairs$estimate) / pair_se, ml3_tol$reference_estimate_se)
  ml3_check(case$name, variance, against, "pair SE (relative)",
            pairs$standard_error / pair_se - 1, tol)
  z <- reference$pairs$statistic
  ml3_check(case$name, variance, against, "pair statistic (relative)",
            (pairs$statistic - z) / pmax(1, abs(z)), 1e-5)
  # p-values to the precision the statistic carries: dp = 2 phi(z) dz, and
  # a BH-adjusted value moves by at most the number of pairs times that.
  dp <- 2 * stats::dnorm(z) * pmax(1, abs(z)) + 1e-300
  ml3_check(case$name, variance, against, "pair p-value (in units of dz)",
            (pairs$p_value - reference$pairs$p_value) / dp, 1e-5)
  ml3_check(case$name, variance, against, "pair p (BH) (in units of dz)",
            (pairs$p_adjusted - reference$pairs$p_adjusted) / (length(z) * max(dp)), 1e-5)
  invisible(NULL)
}

test_that("the independent reference agrees on every kit case", {
  ok <- Filter(\(case) identical(case$expect, "ok"), ml3_kit_cases)
  invisible(lapply(ok, \(case) lapply(ml3_variances, \(variance)
    ml3_against_reference(case, ml3_fits[[case$name]], variance))))
  ml3_assert()
})

test_that("the independent reference agrees on the sweep of 120 datasets", {
  invisible(lapply(ml3_sweep_cases, \(case) lapply(ml3_variances, \(variance) {
    fit <- ml3_fits[[case$name]]
    run <- ml3_run(fit, case, variance = variance, vcov_type = "independent")
    if (!is.null(run$error)) {
      # A class whose outcome variance collapses (one unit, varying
      # variances) has no information: refused, after the boundary warning.
      expect_s3_class(run$error, "latents_singular_information")
      expect_true("latents_boundary" %in% run$warnings, label = case$name)
      return(NULL)
    }
    ml3_against_reference(case, fit, variance, "reference (sweep)")
  })))
  ml3_assert()
})

test_that("on the sweep's multimodal cases latents reaches the best known maximum", {
  # The sweep cases where latents' former single start stopped below the
  # reference's best (finding 1), run through Latent GOLD from its default
  # non-random start and from 100 random start sets (make-random-starts.R).
  local <- sprintf("sweep_%03d", c(18L, 28L, 30L, 83L, 105L, 108L))
  invisible(lapply(ml3_sweep_cases[local], \(case) lapply(ml3_variances, \(variance) {
    v <- paste0(substr(variance, 1L, 1L), "r")
    default <- ml3_lg_result(case$name, v, case$spec$k, variance)
    random <- ml3_lg_result(paste0("x_", case$name), v, case$spec$k, variance)
    inner <- ml3_inner(ml3_fits[[case$name]], case, variance)
    reference <- ref_step3(case$posteriors, case$data$y, variance)
    lg_best <- max(random$log_likelihood, default$log_likelihood)
    # The reference (25 starts) reaches at least Latent GOLD's best maximum.
    ml3_check(case$name, variance, "Latent GOLD random starts", "reference shortfall",
              max(0, lg_best - reference$log_likelihood), ml3_tol$log_likelihood)
    ml3_check(case$name, variance, "best known maximum", "latents log likelihood shortfall",
              max(0, max(lg_best, reference$log_likelihood) - inner$log_likelihood),
              ml3_tol$log_likelihood)
    NULL
  })))
  ml3_assert()
})

test_that("a class mean carried by a single unit does not report certainty", {
  # sweep_041: four profiles on 50 units, one unit alone in profile 3 with
  # posterior 1. Under equal variances that unit's score for the class mean
  # is zero at the estimate, so the sandwich gives a standard error of 0 and a
  # zero-width interval (the model-based error is about 66). latents now
  # reports NA with a latents_small_classes warning: finding 2,
  # repro-zero-se.R.
  case <- ml3_sweep_cases$sweep_041
  run <- ml3_run(ml3_fits$sweep_041, case, variance = "equal", vcov_type = "independent")
  expect_true(length(run$warnings) > 0L || all(run$result$standard_error > 0),
              info = paste(format(run$result$standard_error), collapse = " "))
})

test_that("Latent GOLD end to end: its own step one, then step three", {
  pipeline <- Filter(\(case) isTRUE(case$pipeline), ml3_kit_cases)
  expect_length(pipeline, 6L)
  invisible(lapply(pipeline, \(case) {
    k <- case$spec$k
    step1 <- sprintf("pl_%s_s1", case$name)
    step3 <- sprintf("pl_%s_s3", case$name)
    fit <- ml3_fits[[case$name]]
    lg_step1 <- ml3_lg_step1_loglik(file.path(ml3_returned(), paste0(step1, ".lst")))
    saved <- utils::read.delim(file.path(ml3_returned(), paste0(step1, "_post.txt")),
                               check.names = FALSE)
    expect_identical(as.integer(saved$id), as.integer(case$data$id))
    theirs <- as.matrix(saved[sprintf("Cluster#%d", seq_len(k))])
    # Latent GOLD rescales the posteriors it reads to sum to one.
    theirs <- theirs / rowSums(theirs)
    # Latent GOLD's step one has no variance floor (bayes variances = 0): on
    # lpa_k3_lo_01 its best start is a spike, one class on a single unit
    # (size 1.0000), 5.09 above latents' proper maximum. Such a step one is
    # degenerate and step three is not compared on it.
    degenerate <- min(colSums(theirs)) < 2
    expect_identical(degenerate, identical(case$name, "lpa_k3_lo_01"), label = case$name)
    if (degenerate) {
      expect_gt(lg_step1, fit$log_likelihood)
      return(NULL)
    }
    # Listing to four decimals.
    ml3_check(case$name, "-", "Latent GOLD step one", "step-one log likelihood",
              lg_step1 - fit$log_likelihood, 1e-4)
    matched <- ml3_relabel(unname(fit$subject_posteriors), unname(theirs))
    expect_lt(matched$gap, 1e-4, label = case$name)
    lapply(ml3_variances, \(variance) {
      lg <- ml3_lg_result(step3, paste0(substr(variance, 1L, 1L), "r"), k, variance)
      # On Latent GOLD's own (six-digit) posteriors the reference reproduces
      # its step three as tightly as on latents' posteriors.
      reference <- ref_step3(unname(theirs), saved$y, variance)
      ml3_check(case$name, variance, "Latent GOLD pipeline", "reference log likelihood",
                reference$log_likelihood - lg$log_likelihood, ml3_tol$log_likelihood)
      ml3_check(case$name, variance, "Latent GOLD pipeline", "reference mean (SE units)",
                (reference$classes$estimate - lg$means) / lg$mean_se, ml3_tol$estimate_se)
      ml3_check(case$name, variance, "Latent GOLD pipeline", "reference SE (relative)",
                reference$classes$standard_error / lg$mean_se - 1, ml3_tol$relative)
      # latents end to end: its step one differs from Latent GOLD's by the
      # optimizers' tolerances and the six printed digits (gap < 1e-4), which
      # moves step three by far less than 1e-3 standard errors.
      ours <- ml3_run(fit, case, variance = variance, vcov_type = "independent")$result
      ours <- ours[matched$order, ]
      ml3_check(case$name, variance, "Latent GOLD pipeline", "end-to-end mean (SE units)",
                (ours$estimate - lg$means) / lg$mean_se, 1e-3)
      ml3_check(case$name, variance, "Latent GOLD pipeline", "end-to-end SE (relative)",
                ours$standard_error / lg$mean_se - 1, 1e-3)
      ml3_check(case$name, variance, "Latent GOLD pipeline", "end-to-end variance (relative)",
                ours$variance / rep_len(lg$variances, k) - 1, 1e-3)
      NULL
    })
  }))
  ml3_assert()
})

test_that("Mplus DU3STEP and DE3STEP follow the fixed-size convention without CR1", {
  mplus <- Filter(\(case) isTRUE(case$mplus), ml3_kit_cases)
  expect_length(mplus, 10L)
  invisible(lapply(mplus, \(case) {
    k <- case$spec$k
    m <- ml3_mplus_result(case$name, k)
    expect_false(is.null(m), label = case$name)
    expect_identical(as.integer(m$ids), as.integer(case$data$id))
    fit <- ml3_fits[[case$name]]
    # Step one: Mplus prints its likelihood to three decimals.
    ml3_check(case$name, "-", "Mplus", "step-one log likelihood",
              fit$log_likelihood - m$log_likelihood, 1e-3)
    matched <- ml3_relabel(unname(fit$subject_posteriors), m$posteriors)
    expect_lt(matched$gap, 1e-4)
    n <- nrow(m$posteriors)
    lapply(ml3_variances, \(variance) {
      table <- if (identical(variance, "varying")) m$du3step else m$de3step
      reference <- ref_step3(m$posteriors, case$data$y, variance)
      # Mplus prints three decimals and stops step three at its own (looser)
      # criterion: agreement to 2e-2 standard errors and 5e-4 relative. The
      # free-size model (class sizes estimated) misses by 0.1 to 10 SE.
      ml3_check(case$name, variance, "Mplus", "reference mean (SE units)",
                (reference$classes$estimate - table$means) / table$mean_se, 2e-2)
      # Mplus's MLR sandwich has no N / (N - 1) factor.
      without_cr1 <- reference$classes$standard_error * sqrt((n - 1) / n)
      ml3_check(case$name, variance, "Mplus", "reference SE without CR1 (relative)",
                without_cr1 / table$mean_se - 1, 5e-4)
      # Mplus's "Last Step" class proportions are the average step-three
      # posteriors, latents' effective_n / n; the sizes themselves are fixed.
      ml3_check(case$name, variance, "Mplus", "last-step proportions",
                reference$classes$effective_n / n - table$last_sizes, 5e-5)
      # Compared as z = sqrt(chi-square), on the same 2e-2 scale as the means.
      chi_square <- reference$pairs$statistic^2 * n / (n - 1)
      printed <- if (k > 2L) table$pairs$chi_square else table$overall
      ml3_check(case$name, variance, "Mplus", "pair |z| = sqrt(chi-square)",
                sqrt(chi_square) - sqrt(printed), 2e-2)
      # latents itself, matched to Mplus's class numbering.
      ours <- ml3_run(fit, case, variance = variance, vcov_type = "independent")$result
      ours <- ours[matched$order, ]
      ml3_check(case$name, variance, "Mplus", "latents mean (SE units)",
                (ours$estimate - table$means) / table$mean_se, 2e-2)
      ml3_check(case$name, variance, "Mplus", "latents SE without CR1 (relative)",
                ours$standard_error * sqrt((n - 1) / n) / table$mean_se - 1, 5e-4)
      NULL
    })
  }))
  ml3_assert()
})

test_that("edge cases: transformed outcomes, perfect classification, refusals", {
  fits <- ml3_fits
  base <- ml3_kit_cases$lpa_k2_hi_01
  invisible(lapply(ml3_variances, \(variance) {
    plain <- ml3_run(fits$lpa_k2_hi_01, base, variance = variance)$result
    offset <- ml3_run(fits$e_offset_01, ml3_kit_cases$e_offset_01, variance = variance)$result
    units <- ml3_run(fits$e_units_01, ml3_kit_cases$e_units_01, variance = variance)$result
    ml3_check("e_offset_01", variance, "lpa_k2_hi_01", "mean - 1e6 (SE units)",
              (offset$estimate - 1e6 - plain$estimate) / plain$standard_error, 1e-6)
    ml3_check("e_offset_01", variance, "lpa_k2_hi_01", "SE (relative)",
              offset$standard_error / plain$standard_error - 1, 1e-6)
    ml3_check("e_units_01", variance, "lpa_k2_hi_01", "mean x 1e4 (SE units)",
              (units$estimate * 1e4 - plain$estimate) / plain$standard_error, 1e-6)
    ml3_check("e_units_01", variance, "lpa_k2_hi_01", "SE x 1e4 (relative)",
              units$standard_error * 1e4 / plain$standard_error - 1, 1e-6)
    # With near-perfect classification the ML means are the modal means.
    perfect_case <- ml3_kit_cases$e_perfect_01
    ml <- ml3_run(fits$e_perfect_01, perfect_case, variance = variance)$result
    modal <- as.data.frame(three_step(fits$e_perfect_01, perfect_case$data, "y",
                                      method = "modal"))
    ml3_check("e_perfect_01", variance, "modal method", "ML mean - modal mean (SE units)",
              (ml$estimate - modal$estimate) / ml$standard_error, 1e-6)
    NULL
  }))
  expect_gt(ml3_kit_cases$e_perfect_01$entropy, 0.9999)
  expect_identical(min(ml3_kit_cases$e_empty_01$modal_counts), 0L)
  expect_true(min(ml3_kit_cases$e_few_01$modal_counts) %in% 1:3)
  refusals <- Filter(\(case) !identical(case$expect, "ok"), ml3_kit_cases)
  expect_setequal(names(refusals), c("e_constant_01", "e_na_01", "e_weighted_01"))
  invisible(lapply(ml3_variances, \(variance) {
    run <- function(name) ml3_run(fits[[name]], ml3_kit_cases[[name]], variance = variance)$error
    expect_s3_class(run("e_constant_01"), "latents_bad_outcome")
    expect_s3_class(run("e_weighted_01"), "latents_unsupported_weights")
    # Missing outcomes are a plain contract of three_step(), without a class.
    expect_match(conditionMessage(run("e_na_01")), "must not be missing")
  }))
  # Regular cases produce no warnings.
  regular <- Filter(\(case) !startsWith(case$name, "e_") && identical(case$expect, "ok"),
                    ml3_kit_cases)
  warned <- unlist(lapply(regular, \(case) lapply(ml3_variances, \(variance)
    ml3_run(fits[[case$name]], case, variance = variance)$warnings)))
  expect_length(warned, 0L)
  ml3_assert()
})

test_that("two-level fits: groups as units, and clustered against independent", {
  groups <- Filter(\(case) identical(case$config, "ml_grp"), ml3_kit_cases)
  expect_length(groups, 3L)
  invisible(lapply(groups, \(case) lapply(ml3_variances, \(variance) {
    fit <- ml3_fits[[case$name]]
    clustered <- ml3_run(fit, case, variance = variance)$result
    independent <- ml3_run(fit, case, variance = variance, vcov_type = "independent")$result
    expect_identical(unique(clustered$level), "groups")
    expect_identical(nrow(clustered), case$spec$h)
    # Each group is its own unit at this level.
    expect_identical(clustered$standard_error, independent$standard_error)
    NULL
  })))
  individuals <- Filter(\(case) identical(case$config, "ml_ind"), ml3_kit_cases)
  invisible(lapply(individuals, \(case) {
    fit <- ml3_fits[[case$name]]
    clustered <- ml3_run(fit, case, variance = "varying")$result
    independent <- ml3_run(fit, case, variance = "varying", vcov_type = "independent")$result
    expect_identical(clustered$estimate, independent$estimate)
    # A shared group effect (SD 30) inflates the clustered errors.
    expect_true(all(clustered$standard_error > independent$standard_error), label = case$name)
  }))
})

ml3_summary <- do.call(rbind, ml3_log$rows)
if (!is.null(ml3_summary)) {
  ml3_maxima <- aggregate(difference ~ against + quantity, data = ml3_summary, FUN = max)
  ml3_counts <- aggregate(difference ~ against + quantity, data = ml3_summary, FUN = length)
  ml3_maxima$comparisons <- ml3_counts$difference
  ml3_maxima$failed <- aggregate(I(difference > tolerance) ~ against + quantity, data = ml3_summary,
                                 FUN = sum)[[3L]]
  utils::write.csv(ml3_maxima, file.path(ml3_root(), "maxima.csv"), row.names = FALSE)
  utils::write.csv(ml3_summary, file.path(ml3_root(), "comparisons.csv"), row.names = FALSE)
}
