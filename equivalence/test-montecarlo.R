# Monte Carlo equivalence: 8 configurations x 15 simulated datasets, each fitted
# by genuine Mplus 9 (Demo) and Latent GOLD 6.1 and by this package. Data,
# inputs and outputs are in montecarlo/ (see its README).
#
# Tolerances follow from what each program prints, plus the difference between
# two optimizers stopping at their own convergence criteria, measured on each
# case as the largest posterior difference `d`:
# - Mplus keeps the log likelihood and BIC to 8 significant digits and the
#   posteriors to 12 decimals; Latent GOLD prints the log likelihood to 4
#   decimals and the posteriors to about 6 significant digits.
# - AvePP is a mean of posteriors, so it can differ by the printed rounding
#   (5e-4) plus `d`; a class count is a sum of n posteriors, so by n * d.
# - AvePP also depends on each row's modal class, which jumps. A row whose two
#   largest posteriors are within 2 d of each other can fall in the other
#   class in the two programs; each such row moves a class's AvePP by at most
#   1 / (its modal size). `modal_slack()` adds that term.

montecarlo_case_names <- names(montecarlo_cases())

#' Largest change in an AvePP that rows near a tie between classes can cause
modal_slack <- function(posteriors, d) {
  sorted <- t(apply(posteriors, 1L, sort, decreasing = TRUE))
  near_tie <- sum(sorted[, 1L] - sorted[, 2L] < 2 * d)
  smallest_class <- min(tabulate(max.col(posteriors, ties.method = "first"),
                                 ncol(posteriors)))
  near_tie / max(smallest_class, 1L)
}

# Fit every case once (seed 1), on several cores: the fits are independent.
montecarlo_fits <- montecarlo_fit_all(montecarlo_cases())

test_that("the kit holds 8 configurations of 15 replicates, each with both outputs", {
  cases <- montecarlo_cases()
  expect_identical(length(cases), 120L)
  expect_identical(as.vector(table(vapply(cases, `[[`, character(1), "config"))),
                   rep(15L, 8L))
  returned <- file.path(montecarlo_dir, "returned")
  expect_true(all(file.exists(file.path(returned, paste0(names(cases), ".out")))))
  expect_true(all(file.exists(file.path(returned, paste0(names(cases), ".lst")))))
  expect_false(any(vapply(montecarlo_fits, inherits, logical(1), "try-error")))
})

test_that("every case reproduces Mplus", {
  cases <- montecarlo_cases()
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- montecarlo_fits[[name]]
    m <- montecarlo_mplus(name, case)
    expect_true(m$terminated && m$replicated, label = name)
    expect_identical(as.integer(fit$n_parameters), m$n_parameters, label = name)
    expect_lt(abs(fit$log_likelihood - m$log_likelihood), 5e-4, label = paste(name, "LL"))
    # Mplus analyses the rows with an observed value, as this package counts.
    observed <- rowSums(!is.na(case$data[case$spec$vars])) > 0L
    expect_identical(fit$n_informative, m$n_observations, label = name)
    expect_false(anyDuplicated(m$id) > 0L)
    expect_setequal(m$id, case$data$id[observed])
    expect_lt(abs(fit$bic_individual - m$bic), 1e-3, label = paste(name, "BIC"))
    rows <- match(m$id, case$data$id)
    matched <- montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE],
                                  m$posteriors$profiles)
    d <- matched$gap
    expect_lt(d, 1e-3, label = paste(name, "posteriors"))
    individuals <- latents_classification(fit, "individuals")
    reference <- m$classification$individuals
    # Single level: Mplus's results file keeps the entropy to 7 decimals.
    if (!is.na(m$entropy_precise)) {
      expect_lt(abs(individuals$entropy - m$entropy_precise), 1e-5 + d,
                label = paste(name, "entropy"))
    } else {
      expect_lt(abs(individuals$entropy - reference$entropy), 5e-4 + d,
                label = paste(name, "entropy"))
    }
    expect_lt(max(abs(individuals$avepp[matched$order] - reference$avepp)),
              5e-4 + d + modal_slack(m$posteriors$profiles, d), label = paste(name, "AvePP"))
    expect_lt(max(abs(fit$effective_profile_counts[matched$order] - m$counts)),
              5e-5 + m$n_observations * d, label = paste(name, "counts"))
    # Profile means: the reference's own posteriors imply them (EM fixed point).
    if (!isTRUE(case$spec$categorical)) {
      values <- as.matrix(case$data[rows, case$spec$vars])
      implied <- montecarlo_implied_means(m$posteriors$profiles, values)
      expect_lt(max(abs(fit$means[matched$order, , drop = FALSE] - implied)),
                5e-3, label = paste(name, "means"))
    }
  }))
})

test_that("two-level cases reproduce Mplus at the group level", {
  cases <- montecarlo_cases()
  two <- names(cases)[vapply(cases, \(case) identical(case$spec$level, "two"), logical(1))]
  expect_identical(length(two), 30L)
  invisible(lapply(two, function(name) {
    case <- cases[[name]]
    fit <- montecarlo_fits[[name]]
    m <- montecarlo_mplus(name, case)
    expect_identical(fit$n_informative_groups, m$n_clusters, label = name)
    # Mplus averages the group level over observed rows, each carrying its
    # group's posterior; applied to this package's group posteriors that gives
    # Mplus's printed figures. The groups are matched by the group posteriors.
    per_row <- fit$group_posteriors[fit$group_index[fit$informative], , drop = FALSE]
    rows <- match(m$id, case$data$id)
    mine <- fit$group_posteriors[fit$group_index[rows], , drop = FALSE]
    matched <- montecarlo_relabel(mine, m$posteriors$groups)
    expect_lt(matched$gap, 1e-3, label = paste(name, "group posteriors"))
    printed <- m$classification$groups
    expect_lt(abs(.multilpa_relative_entropy(per_row) - printed$entropy),
              5e-4 + matched$gap, label = paste(name, "group entropy"))
    expect_lt(max(abs(diag(.multilpa_average_posterior_matrix(per_row))[matched$order] -
                        printed$avepp)),
              5e-4 + matched$gap + modal_slack(m$posteriors$groups, matched$gap),
              label = paste(name, "group AvePP"))
  }))
})

test_that("every case reproduces Latent GOLD, rows with nothing observed included", {
  cases <- montecarlo_cases()
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- montecarlo_fits[[name]]
    lg <- montecarlo_latentgold(name, case)
    expect_identical(as.integer(fit$n_parameters), lg$n_parameters, label = name)
    expect_lt(abs(fit$log_likelihood - lg$log_likelihood), 5e-5 + 2e-4,
              label = paste(name, "LL"))
    expect_false(anyDuplicated(lg$id) > 0L)
    expect_setequal(lg$id, case$data$id)
    rows <- match(lg$id, case$data$id)
    expect_lt(montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE],
                                 lg$profiles)$gap, 1e-3, label = paste(name, "posteriors"))
  }))
})

test_that("no reference program finds a higher maximum than the package", {
  cases <- montecarlo_cases()
  gaps <- vapply(names(cases), function(name) {
    reference <- max(montecarlo_mplus(name, cases[[name]])$log_likelihood,
                     montecarlo_latentgold(name, cases[[name]])$log_likelihood)
    montecarlo_fits[[name]]$log_likelihood - reference
  }, numeric(1))
  expect_true(all(gaps > -5e-4), info = paste(names(gaps)[gaps <= -5e-4], collapse = ", "))
})

test_that("the retained global-maximum audit covers every case and found no lower maximum", {
  # montecarlo/audit-maxima.R: every case from 5 seeds, with 10 and 20 starts.
  audit <- utils::read.csv(file.path(montecarlo_dir, "audit-maxima.csv"))
  expect_setequal(unique(audit$case), names(montecarlo_cases()))
  expect_identical(nrow(audit), 120L * 5L * 2L)
  expect_true(all(audit$reached))
  expect_true(all(audit$gap > -5e-4))
})

test_that("criteria, the AvePP matrix and OCC follow from Mplus's own quantities", {
  # The criteria the package reports are checked against their definitions
  # evaluated at Mplus's log likelihood, parameter count, sample size and
  # posteriors. Entropy is already verified to 1e-5 + d, which bounds the
  # entropy sum EN, and through it ICL, AWE and CLC, by (1e-5 + d) n log K.
  cases <- montecarlo_cases()
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- montecarlo_fits[[name]]
    m <- montecarlo_mplus(name, case)
    rows <- match(m$id, case$data$id)
    matched <- montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE],
                                  m$posteriors$profiles)
    d <- matched$gap
    n <- m$n_observations
    q <- m$n_parameters
    k <- case$spec$k
    p <- m$posteriors$profiles
    entropy_sum <- -sum(p * log(pmax(p, .Machine$double.xmin)))
    deviance <- -2 * m$log_likelihood
    expected <- c(aic = deviance + 2 * q, kic = deviance + 3 * (q + 1),
                  bic_individual = deviance + q * log(n),
                  sabic_individual = deviance + q * log((n + 2) / 24),
                  caic_individual = deviance + q * (log(n) + 1),
                  icl_individual = deviance + q * log(n) + 2 * entropy_sum,
                  awe_individual = deviance + 2 * entropy_sum + 2 * q * (1.5 + log(n)),
                  clc_individual = deviance + 2 * entropy_sum)
    reported <- get_results(fit, "information_criteria", format = "wide")
    tolerance <- 2e-3 + 2 * (1e-5 + d) * n * log(k)
    invisible(lapply(names(expected), function(criterion) {
      expect_lt(abs(reported[[criterion]] - expected[[criterion]]), tolerance,
                label = paste(name, criterion))
    }))
    # The whole AvePP matrix, rows and columns in Mplus's labels.
    path <- file.path(montecarlo_dir, "returned", paste0(name, ".out"))
    lines <- read_mplus_lines(path)
    anchor <- if (identical(case$spec$level, "two"))
      grep("Classification Quality for CW\\s*$", lines)[1L] else
        grep("CLASSIFICATION QUALITY", lines)[1L]
    printed <- mplus_avepp_matrix_after(lines, anchor, k)
    table <- subset(get_results(fit, "average_posteriors"), level == "individuals")
    mine <- matrix(table$average_posterior, k, k, byrow = TRUE)
    expect_lt(max(abs(mine[matched$order, matched$order] - printed)),
              5e-4 + d + modal_slack(p, d), label = paste(name, "AvePP matrix"))
    # OCC from Mplus's printed AvePP and its posterior-based class proportions.
    shares <- m$counts / sum(m$counts)
    occ <- (diag(printed) / (1 - diag(printed))) / (shares / (1 - shares))
    classification <- subset(get_results(fit, "classification"), level == "individuals")
    relative <- abs(classification$odds_correct_classification[matched$order] - occ) / occ
    # A relative error: OCC magnifies the AvePP rounding near 1.
    expect_lt(max(relative), 0.05, label = paste(name, "OCC"))
  }))
})

test_that("predict() on the data reproduces Mplus's posteriors", {
  cases <- montecarlo_cases()
  single <- names(cases)[vapply(cases, \(case) identical(case$spec$level, "single"),
                                logical(1))]
  invisible(lapply(single, function(name) {
    case <- cases[[name]]
    m <- montecarlo_mplus(name, case)
    rows <- match(m$id, case$data$id)
    # Prediction for membership-covariate models is not implemented, and says so.
    if (length(case$spec$covariates)) {
      expect_error(predict(montecarlo_fits[[name]], newdata = case$data),
                   class = "latents_unsupported_prediction")
      return(invisible(NULL))
    }
    predicted <- predict(montecarlo_fits[[name]], newdata = case$data[rows, , drop = FALSE],
                         type = "posterior")
    probabilities <- as.matrix(predicted[grep("^probability_", names(predicted))])
    expect_lt(montecarlo_relabel(probabilities, m$posteriors$profiles)$gap, 1e-3,
              label = paste(name, "predict"))
  }))
})

test_that("drop_empty = FALSE counts every row, as Latent GOLD does", {
  cases <- montecarlo_cases()
  missing_cases <- grep("^lpa_missing_", names(cases), value = TRUE)
  invisible(lapply(missing_cases, function(name) {
    case <- cases[[name]]
    fit <- suppressWarnings(lpa(case$data, case$spec$vars, n_profiles = case$spec$k,
                                model = "VVI", missing = "fiml", n_starts = 20,
                                seed = 1, tol = 1e-10, max_iter = 5000,
                                drop_empty = FALSE))
    lg <- montecarlo_latentgold(name, case)
    rows <- match(lg$id, case$data$id)
    matched <- montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE], lg$profiles)
    expect_equal(sum(fit$effective_profile_counts), nrow(case$data), tolerance = 1e-8)
    expect_lt(max(abs(fit$effective_profile_counts[matched$order] - colSums(lg$profiles))),
              1e-4 + nrow(case$data) * matched$gap, label = paste(name, "counts"))
  }))
})
