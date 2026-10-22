# Direct effects of covariates on continuous indicators (`direct_effects`)
# against Latent GOLD 6.1, the Mplus 9 Demo and an independent implementation.
#
# The kit (features/direct-effects/, built by make-kit.R from simulab data and
# run by features/run.sh) holds 67 cases; each is fitted here from three seeds
# and compared with both programs' retained outputs: log likelihood, parameter
# count, BIC, every mean, variance, effect and class proportion with its
# standard error, row-by-row posteriors and predict() on the training data.
# Classes are aligned once from the posteriors. Tolerances are documented in
# de_compare() (features/direct-effects/helpers.R). The independent
# implementation (reference.R) is run on 84 further simulab datasets.

de_root <- file.path("features", "direct-effects")
sys.source(file.path(de_root, "helpers.R"), envir = environment())
sys.source(file.path(de_root, "reference.R"), envir = environment())
de_state <- new.env()
de_results <- function() {
  if (is.null(de_state$external)) {
    de_state$external <- withr::with_seed(20261020, de_compare_all(de_root))
  }
  de_state$external
}

# One expectation per compared row, labelled so a failure names its case.
de_expect_rows <- function(rows) {
  invisible(lapply(seq_len(nrow(rows)), \(i) {
    row <- rows[i, ]
    expect_lte(row$difference, row$tolerance,
               label = sprintf("%s %s %s %s: |%.10g - %.10g|", row$case, row$program,
                               row$quantity, row$item, row$latents, row$reference))
  }))
}

test_that("every kit case has both programs' outputs, and both terminated", {
  cases <- de_cases(de_root)
  expect_length(cases, 67L)
  returned <- file.path(de_root, "returned")
  invisible(lapply(names(cases), \(name) {
    expect_true(all(file.exists(file.path(returned, c(
      paste0(name, c(".lst", "_write.txt", "_posteriors.txt", ".out", "-results.dat",
                     "-posteriors.dat")))))), label = name)
  }))
  expect_length(de_results()$failed, 0L)
})

test_that("log likelihood, parameter count and BIC agree with both programs", {
  rows <- subset(de_results()$comparison,
                 quantity %in% c("log_likelihood", "n_parameters", "bic"))
  expect_setequal(unique(rows$program), c("latentgold", "mplus"))
  de_expect_rows(rows)
})

test_that("posteriors and predict() on the training data agree row by row", {
  de_expect_rows(subset(de_results()$comparison, quantity %in% c("posterior", "predict")))
})

test_that("means, variances, effects and class proportions agree", {
  rows <- subset(de_results()$comparison, endsWith(quantity, "_estimate"))
  expect_true(all(c("mean_estimate", "variance_estimate", "effect_estimate",
                    "proportion_estimate") %in% rows$quantity))
  de_expect_rows(rows)
})

test_that("their standard errors agree", {
  rows <- subset(de_results()$comparison, endsWith(quantity, "_se"))
  # Latent GOLD reports no standard error for the one class of k = 1 (-1).
  rows <- subset(rows, reference >= 0)
  expect_true(all(c("mean_se", "variance_se", "effect_se") %in% rows$quantity))
  de_expect_rows(rows)
})

test_that("latents reaches the best maximum either program found", {
  maxima <- de_results()$maxima
  mine <- subset(maxima, program == "latents")
  # The same tolerance as the log-likelihood comparison: Mplus's own value
  # drifts from the others by up to 3.5e-8 relative at agreeing estimates.
  invisible(lapply(seq_len(nrow(mine)), \(i) {
    expect_gte(mine$gap_to_best[[i]], -1e-6 * abs(mine$best_reference[[i]]),
               label = sprintf("%s latents %.6f vs best reference %.6f", mine$case[[i]],
                               mine$log_likelihood[[i]], mine$best_reference[[i]]))
  }))
  # Every seed of latents reaches the same maximum.
  expect_true(all(mine$seed_spread < 1e-6))
  # Every case is compared quantity by quantity with at least one program.
  compared <- aggregate(compared ~ case, data = subset(maxima, program != "latents"), FUN = any)
  expect_true(all(compared$compared), info = paste(subset(compared, !compared)$case,
                                                   collapse = ", "))
})

test_that("one profile is a linear regression of each indicator", {
  cases <- Filter(\(case) case$spec$k == 1L, de_cases(de_root))
  expect_length(cases, 2L)
  invisible(lapply(cases, \(case) {
    fit <- de_results()$fits[[case$name]]
    models <- lapply(case$spec$vars, \(indicator) {
      covariates <- case$spec$effects$covariate[case$spec$effects$indicator == indicator]
      stats::lm(stats::reformulate(if (length(covariates)) covariates else "1", indicator),
                data = case$data)
    })
    expect_equal(unname(fit$means[1, ]),
                 vapply(models, \(m) stats::coef(m)[["(Intercept)"]], numeric(1)),
                 tolerance = 1e-9)
    expected <- unlist(lapply(seq_len(nrow(case$spec$effects)), \(e) {
      stats::coef(models[[match(case$spec$effects$indicator[[e]], case$spec$vars)]])[[
        case$spec$effects$covariate[[e]]]]
    }))
    expect_equal(unname(fit$direct_effects), expected, tolerance = 1e-9)
    expect_equal(unname(fit$variances[1, ]),
                 vapply(models, \(m) mean(stats::residuals(m)^2), numeric(1)),
                 tolerance = 1e-9)
    expect_equal(fit$log_likelihood,
                 sum(vapply(models, \(m) as.numeric(stats::logLik(m)), numeric(1))),
                 tolerance = 1e-12)
  }))
})

test_that("the estimated effects recover the generating effects", {
  # The collinear cases draw z from the class (z = 2 [class 2] + noise), so
  # P(class | z) depends on z while the model holds it at pi: misspecified by
  # construction, and all three programs agree on b = 1.36 and 1.48 there
  # (truth 0.8). They are compared with the programs above, not with truth.
  cases <- Filter(\(case) !identical(case$truth$covariates$kind, "collinear"),
                  de_cases(de_root))
  rows <- do.call(rbind, lapply(cases, \(case) {
    table <- get_results(de_results()$fits[[case$name]], "direct_effects")
    truth <- merge(table, case$truth$effects, by = c("indicator", "covariate"))
    data.frame(case = case$name, z = (truth$estimate - truth$b) / truth$standard_error)
  }))
  expect_equal(nrow(rows), sum(vapply(cases, \(case) nrow(case$truth$effects), integer(1))))
  # About 80 effects: a |z| above 4.5 among them has probability below 1e-3.
  expect_true(all(abs(rows$z) < 4.5), info = paste(rows$case[abs(rows$z) >= 4.5],
                                                   collapse = ", "))
})

test_that("an independent implementation agrees on 84 more datasets", {
  design <- ref_design()
  expect_gte(nrow(design$grid), 80L)
  battery <- withr::with_seed(20261020, do.call(rbind, lapply(
    seq_len(nrow(design$grid)), \(i) de_quietly(ref_check(design$grid[i, ], design$effect_sets)))))
  invisible(lapply(seq_len(nrow(battery)), \(i) {
    row <- battery[i, ]
    # (a) The same formula at the same point: rounding over at most 2000 rows.
    expect_lte(abs(row$reference_at_latents - row$latents), 1e-12 * abs(row$latents),
               label = sprintf("%s likelihood at latents' estimates", row$name))
    # (b) latents' maximum is not below the independent one beyond 10 times
    # latents' relative stopping tolerance.
    expect_gte(row$latents - row$reference_max, -1e-9 * abs(row$latents),
               label = sprintf("%s latents %.10f vs independent %.10f", row$name,
                               row$latents, row$reference_max))
    if (!row$boundary) {
      # (c) numDeriv's Hessian (Richardson, steps of 0.1 to 0.0125 standard
      # errors) against latents' differences of its analytic score, at points
      # that differ by EM's stopping error.
      expect_lte(row$max_relative_se, 1e-4,
                 label = sprintf("%s standard error of %s", row$name, row$worst_se))
      # (d) The score per standard error at latents' estimate.
      expect_lte(row$max_scaled_score, 1e-3, label = sprintf("%s score", row$name))
    }
  }))
  expect_lte(sum(battery$boundary), 2L)
})
