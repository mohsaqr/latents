# Local dependence between categorical items (`local_dependence = "a:b"`):
# the joint-items kit (features/joint-items/) against Latent GOLD 6.1
# (`a <-> b | Cluster`, nominal items), the Mplus 9 Demo (the joined item as a
# NOMINAL variable), poLCA (the joined item recoded) and an independent
# reference written from the formula (features/joint-items/reference.R).
# Design, conventions, tolerances and the measured differences:
# features/joint-items/make-kit.R (design), compare.R (tolerances) and
# ji_table (every compared quantity); the findings: findings/*.R there.
# Every seed first uses the original recipe, then resumes its own selected
# parameters at tol <= 1e-13 for verification. Original fits and their gaps are
# retained separately; no production convergence guarantee is strengthened.

ji_env <- environment()
invisible(lapply(c("helpers.R", "reference.R", "compare.R", "generators.R"), \(file)
  sys.source(file.path("features", "joint-items", file), envir = ji_env)))

ji_all_cases <- ji_cases()
ji_all_fits <- ji_fit_all(ji_all_cases)
ji_compared <- Map(ji_compare_case, ji_all_cases, ji_all_fits)
ji_table <- do.call(rbind, lapply(ji_compared, `[[`, "rows"))
ji_notes <- unlist(lapply(ji_compared, `[[`, "notes"))
ji_polished <- ji_parallel(Map(list, ji_all_cases, ji_all_fits), \(x) ji_polish(x[[1]], x[[2]]))
ji_gaps <- do.call(rbind, lapply(ji_polished, `[[`, "row"))
ji_original_diagnostics <- do.call(rbind, Map(ji_verification_diagnostics,
                                            ji_all_cases, ji_all_fits, ji_polished))
equivalence_write_tables("joint-items", list(comparison = ji_table,
  "strict-polish-gaps" = ji_gaps, "original-seed-diagnostics" = ji_original_diagnostics))
# Keep every compared quantity in ji_table without rewriting stored reports.

expect_within <- function(rows) {
  expect_gt(nrow(rows), 0L)
  failing <- rows[!(rows$difference <= rows$tolerance), , drop = FALSE]
  expect_identical(nrow(failing), 0L,
                   info = paste(utils::capture.output(print(failing)), collapse = "\n"))
}

test_that("every case of the kit has its Latent GOLD and Mplus output", {
  expect_length(ji_all_cases, 54L)
  runs <- unlist(lapply(ji_all_cases, \(case) case$spec$latentgold))
  expect_length(runs, 62L)
  expect_true(all(file.exists(file.path(ji_root, "returned", paste0(runs, ".lst")))))
  expect_true(all(file.exists(file.path(ji_root, "returned", paste0(runs, "_posteriors.txt")))))
  mplus <- names(Filter(\(case) case$spec$mplus, ji_all_cases))
  expect_length(mplus, 51L)
  expect_true(all(vapply(mplus, \(name) ji_mplus(name, ji_all_cases[[name]]$spec$k)$terminated,
                         logical(1))))
})

test_that("two seeds of latents reach the same maximum, and no reference a higher one", {
  spread <- vapply(ji_all_fits, \(fit) diff(range(fit$ji_seed_maxima)) /
                     ji_latents_loglik(fit$log_likelihood), numeric(1))
  expect_true(all(spread <= 2), info = paste(names(spread)[spread > 2], collapse = ", "))
  expect_true(all(vapply(ji_all_fits, \(fit) isTRUE(fit$converged), logical(1))))
  expect_within(subset(ji_table, quantity == "log_likelihood"))
  # A reference that stopped lower than latents is reported, not compared.
  expect_identical(ji_notes, character())
})

test_that("the reference likelihood at latents' estimates is latents' likelihood", {
  expect_identical(nrow(ji_gaps), 52L)
  expect_true(all(ji_gaps$relative_at_estimates < 1e-12),
              info = format(max(ji_gaps$relative_at_estimates)))
  expect_true(all(ji_gaps$posteriors_at_estimates < 1e-12),
              info = format(max(ji_gaps$posteriors_at_estimates)))
  expect_true(all(ji_gaps$polished_converged))
})

test_that("the strict own-parameter verification recipe reaches the independent maximum", {
  # The numerical oracle is unchanged. Per-step convergence alone does not
  # bound this remaining gap, so the verification recipe uses tol <= 1e-13.
  short <- ji_gaps[ji_gaps$relative_loglik_gap > 1e-9, c("case", "loglik_gap",
                                                       "relative_loglik_gap",
                                                       "posterior_gap")]
  expect_identical(nrow(short), 0L,
                   info = paste(utils::capture.output(print(short)), collapse = "\n"))
  expect_true(all(ji_gaps$loglik_gap > -1e-9))
})

test_that("every original joint-item seed remains visible beside its verification", {
  expect_equal(nrow(ji_original_diagnostics), 2L * length(ji_all_cases))
  expect_true(all(ji_original_diagnostics$original_n_starts == 15L))
  expect_true(all(ji_original_diagnostics$original_tol == 1e-10))
  expect_true(all(ji_original_diagnostics$verification_n_starts == 1L))
  expect_true(all(ji_original_diagnostics$verification_tol <= 1e-13))
  known <- subset(ji_original_diagnostics, case == "twopair_mod_05")
  expect_true(all(known$original_gap_to_reference > 1e-4))
  expect_true(all(known$original_posterior_gap > .001))
})

test_that("the parameter counts agree under Latent GOLD's full-grid convention", {
  # Latent GOLD's `a <-> b | Cluster` keeps C1 C2 - 1 probabilities per class
  # even for combinations never observed (they go to zero), as does latents.
  # The Mplus kit and poLCA receive only observed joined categories. Their
  # reported counts are increased by k x (unobserved combinations), and their
  # BIC by that count difference x log(N), before comparing full-grid models.
  expect_within(subset(ji_table, quantity == "n_parameters"))
  unobserved <- vapply(ji_all_cases, ji_unobserved_cells, numeric(1))
  expect_true(all(unobserved[grepl("^(sparse|perfect)_", names(unobserved))] > 0))
  expect_identical(unname(unobserved[["perfect_01"]]), 6)
})

test_that("BIC agrees with Latent GOLD and Mplus", {
  expect_within(subset(ji_table, quantity %in% c("bic", "bic_individual")))
})

test_that("class sizes, item and joint probabilities agree with Latent GOLD", {
  expect_within(subset(ji_table, reference == "latentgold" &
                         quantity %in% c("class_sizes", "item_probabilities",
                                         "joint_probabilities")))
  # Unobserved cells printed by Latent GOLD are zero. Whether or not a
  # reference prints them, latents must retain every one at its exact floor.
  expect_within(subset(ji_table, quantity %in%
                         c("joint_unobserved_zero", "joint_unobserved_floor")))
  floor_rows <- subset(ji_table, quantity == "joint_unobserved_floor")
  expected_cells <- sum(vapply(ji_all_cases,
    \(case) case$spec$k * ji_unobserved_cells(case), numeric(1)))
  expect_equal(sum(floor_rows$n), expected_cells)
})

test_that("posteriors agree row by row with Latent GOLD and Mplus", {
  expect_within(subset(ji_table, quantity == "posteriors"))
})

test_that("the probabilities implied by Mplus's posteriors agree", {
  expect_within(subset(ji_table, reference == "mplus" &
                         quantity %in% c("class_sizes", "joint_probabilities")))
})

test_that("a partially observed pair loses the pair, unlike Latent GOLD's includeall", {
  partial <- ji_all_cases[grepl("^partial_", names(ji_all_cases))]
  expect_length(partial, 3L)
  invisible(lapply(names(partial), function(name) {
    case <- partial[[name]]
    fit <- ji_all_fits[[name]]
    expect_true("latents_partial_pair" %in% fit$ji_warnings)
    # latents was compared above with Latent GOLD on data whose partial pairs
    # were made missing. On the raw data Latent GOLD keeps the observed item's
    # marginal, sum_b P(a, b | class); the reference reproduces that maximum.
    raw <- ji_latentgold(case$spec$latentgold[["main"]], case$spec$k)
    marginal <- ref_marginal_em(case$data, c("p1", "p2"), c("u1", "u2", "u3"), case$spec$k,
                                n_starts = 10L, seed = case$seed)
    expect_lt(abs(marginal$log_likelihood - raw$log_likelihood), ji_tol$lg_printed + 1e-6)
    expect_identical(raw$n_parameters, as.integer(fit$n_parameters))
    # Latent GOLD's raw run also scores the lone items, so its likelihood is lower.
    expect_gt(fit$log_likelihood, raw$log_likelihood + 10)
  }))
})

test_that("a truly independent pair gives Latent GOLD's likelihood ratio", {
  independent <- ji_all_cases[grepl("^indep_", names(ji_all_cases))]
  plain_fits <- lapply(independent, ji_fit, local_dependence = FALSE)
  invisible(lapply(names(independent), function(name) {
    case <- independent[[name]]
    plain <- ji_latentgold(case$spec$latentgold[["plain"]], case$spec$k)
    joint <- ji_latentgold(case$spec$latentgold[["main"]], case$spec$k)
    tolerance <- ji_tol$lg_printed + ji_latents_loglik(plain$log_likelihood)
    expect_lt(abs(plain_fits[[name]]$log_likelihood - plain$log_likelihood), tolerance)
    expect_identical(as.integer(plain_fits[[name]]$n_parameters), plain$n_parameters)
    ratio <- 2 * (ji_all_fits[[name]]$log_likelihood - plain_fits[[name]]$log_likelihood)
    expect_lt(abs(ratio - 2 * (joint$log_likelihood - plain$log_likelihood)), 4 * tolerance)
    expect_gte(ratio, 0)
    # k (C1 - 1)(C2 - 1) = 2 association parameters.
    expect_identical(as.integer(ji_all_fits[[name]]$n_parameters -
                                  plain_fits[[name]]$n_parameters), 2L)
  }))
})

test_that("integer weights are duplicated rows, up to latents' rescaling", {
  weighted <- ji_all_cases[grepl("^weights_", names(ji_all_cases))]
  expect_length(weighted, 2L)
  invisible(lapply(names(weighted), function(name) {
    case <- weighted[[name]]
    duplicated_run <- ji_latentgold(case$spec$latentgold[["duplicated"]], case$spec$k)
    weighted_run <- ji_latentgold(case$spec$latentgold[["main"]], case$spec$k)
    # Latent GOLD: a case weight is a repeated row, N included.
    expect_lt(abs(duplicated_run$log_likelihood - weighted_run$log_likelihood),
              2 * ji_tol$lg_printed)
    expect_identical(duplicated_run$n_cases, weighted_run$n_cases)
    rows <- case$data[rep(seq_len(nrow(case$data)), case$data$w),
                      setdiff(names(case$data), "w")]
    duplicated_case <- case
    duplicated_case$data <- rows
    duplicated_case$spec$weights <- NULL
    fit <- ji_fit(duplicated_case)
    expect_lt(abs(fit$log_likelihood - duplicated_run$log_likelihood),
              ji_tol$lg_printed + ji_latents_loglik(fit$log_likelihood))
    # latents' weights sum to the number of rows: the same maximum, rescaled.
    expect_lt(abs(fit$log_likelihood * nrow(case$data) / nrow(rows) -
                    ji_all_fits[[name]]$log_likelihood),
              2 * ji_latents_loglik(fit$log_likelihood))
  }))
})

test_that("the two-level model with a joined pair is compared in full", {
  multilevel <- subset(ji_table, startsWith(case, "multilevel"))
  expect_setequal(unique(multilevel$quantity),
                  c("log_likelihood", "n_parameters", "bic_individual", "posteriors",
                    "class_sizes", "item_probabilities", "joint_probabilities",
                    "joint_unobserved_zero", "joint_unobserved_floor"))
  expect_within(multilevel)
})

test_that("poLCA on the joined item reaches latents' maximum and estimates", {
  skip_if_not_installed("poLCA")
  eligible <- Filter(\(case) is.null(case$spec$h) && is.null(case$spec$weights) &&
                       !identical(case$config, "huge"), ji_all_cases)
  expect_length(eligible, 49L)
  compared <- ji_parallel(Map(list, eligible, ji_all_fits[names(eligible)]),
                          \(x) ji_polca(x[[1]], x[[2]]))
  rows <- do.call(rbind, lapply(compared, `[[`, "rows"))
  original_rows <- do.call(rbind, lapply(compared, \(x) x$original$rows))
  diagnostics <- do.call(rbind, lapply(compared, `[[`, "diagnostics"))
  attempts <- do.call(rbind, lapply(compared, `[[`, "attempts"))
  equivalence_write_tables("joint-items", list(
    "polca-comparison" = rows, "polca-original-cap-comparison" = original_rows,
    "polca-cap-diagnostics" = diagnostics, "polca-attempts" = attempts))
  expect_true(all(diagnostics$maxiter[diagnostics$stage == "original"] == 20000L))
  expect_true(all(diagnostics$maxiter[diagnostics$stage == "verification"] == 200000L))
  expect_true(all(diagnostics$nrep == 10L & diagnostics$tol == 1e-12))
  expect_within(rows)
})

test_that("the independent reference is independent of latents", {
  reference <- new.env()
  sys.source(file.path(ji_root, "reference.R"), envir = reference)
  used <- unique(unlist(lapply(ls(reference), \(name)
    codetools::findGlobals(get(name, reference), merge = TRUE))))
  package <- setdiff(ls(asNamespace("latents"), all.names = TRUE),
                     c(ls(baseenv(), all.names = TRUE), "%||%"))
  expect_identical(intersect(used, package), character())
})

test_that("the reference study: 140 more datasets against the independent EM", {
  study <- ji_reference_study()
  equivalence_write_tables("joint-items", list("reference-study" = study))
  # The study remains in memory; running tests must not rewrite its report.
  expect_identical(nrow(study), 140L)
  # (a) the reference likelihood at latents' estimates is latents' likelihood
  expect_true(all(study$relative_at_estimates < 1e-9),
              info = format(max(study$relative_at_estimates)))
  # (b) latents reaches the reference EM's maximum
  expect_true(all(study$latents_minus_reference > -ji_latents_loglik(study$log_likelihood)),
              info = format(min(study$latents_minus_reference)))
  # (c) the dependence table recomputed from the estimates
  expect_true(all(study$table_difference < 1e-12), info = format(max(study$table_difference)))
  # (d) the joined model contains the plain one
  expect_true(all(study$joined_minus_plain > -ji_latents_loglik(study$log_likelihood)),
              info = format(min(study$joined_minus_plain)))
  expect_true(all(study$n_parameters_match))
  # A truly independent pair: the likelihood ratio is chi-square with
  # k (C1 - 1)(C2 - 1) = 2 degrees of freedom (40 datasets, n = 400).
  ratio <- 2 * pmax(subset(study, config == "indep")$joined_minus_plain, 0)
  expect_gt(stats::ks.test(ratio, "pchisq", df = 2)$p.value, 0.01)
  expect_lt(abs(mean(ratio) - 2), 3 * 2 / sqrt(length(ratio)))
})
