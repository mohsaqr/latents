# The Vuong-Lo-Mendell-Rubin test (get_results(grid, "lmr")) against Latent
# GOLD 6.1, the Mplus 9 Demo's TECH11 and an independent R reference, on the
# simulab kit of features/vlmr/ (make-kit.R; summarise.R writes the measured
# differences to summary.csv).
#
# - Latent GOLD (`latent Cluster nominal 1:4`) prints each model's log
#   likelihood and, per neighbouring pair, LR, the adjusted LR, their
#   p-values and the reference mean and SD, to 4 decimals. latents must first
#   reach Latent GOLD's maximum (unless Latent GOLD's is degenerate: a share
#   below 1 / n or a zero variance). LR then agrees to the printing plus
#   twice the two log-likelihood gaps. The weights of over-extracted models
#   move by up to 0.02 when estimates move within their printed 4 decimals
#   (the reference evaluated at Latent GOLD's printed estimates shows it), so
#   the moments and p-values are allowed twice that measured movement on top
#   of the printing. Two conventions differ and are tested as such: Latent
#   GOLD's `caseweight` is a frequency weight (B = sum w s s'), latents'
#   weights are sampling weights (B = sum w^2 s s'); and Latent GOLD counts
#   rows with nothing observed in the adjustment's n, latents does not.
# - Mplus TECH11 (estimator MLR) prints the H0 log likelihood, LR, the
#   parameter difference and the adjusted LR to 3 decimals. Its H0 is
#   re-estimated from the K-class model with its first class deleted. Its
#   mean, SD and p-values follow another construction of the weights, not
#   reproduced; they are reported by summarise.R, not tested.
# - The independent reference (features/vlmr/reference.R) recomputes the
#   likelihoods in its own coordinates, the scores and Hessians by numDeriv,
#   Vuong's matrix from its published formula, and the tail probability by
#   CompQuadForm (Imhof, Davies, Farebrother) and by 1e6 Monte Carlo draws.

sys.source(file.path("features", "vlmr", "helpers.R"), envir = environment())
sys.source(file.path("features", "vlmr", "reference.R"), envir = environment())

vlmr_all_cases <- vlmr_cases()
vlmr_external <- names(vlmr_all_cases)[vapply(vlmr_all_cases, `[[`, logical(1), "external")]
vlmr_run <- vlmr_compute_all(names(vlmr_all_cases), cores = 2L)
vlmr_table <- vlmr_run$table

# Printed precision: Latent GOLD 4 decimals, Mplus 3.
lg_print <- 5e-5
mplus_print <- 5e-4
# The reference's numDeriv derivatives (Richardson, in coordinates shifted to
# one at the estimate) agree with latents to 9e-7 in the weights, 1.7e-6 in
# the mean and 1.2e-7 in the p-values over the 45 reference-only datasets;
# repeated with other numDeriv steps they move by ~1e-7. The tolerances are
# ten times the largest measured difference. Stationarity: EM stops on the
# relative likelihood change (1e-12), which leaves summed scores up to ~1e-3
# on slow, flat fits; the weights are compared at the same point, so this
# only screens out fits far from a maximum.
weight_tolerance <- 1e-5
moment_tolerance <- 1e-5
p_tolerance <- 1e-6
score_per_unit <- 1e-5

#' The derivative of P(Q > q) at q, from latents' own weights
vlmr_density <- function(q, lambda, h = 1e-3) {
  abs(.latents_imhof_upper(max(q - h, 0), lambda) - .latents_imhof_upper(q + h, lambda)) /
    (q + h - max(q - h, 0))
}

vlmr_latents_lambda <- function(pair) {
  vlmr_run$lambda[[sprintf("%s:%d", pair$case, pair$null_k)]]$latents
}

test_that("the kit holds 107 simulab datasets, 62 of them run through both programs", {
  expect_length(vlmr_all_cases, 107L)
  expect_length(vlmr_external, 62L)
  returned <- file.path(vlmr_dir, "returned")
  expect_true(all(file.exists(file.path(returned, paste0(vlmr_external, ".lst")))))
  mplus <- unlist(lapply(vlmr_external, function(name) {
    sprintf("%s_k%d.out", name, setdiff(vlmr_all_cases[[name]]$spec$range, 1L))
  }))
  expect_true(all(file.exists(file.path(returned, mplus))))
  expect_length(vlmr_run$errors, 0L)
  # Every listing prints every model of its range, and each Mplus run TECH11.
  printed <- vapply(vlmr_external, function(name) {
    listing <- vlmr_read_lg(name)
    identical(listing$k, vlmr_all_cases[[name]]$spec$range) && !anyNA(listing$ll) &&
      !anyNA(listing$lr[-1L])
  }, logical(1))
  expect_true(all(printed), label = paste(vlmr_external[!printed], collapse = ", "))
  expect_false(anyNA(subset(vlmr_table, case %in% vlmr_external)$mp_lr))
  expect_setequal(unique(vlmr_table$case), names(vlmr_all_cases))
})

test_that("every latents candidate reaches Latent GOLD's non-degenerate maximum", {
  external <- subset(vlmr_table, !is.na(lg_ll1))
  regular <- subset(external, !(lg_degenerate %in% TRUE))
  # Each model appears as a null and as an alternative; check both columns.
  gaps <- c(regular$ll0 - regular$lg_ll0, regular$ll1 - regular$lg_ll1)
  labels <- c(sprintf("%s k=%d", regular$case, regular$null_k),
              sprintf("%s k=%d", regular$case, regular$null_k + 1L))
  invisible(Map(function(gap, label) {
    expect_gt(gap, -1e-4, label = label)
  }, gaps, labels))
  # Latent GOLD's degenerate solutions (a printed share below 1 / n or a zero
  # variance) are VVI single-unit spikes, whose likelihood each program's
  # variance floor sets, so they are not compared. Where latents' maximum is
  # the lower, its floor excluded Latent GOLD's spike; where it is the higher
  # (vvi_small_01-03 at three and four profiles under the relative floor of
  # 2026-10-26, by 0.14 to 3.0), latents' fit is itself bound-active, a spike
  # on its own floor.
  degenerate <- subset(external, lg_degenerate %in% TRUE)
  expect_true(all(degenerate$model == "VVI"))
  above <- c(degenerate$ll0 - degenerate$lg_ll0, degenerate$ll1 - degenerate$lg_ll1) > 1e-4
  on_floor <- c(degenerate$boundary0, degenerate$boundary1)
  expect_true(all(on_floor[above]),
              label = paste(rep(degenerate$case, 2L)[above & !on_floor], collapse = ", "))
})

test_that("LR, adjusted LR, moments and p-values reproduce Latent GOLD", {
  pairs <- subset(vlmr_table, !is.na(lg_lr) & !(lg_degenerate %in% TRUE) &
                    abs(ll0 - lg_ll0) < 1e-4 & abs(ll1 - lg_ll1) < 1e-4)
  expect_gt(nrow(pairs), 130L)
  invisible(lapply(seq_len(nrow(pairs)), function(index) {
    pair <- pairs[index, ]
    label <- sprintf("%s %d vs %d", pair$case, pair$null_k, pair$null_k + 1L)
    # latents refuses a pair whose candidate has bound-active probabilities
    # or a rank-deficient score (documented); Latent GOLD still prints one.
    if (is.na(pair$lr)) {
      expect_match(pair$note, "bound-active|rank deficient", label = paste(label, "refusal"))
      return(invisible(NULL))
    }
    lr_tolerance <- lg_print + 2 * (abs(pair$ll0 - pair$lg_ll0) + abs(pair$ll1 - pair$lg_ll1) +
                                      2 * lg_print)
    expect_lt(abs(pair$lr - pair$lg_lr), lr_tolerance, label = paste(label, "LR"))
    # Latent GOLD's n counts every row read.
    expect_lt(abs(pair$lr / (1 + 1 / (pair$parameters_added * log(pair$n))) -
                    pair$lg_adjusted_lr), lr_tolerance, label = paste(label, "adjusted LR"))
    weighted <- !is.na(pair$ref_freq_mean)
    ours <- if (weighted) {
      c(mean = pair$ref_freq_mean, sd = pair$ref_freq_sd, p = pair$ref_freq_p,
        almr_p = pair$ref_freq_almr_p)
    } else {
      c(mean = pair$reference_mean, sd = pair$reference_sd, p = pair$vlmr_p,
        almr_p = pair$almr_p)
    }
    # Movement of each quantity between latents' estimates and Latent GOLD's
    # printed ones, or under rounding latents' estimates to 4 decimals,
    # whichever is larger (categorical listings are not parsed: rounding only).
    printed <- c(mean = pair$lgpt_mean, sd = pair$lgpt_sd, p = pair$lgpt_p,
                 almr_p = pair$lgpt_almr_p)
    rounding <- c(mean = pair$round_mean, sd = pair$round_sd, p = pair$round_p,
                  almr_p = pair$round_almr_p)
    movement <- pmax(ifelse(is.na(printed), 0, 2 * abs(printed - ours)),
                     ifelse(is.na(rounding), 0, rounding))
    lambda <- vlmr_latents_lambda(pair)
    expect_lt(abs(ours[["mean"]] - pair$lg_mean), lg_print + movement[["mean"]],
              label = paste(label, "mean"))
    expect_lt(abs(ours[["sd"]] - pair$lg_sd), lg_print + movement[["sd"]],
              label = paste(label, "sd"))
    expect_lt(abs(ours[["p"]] - pair$lg_vlmr_p),
              lg_print + vlmr_density(pair$lr, lambda) * lr_tolerance + movement[["p"]],
              label = paste(label, "VLMR p"))
    if (pair$n_units == pair$n) {
      expect_lt(abs(ours[["almr_p"]] - pair$lg_almr_p),
                lg_print + vlmr_density(pair$adjusted_lr, lambda) * lr_tolerance +
                  movement[["almr_p"]], label = paste(label, "aLMR p"))
    }
  }))
})

test_that("Latent GOLD's caseweight is a frequency weight, latents' a sampling weight", {
  weighted <- subset(vlmr_table, config == "weights" & !is.na(ref_freq_p) & !is.na(lg_vlmr_p))
  expect_identical(nrow(weighted), 9L)
  # latents agrees with the sampling-weight reference (tested below) and
  # differs from Latent GOLD wherever the weights matter.
  expect_gt(max(abs(weighted$reference_mean - weighted$lg_mean)), 0.5)
  expect_lt(max(abs(weighted$ref_freq_p - weighted$lg_vlmr_p)), 1e-3)
})

test_that("TECH11's log likelihoods, LR, parameter difference and adjusted LR match", {
  # Pairs at the maxima latents and Latent GOLD share; there Mplus's H0 and
  # H1 must not exceed them, and where Mplus reaches both, LR must agree.
  pairs <- subset(vlmr_table, !is.na(mp_lr) & !is.na(lr) & !(lg_degenerate %in% TRUE) &
                    abs(ll0 - lg_ll0) < 1e-4 & abs(ll1 - lg_ll1) < 1e-4)
  expect_gt(nrow(pairs), 110L)
  slack <- mplus_print + 1e-4
  invisible(lapply(seq_len(nrow(pairs)), function(index) {
    pair <- pairs[index, ]
    label <- sprintf("%s %d vs %d", pair$case, pair$null_k, pair$null_k + 1L)
    expect_true(pair$mp_terminated, label = label)
    expect_identical(as.integer(pair$mp_df), as.integer(pair$parameters_added), label = label)
    expect_lt(pair$mp_h0_ll - pair$ll0, slack, label = paste(label, "H0"))
    expect_lt(pair$mp_h1_ll - pair$ll1, slack, label = paste(label, "H1"))
    if (abs(pair$mp_h0_ll - pair$ll0) < slack && abs(pair$mp_h1_ll - pair$ll1) < slack) {
      lr_tolerance <- mplus_print + 2 * 2 * slack
      expect_lt(abs(pair$lr - pair$mp_lr), lr_tolerance, label = paste(label, "LR"))
      expect_lt(abs(pair$adjusted_lr - pair$mp_adjusted_lr), lr_tolerance,
                label = paste(label, "adjusted LR"))
    }
  }))
  reached <- with(pairs, abs(mp_h0_ll - ll0) < slack & abs(mp_h1_ll - ll1) < slack)
  expect_gte(sum(reached), nrow(pairs) - 2L)
})

test_that("the independent reference reproduces the weights, moments and p-values", {
  pairs <- subset(vlmr_table, !is.na(ref_lr))
  expect_gt(nrow(pairs), 230L)
  expect_gte(length(unique(pairs$case)), 100L)
  invisible(lapply(seq_len(nrow(pairs)), function(index) {
    pair <- pairs[index, ]
    label <- sprintf("%s %d vs %d", pair$case, pair$null_k, pair$null_k + 1L)
    expect_equal(pair$lr, pair$ref_lr, tolerance = 1e-9, label = paste(label, "LR"))
    expect_equal(pair$adjusted_lr, pair$ref_adjusted_lr, tolerance = 1e-9,
                 label = paste(label, "adjusted LR"))
    expect_identical(as.integer(pair$parameters_added), as.integer(pair$ref_parameters_added),
                     label = label)
    expect_identical(as.integer(pair$n_units), as.integer(pair$ref_n_units), label = label)
    expect_identical(pair$n_weights, pair$n_latents_weights, label = label)
    expect_lt(pair$score_sum, score_per_unit * pair$n_units, label = paste(label, "stationary"))
    expect_lt(pair$lambda_gap, weight_tolerance * max(1, pair$lambda_max),
              label = paste(label, "weights"))
    expect_lt(abs(pair$reference_mean - pair$ref_mean), moment_tolerance,
              label = paste(label, "mean"))
    expect_lt(abs(pair$reference_sd - pair$ref_sd), moment_tolerance, label = paste(label, "sd"))
    # The reference p-values carry the reference weights' numerical error.
    expect_lt(abs(pair$vlmr_p - pair$ref_imhof_p), p_tolerance, label = paste(label, "Imhof"))
    expect_lt(abs(pair$vlmr_p - pair$ref_davies_p), p_tolerance, label = paste(label, "Davies"))
    expect_lt(abs(pair$almr_p - pair$ref_davies_almr_p), p_tolerance,
              label = paste(label, "aLMR Davies"))
    if (!is.na(pair$ref_farebrother_p)) {
      expect_lt(abs(pair$vlmr_p - pair$ref_farebrother_p), p_tolerance,
                label = paste(label, "Farebrother"))
    }
    # Monte Carlo: 4.5 binomial standard errors of 1e6 draws.
    margin <- 4.5 * sqrt(max(pair$vlmr_p * (1 - pair$vlmr_p), 1e-6) / pair$ref_draws)
    expect_lt(abs(pair$vlmr_p - pair$ref_mc_p), margin + p_tolerance,
              label = paste(label, "Monte Carlo"))
    margin <- 4.5 * sqrt(max(pair$almr_p * (1 - pair$almr_p), 1e-6) / pair$ref_draws)
    expect_lt(abs(pair$almr_p - pair$ref_mc_almr_p), margin + p_tolerance,
              label = paste(label, "aLMR Monte Carlo"))
  }))
})

test_that("Imhof matches exact probabilities on 460 stressed weight vectors", {
  stress <- vlmr_imhof_stress()
  expect_identical(nrow(stress), 460L)
  exact <- subset(stress, reference_kind == "exact")
  expect_lt(max(exact$abs_error), 1e-10)
  expect_lt(max(subset(stress, reference_kind == "davies")$abs_error), 1e-9)
  expect_lt(max(subset(stress, reference_kind == "unscaled latents")$abs_error), 1e-10)
  expect_true(all(stress$latents >= 0 & stress$latents <= 1))
})

test_that("kit p-values do not move with the indicators' units and origin", {
  columns <- c("lr", "vlmr_p", "adjusted_lr", "almr_p", "reference_mean", "reference_sd")
  invisible(lapply(c("eei_k3poor_01", "vvi_k3_01", "eee_k2_01"), function(name) {
    case <- vlmr_all_cases[[name]]
    fit <- function(data, size = 1) {
      grid <- vlmr_quietly(enumerate_lpa(data, case$spec$vars, n_profiles = 1:3,
                                         model = case$spec$model, n_starts = 20, seed = 1,
                                         tol = 1e-12, max_iter = 5000,
                                         min_variance = 1e-6 * size^2))
      get_results(grid, "lmr")
    }
    base <- fit(case$data)
    vars <- case$spec$vars
    small <- case$data
    small[vars] <- small[vars] * 1e-4
    large <- case$data
    large[vars] <- large[vars] + 1e6
    shifted <- fit(large)
    expect_true(all(is.na(shifted$note)), label = name)
    expect_equal(shifted[columns], base[columns], tolerance = 1e-7, label = name)
    # Rescaled, EM takes another path and stops (relative change 1e-12) with
    # estimates 5e-6 apart (relative); the weights follow them by up to 3e-5,
    # while the independent reference agrees with latents at each endpoint.
    rescaled <- fit(small, 1e-4)
    expect_true(all(is.na(rescaled$note)), label = name)
    expect_equal(rescaled[columns], base[columns], tolerance = 1e-4, label = name)
  }))
})

test_that("weighted and incomplete cases count their informative units", {
  weighted <- subset(vlmr_table, config %in% c("weights", "ref_weights_vvi") & !is.na(lr))
  expect_true(all(weighted$n_units == weighted$n))
  incomplete <- subset(vlmr_table, case == "missing_03" & !is.na(lr))
  expect_true(all(incomplete$n_units == 480L))
  expect_equal(incomplete$adjusted_lr,
               incomplete$lr / (1 + 1 / (incomplete$parameters_added * log(480))),
               tolerance = 1e-14)
})
