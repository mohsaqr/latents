# The comparisons of the joint-items kit: latents against Latent GOLD, Mplus,
# poLCA and the independent reference (reference.R). Each function returns
# tidy rows (case, reference, quantity, n, difference, tolerance) that the
# test asserts; equivalence/test-features-joint-items.R sources this file
# after helpers.R and reference.R.
#
# Tolerances = what the reference prints + what the reference's own optimizer
# leaves + what latents' optimizer may leave, each measured over this kit
# against the maximum, i.e. latents' estimates polished by the reference EM
# to an absolute change of 1e-12 (`ji_polish`):
# - Latent GOLD prints 4 decimals (within 5e-5 of its estimate) and its
#   posterior file 6 significant digits (5e-7). Its printed probabilities were
#   at most 5.3e-5 and its posteriors 4.8e-5 from the maximum: 6e-5 covers
#   rounding plus its own optimizer.
# - Mplus's RESULTS file holds 8 significant digits (half a unit in the last
#   one, `ji_mplus_rounding`); beyond that its log likelihood was at most
#   4.3e-7 from the maximum (LOGCRITERION 1e-8): 1e-6. Its posteriors (12
#   decimals) were at most 5.0e-5 away (MCONVERGENCE 1e-6): 6e-5.
# - poLCA (tol 1e-12) was within 1.8e-8 of the maximum in log likelihood,
#   2.2e-5 in its probabilities and 6.7e-5 in its posteriors: 1e-7, 3e-5, 1e-4.
# - latents stops at a relative change of 1e-10 per step. Its log likelihood
#   may be 1e-9 (1 + |LL|) below the maximum, ten times that change (the
#   measured worst, outside the one finding, was 5.3e-10); its probabilities
#   were at most 8.7e-5 and its posteriors 4.1e-4 from the maximum (EM moves
#   slowly along flat directions of the likelihood): 1e-4 and 5e-4. The one
#   case beyond these, twopair_mod_05, is a finding (findings/); its
#   comparisons are expected to fail.

ji_tol <- list(lg_printed = 5e-5, lg_estimate = 6e-5, lg_posterior = 6e-5,
               mplus_loglik = 1e-6, mplus_estimate = 6e-5, mplus_posterior = 6e-5,
               polca_loglik = 1e-7, polca_estimate = 3e-5, polca_posterior = 1e-4,
               latents_estimate = 1e-4, latents_posterior = 5e-4)

#' latents' allowed distance from the maximum in log likelihood
ji_latents_loglik <- function(log_likelihood) 1e-9 * (1 + abs(log_likelihood))

#' Half a unit in the 8th significant digit of a value Mplus saved
ji_mplus_rounding <- function(value) 0.5 * 10^(floor(log10(abs(value))) + 1 - 8)

#' Compare one case with Latent GOLD and Mplus
#'
#' Classes are aligned once per reference, from the posteriors, and the
#' alignment is reused for every class-level quantity of that reference.
#' A reference that stopped at a lower maximum than latents is reported in
#' `notes` and its estimates are not compared.
ji_compare_case <- function(case, fit) {
  spec <- case$spec
  k <- spec$k
  latents_ll <- ji_latents_loglik(fit$log_likelihood)
  rows <- list()
  notes <- character()
  add <- function(...) rows[[length(rows) + 1L]] <<- ji_row(case$name, ...)
  # Latent GOLD's main run: the mask run for partial pairs (latents' convention).
  run <- if ("mask" %in% names(spec$latentgold)) spec$latentgold[["mask"]] else
    spec$latentgold[["main"]]
  lg <- ji_latentgold(run, k)
  at <- match(case$data$id, lg$id)
  stopifnot("Latent GOLD lost rows" = !anyNA(at))
  # Case weights: Latent GOLD sums raw weights (N = their sum); latents
  # rescales them to sum to the number of rows.
  scale <- if (is.null(spec$weights)) 1 else sum(case$data$w) / nrow(case$data)
  # The kit codes Mplus's joined nominal items using observed combinations
  # only. Each omitted cell accounts for one probability per class; latents
  # and Latent GOLD both count those probabilities in their full grids.
  offset <- ji_unobserved_cells(case) * k
  observed_levels <- ref_code(case$data, spec$vars, spec$pairs)$levels
  ll_tol <- ji_tol$lg_printed + scale * latents_ll
  higher <- fit$log_likelihood * scale - lg$log_likelihood > ll_tol
  if (higher) notes <- c(notes, sprintf("%s: Latent GOLD stopped lower (%.4f vs %.4f)",
                                        case$name, lg$log_likelihood,
                                        fit$log_likelihood * scale))
  add("latentgold", "log_likelihood", fit$log_likelihood * scale, lg$log_likelihood,
      if (higher) Inf else ll_tol)
  add("latentgold", "n_parameters", fit$n_parameters, lg$n_parameters, 0)
  # Latent GOLD's BIC counts every case (individuals in a two-level model, the
  # sum of the weights with case weights) and its full-grid parameters.
  if (is.null(spec$h)) {
    expected_bic <- if (scale == 1) lg$bic else
      -2 * lg$log_likelihood / scale + lg$n_parameters * log(nrow(case$data))
    add("latentgold", "bic", fit$bic, expected_bic,
        if (higher) Inf else 2 * (ji_tol$lg_printed + latents_ll) + ji_tol$lg_printed)
  } else {
    add("latentgold", "bic_individual", fit$bic_individual, lg$bic,
        if (higher) Inf else 2 * (ji_tol$lg_printed + latents_ll) + ji_tol$lg_printed)
  }
  if (!higher) {
    relabel <- montecarlo_relabel(fit$subject_posteriors, lg$posteriors[at, , drop = FALSE])
    mine <- ji_latents(fit, relabel$order)
    add("latentgold", "posteriors", mine$posteriors, lg$posteriors[at, ],
        ji_tol$lg_posterior + ji_tol$latents_posterior)
    add("latentgold", "class_sizes", mine$sizes, lg$sizes,
        ji_tol$lg_estimate + ji_tol$latents_estimate)
    others <- setdiff(names(lg$estimated), names(mine$joint))
    add("latentgold", "item_probabilities",
        unlist(lapply(others, \(v) mine$items[[v]][, colnames(lg$estimated[[v]])])),
        unlist(lapply(others, \(v) lg$estimated[[v]])),
        ji_tol$lg_estimate + ji_tol$latents_estimate)
    add("latentgold", "joint_probabilities",
        unlist(lapply(names(mine$joint), \(p)
          mine$joint[[p]][, colnames(lg$estimated[[p]]), drop = FALSE])),
        unlist(lapply(names(mine$joint), \(p) lg$estimated[[p]])),
        ji_tol$lg_estimate + ji_tol$latents_estimate)
    # Preserve the zero check on any unobserved cells Latent GOLD prints.
    # Even when a reference omits a cell, latents must estimate it at the
    # floor. Identify missing cells from the data, not table set differences,
    # which would make the zero check vacuous now that latents has a full grid.
    unseen <- ji_unobserved_probabilities(case, lg$estimated)
    add("latentgold", "joint_unobserved_zero", numeric(length(unseen)), unseen,
        ji_tol$lg_printed)
    unseen <- ji_unobserved_probabilities(case, mine$joint)
    add("latentgold", "joint_unobserved_floor", unseen,
        rep(fit$min_probability, length(unseen)), 0)
  }
  if (isTRUE(spec$mplus)) {
    mp <- ji_mplus(case$name, k)
    at <- match(case$data$id, mp$id)
    stopifnot("Mplus lost rows" = !anyNA(at), "Mplus did not terminate normally" = mp$terminated)
    ll_tol <- ji_mplus_rounding(mp$log_likelihood) + ji_tol$mplus_loglik + latents_ll
    higher <- fit$log_likelihood - mp$log_likelihood > ll_tol
    if (higher) notes <- c(notes, sprintf("%s: Mplus stopped lower (%.4f vs %.4f)", case$name,
                                          mp$log_likelihood, fit$log_likelihood))
    add("mplus", "log_likelihood", fit$log_likelihood, mp$log_likelihood,
        if (higher) Inf else ll_tol)
    # Translate Mplus's observed-grid count/BIC to the full-grid convention;
    # the likelihood and all likelihood/posterior tolerances are unchanged.
    add("mplus", "n_parameters", fit$n_parameters, mp$n_parameters + offset, 0)
    add("mplus", "bic", fit$bic, mp$bic + offset * log(nrow(case$data)),
        if (higher) Inf else ji_mplus_rounding(mp$bic) + 2 * ll_tol)
    if (!higher) {
      posteriors <- mp$posteriors[at, , drop = FALSE]
      relabel <- montecarlo_relabel(fit$subject_posteriors, posteriors)
      mine <- ji_latents(fit, relabel$order)
      add("mplus", "posteriors", mine$posteriors, posteriors,
          ji_tol$mplus_posterior + ji_tol$latents_posterior)
      weights <- if (is.null(spec$weights)) rep(1, nrow(case$data)) else case$data$w
      add("mplus", "class_sizes", mine$sizes,
          colSums(posteriors * weights) / sum(weights),
          ji_tol$mplus_estimate + ji_tol$latents_estimate)
      # The joint probabilities Mplus's posteriors imply (EM's fixed point).
      implied <- unlist(lapply(names(mine$joint), function(p) {
        items <- strsplit(p, ":", fixed = TRUE)[[1L]]
        label <- paste(case$data[[items[1]]], case$data[[items[2]]], sep = ":")
        observed <- !is.na(case$data[[items[1]]]) & !is.na(case$data[[items[2]]])
        vapply(observed_levels[[p]], function(cell) {
          colSums(posteriors[observed, , drop = FALSE] * weights[observed] *
                    (label[observed] == cell)) /
            colSums(posteriors[observed, , drop = FALSE] * weights[observed])
        }, numeric(k))
      }))
      add("mplus", "joint_probabilities",
          unlist(lapply(names(mine$joint), \(p)
            mine$joint[[p]][, observed_levels[[p]], drop = FALSE])), implied,
          ji_tol$mplus_estimate + ji_tol$latents_estimate)
    }
  }
  list(rows = do.call(rbind, rows), notes = notes)
}

#' Polish a latents fit with the reference EM and measure how far it moves
#'
#' The reference EM starts at latents' estimates and runs to an absolute
#' change of 1e-12. The distance between latents' reported quantities and the
#' polished ones is latents' optimizer gap, per quantity; it enters the
#' tolerances of the comparisons with the other programs. Single-level fits
#' only (the reference has no group level).
#' @return A list: `row` (one summary row; NULL for a two-level fit),
#'   `gap` (per quantity) and `polished` (priors, named probabilities,
#'   posteriors, log likelihood), classes in latents' order.
ji_polish <- function(case, fit) {
  if (!isTRUE(fit$single_level)) {
    # No group level in the reference: the gap is measured against latents
    # itself refitted from the same seeds to a change of 1e-14.
    tight <- ji_fit(case, tol = 1e-14, max_iter = 200000)
    order <- montecarlo_relabel(tight$subject_posteriors, fit$subject_posteriors)$order
    probability_gap <- function(set) max(abs(unlist(lapply(set, \(v)
      tight$response_probabilities[[v]][order, , drop = FALSE] - fit$response_probabilities[[v]]))))
    pairs <- grep(":", fit$vars, fixed = TRUE, value = TRUE)
    sizes <- \(f, o = seq_len(f$n_profiles)) c(f$group_probabilities %*% f$profile_probabilities)[o]
    gap <- list(log_likelihood = tight$log_likelihood - fit$log_likelihood,
                sizes = max(abs(sizes(tight, order) - sizes(fit))),
                items = probability_gap(setdiff(fit$vars, pairs)), joint = probability_gap(pairs),
                posteriors = max(abs(tight$subject_posteriors[, order] - fit$subject_posteriors)))
    return(list(row = NULL, gap = gap, polished = NULL))
  }
  spec <- case$spec
  coded <- ref_code(case$data, spec$vars, spec$pairs)
  priors <- c(fit$profile_probabilities[1L, ])
  indicators <- colnames(coded$codes)
  probabilities <- lapply(indicators,
                          \(v) fit$response_probabilities[[v]][, coded$levels[[v]], drop = FALSE])
  weights <- if (is.null(spec$weights)) NULL else
    case$data$w * nrow(case$data) / sum(case$data$w)
  at_fit <- ref_loglik(coded$codes, priors, probabilities, weights)
  joint_at_fit <- ref_log_joint(coded$codes, priors, probabilities)
  posteriors_at_fit <- exp(joint_at_fit - ref_log_sum_exp(joint_at_fit))
  polished <- ref_em(coded$codes, coded$levels, spec$k, tol = 1e-12, max_iter = 200000L,
                     weights = weights, start = list(priors = priors, probabilities = probabilities))
  names(polished$probabilities) <- indicators
  invisible(lapply(indicators, \(v) colnames(polished$probabilities[[v]]) <<- coded$levels[[v]]))
  pairs <- indicators[grepl(":", indicators, fixed = TRUE)]
  difference <- function(set) {
    if (length(set) == 0L) return(0)
    max(abs(unlist(lapply(set, \(v) polished$probabilities[[v]] -
                            fit$response_probabilities[[v]][, coded$levels[[v]], drop = FALSE]))))
  }
  gap <- list(log_likelihood = polished$log_likelihood - fit$log_likelihood,
              sizes = max(abs(polished$priors - priors)),
              items = difference(setdiff(indicators, pairs)),
              joint = difference(pairs),
              posteriors = max(abs(polished$posteriors - fit$subject_posteriors)))
  row <- data.frame(case = case$name, log_likelihood = fit$log_likelihood,
                    relative_at_estimates = abs(at_fit - fit$log_likelihood) /
                      abs(fit$log_likelihood),
                    posteriors_at_estimates = max(abs(posteriors_at_fit - fit$subject_posteriors)),
                    loglik_gap = gap$log_likelihood,
                    relative_loglik_gap = gap$log_likelihood / (1 + abs(fit$log_likelihood)),
                    parameter_gap = max(unlist(gap[c("sizes", "items", "joint")])),
                    posterior_gap = gap$posteriors, polished_converged = polished$converged)
  list(row = row, gap = gap,
       polished = list(log_likelihood = polished$log_likelihood, priors = polished$priors,
                       probabilities = polished$probabilities,
                       posteriors = polished$posteriors))
}

#' poLCA on the data with the pair replaced by the reference's own joined item
ji_polca <- function(case, fit, polish = NULL, nrep = 10L) {
  spec <- case$spec
  coded <- ref_code(case$data, spec$vars, spec$pairs)
  frame <- as.data.frame(coded$codes)
  names(frame) <- paste0("V", seq_len(ncol(frame)))
  formula <- stats::as.formula(sprintf("cbind(%s) ~ 1", paste(names(frame), collapse = ", ")))
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(case$seed)
  captured <- utils::capture.output(polca <- poLCA::poLCA(
    formula, frame, nclass = spec$k, maxiter = 20000, tol = 1e-12, na.rm = FALSE,
    nrep = nrep, verbose = FALSE, calc.se = FALSE))
  relabel <- montecarlo_relabel(fit$subject_posteriors, polca$posterior)
  order <- relabel$order
  mine <- unlist(lapply(colnames(coded$codes), \(v)
    fit$response_probabilities[[v]][order, coded$levels[[v]], drop = FALSE]))
  theirs <- unlist(lapply(polca$probs, identity))
  latents_ll <- ji_latents_loglik(fit$log_likelihood)
  # ref_code supplies contiguous observed-category codes to poLCA, just as
  # the stored Mplus kit does. Translate its count and BIC to the full grid;
  # compare its probabilities only on those observed codes and separately
  # require every omitted latents cell to be exactly at its probability floor.
  offset <- ji_unobserved_cells(case) * spec$k
  unseen <- ji_unobserved_probabilities(case, fit$response_probabilities)
  # poLCA's own distance from the maximum, for the README (not a tolerance).
  own <- if (!is.null(polish$polished)) {
    p <- polish$polished
    polished <- unlist(lapply(colnames(coded$codes), \(v) p$probabilities[[v]][order, ]))
    data.frame(case = case$name, loglik = abs(polca$llik - p$log_likelihood),
               estimate = max(abs(c(polished - theirs, p$priors[order] - polca$P))),
               posterior = max(abs(p$posteriors[, order] - polca$posterior)))
  }
  list(log_likelihood = polca$llik, own = own,
       rows = rbind(
         ji_row(case$name, "polca", "n_parameters", fit$n_parameters,
                polca$npar + offset, 0),
         ji_row(case$name, "polca", "bic", fit$bic,
                polca$bic + offset * log(nrow(case$data)),
                2 * (latents_ll + ji_tol$polca_loglik)),
         ji_row(case$name, "polca", "joint_unobserved_floor", unseen,
                rep(fit$min_probability, length(unseen)), 0),
         ji_row(case$name, "polca", "log_likelihood", fit$log_likelihood, polca$llik,
                latents_ll + ji_tol$polca_loglik),
         ji_row(case$name, "polca", "class_sizes", fit$profile_probabilities[1L, order],
                polca$P, ji_tol$latents_estimate + ji_tol$polca_estimate),
         ji_row(case$name, "polca", "probabilities", mine, theirs,
                ji_tol$latents_estimate + ji_tol$polca_estimate),
         ji_row(case$name, "polca", "posteriors", fit$subject_posteriors[, order],
                polca$posterior, ji_tol$latents_posterior + ji_tol$polca_posterior)))
}

#' The reference study design: 140 cheap simulab datasets
ji_study_design <- function() {
  grid <- rbind(
    data.frame(config = "bin_small", n = 80L, k = 2L, replicate = 1:20),
    data.frame(config = "bin_small", n = 300L, k = 2L, replicate = 1:20),
    data.frame(config = "cat3_mod", n = 500L, k = 3L, replicate = 1:20),
    data.frame(config = "bin4_mod", n = 400L, k = 2L, replicate = 1:20),
    data.frame(config = "twopair_mod", n = 500L, k = 3L, replicate = 1:20),
    data.frame(config = "indep", n = 400L, k = 2L, replicate = 1:40))
  grid$seed <- 50000L + seq_len(nrow(grid))
  grid$pairs <- ifelse(grid$config == "twopair_mod", "p1:p2,q1:q2", "p1:p2")
  grid
}

#' One dataset of the reference study
#'
#' latents (joined and plain, 8 starts each) against the reference EM
#' (8 starts): (a) the reference likelihood at latents' estimates, (b) the
#' two maxima, (c) the dependence table recomputed from latents' estimates,
#' (d) joined against plain.
ji_study_one <- function(row) {
  drawn <- generate[[row$config]](row$seed, n = row$n)
  data <- drawn$items
  vars <- names(data)
  pairs <- strsplit(row$pairs, ",", fixed = TRUE)[[1L]]
  quiet <- function(expr) withCallingHandlers(expr, message = \(m) invokeRestart("muffleMessage"))
  joined <- quiet(lca(data, vars, row$k, local_dependence = pairs, n_starts = 8L, seed = 1L,
                      tol = 1e-10, max_iter = 5000))
  plain <- quiet(lca(data, vars, row$k, n_starts = 8L, seed = 1L, tol = 1e-10, max_iter = 5000))
  coded <- ref_code(data, vars, pairs)
  priors <- c(joined$profile_probabilities[1L, ])
  probabilities <- lapply(colnames(coded$codes),
                          \(v) joined$response_probabilities[[v]][, coded$levels[[v]], drop = FALSE])
  at_estimates <- ref_loglik(coded$codes, priors, probabilities)
  reference <- ref_em(coded$codes, coded$levels, row$k, n_starts = 8L, seed = row$seed,
                      tol = 1e-10)
  table <- get_results(joined, "local_dependence")
  recomputed <- do.call(rbind, lapply(pairs, function(p) {
    out <- ref_dependence_table(joined$response_probabilities[[p]])
    out$pair <- p
    out
  }))
  key <- \(t) paste(t$pair, t$profile, t$category_1, t$category_2)
  matched <- recomputed[match(key(table), key(recomputed)), ]
  stopifnot("tables differ in their rows" = nrow(table) == nrow(recomputed) &&
              !anyNA(matched$probability))
  data.frame(config = row$config, n = row$n, seed = row$seed,
             log_likelihood = joined$log_likelihood,
             relative_at_estimates = abs(at_estimates - joined$log_likelihood) /
               abs(joined$log_likelihood),
             latents_minus_reference = joined$log_likelihood - reference$log_likelihood,
             table_difference = max(abs(c(table$probability - matched$probability,
                                          table$independent - matched$independent,
                                          table$ratio - matched$ratio))),
             joined_minus_plain = joined$log_likelihood - plain$log_likelihood,
             # The independent EM uses observed cells (zero omitted mass).
             # Its model count needs one parameter per omitted cell/class
             # to match the full-grid model; its likelihood is unchanged.
             n_parameters_match = identical(as.integer(joined$n_parameters),
               as.integer(ref_n_parameters(coded$levels, row$k) + row$k *
                 ji_unobserved_cells(list(data = data, spec = list(pairs = pairs))))))
}

ji_reference_study <- function() {
  design <- ji_study_design()
  rows <- split(design, seq_len(nrow(design)))
  do.call(rbind, ji_parallel(rows, function(row) ji_study_one(row)))
}
