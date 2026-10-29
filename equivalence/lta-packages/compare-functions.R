# Functions for compare.R: lta() on LMest's RLMSdat (one five-category item,
# seven occasions, 1,718 people; job satisfaction in the Russia Longitudinal
# Monitoring Survey) against the current LMest and seqHMM, and against the
# published analysis of Bartolucci, Pandolfi and Pennoni (2017, JSS 81(4),
# section 3.1): a basic latent Markov model with time-homogeneous transitions
# and the same response probabilities at every occasion.
#
#
# And LMest's data_market_sim (200 customers, five occasions; brand, four
# categories, and price, five ranges, as two categorical items) against LMest
# and the hidden Markov models of Pennoni, Pandolfi and Bartolucci (2024, The
# R Journal, RJ-2024-036, section 3.3: lmestSearch(), time-homogeneous
# transitions, two to four states; the one-state model is not a transition
# model and lta() refuses it).
#
# Tolerances, fixed before the comparison: software, the log likelihood to
# 1e-6 plus 1e-10 of its size and every probability to 1e-4; published
# values, half a unit of the last printed digit.

rlms_long <- function() {
  wide <- get(utils::data(list = "RLMSdat", package = "LMest", envir = environment()),
              envir = environment())
  data.frame(id = rep(seq_len(nrow(wide)), each = ncol(wide)),
             time = rep(seq_len(ncol(wide)), times = nrow(wide)),
             y = as.vector(t(as.matrix(wide))))
}

#' The latents fit and its probabilities: initial (k), transition (k x k,
#' rows the origin) and response (k x categories)
fit_latents_lta <- function(data, states) {
  fit <- lta(data, "y", "id", states, time = "time", categorical = "y", n_starts = 50L,
             seed = 1L, tol = 1e-12)
  theta <- stats::coef(fit)
  categories <- as.character(sort(unique(data$y)))
  initial <- vapply(seq_len(states), function(j)
    theta[[sprintf("profile.initial_probability.profile_%d.group_class_1", j)]], numeric(1))
  transition <- t(vapply(seq_len(states), function(from) vapply(seq_len(states), function(to)
    theta[[sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from)]],
    numeric(1)), numeric(states)))
  response <- t(vapply(seq_len(states), function(j) vapply(categories, function(c)
    theta[[sprintf("measurement.response.profile_%d.y:%s", j, c)]], numeric(1)),
    numeric(length(categories))))
  list(fit = fit, log_likelihood = as.numeric(stats::logLik(fit)), n_parameters = fit$n_parameters,
       bic = stats::BIC(fit), initial = initial, transition = transition, response = response)
}

fit_lmest <- function(data, states) {
  fits <- lapply(seq_len(20L), function(s) {
    printed <- utils::capture.output(fitted <- LMest::lmest(
      responsesFormula = y ~ NULL, index = c("id", "time"), data = data, k = states,
      modBasic = 1, start = if (s == 1L) 0 else 1, seed = s, tol = 1e-12, maxit = 100000L,
      out_se = FALSE))
    fitted
  })
  best <- fits[[which.max(vapply(fits, function(f) f$lk, numeric(1)))]]
  list(log_likelihood = best$lk, initial = as.numeric(best$piv),
       transition = unname(best$Pi[, , 2L]), response = unname(t(best$Psi[, , 1L])))
}

fit_seqhmm <- function(data, states) {
  wide <- matrix(data$y, ncol = max(data$time), byrow = TRUE)
  categories <- as.character(sort(unique(data$y)))
  sequences <- TraMineR::seqdef(as.data.frame(wide), alphabet = categories, cpal = rep("#999999", length(categories)))
  model <- withr::with_seed(1L, seqHMM::build_hmm(sequences, n_states = states))
  printed <- utils::capture.output(fit <- withr::with_seed(1L, seqHMM::fit_model(
    model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
    control_em = list(restart = list(times = 20L), reltol = 1e-14, maxeval = 100000L))))
  m <- fit$model
  list(log_likelihood = as.numeric(stats::logLik(m)), initial = unname(as.numeric(m$initial_probs)),
       transition = unname(m$transition_probs), response = unname(m$emission_probs))
}

market_long <- function() {
  get(utils::data(list = "data_market_sim", package = "LMest", envir = environment()),
      envir = environment())[c("id", "time", "brand", "price")]
}

#' The latents fit of the market data: initial, transition and the response
#' probabilities of both items side by side (states x 4 + 5)
fit_latents_market <- function(data, states) {
  items <- c("brand", "price")
  fit <- lta(data, items, "id", states, time = "time", categorical = items, n_starts = 50L,
             seed = 1L, tol = 1e-12)
  theta <- stats::coef(fit)
  initial <- vapply(seq_len(states), function(j)
    theta[[sprintf("profile.initial_probability.profile_%d.group_class_1", j)]], numeric(1))
  transition <- t(vapply(seq_len(states), function(from) vapply(seq_len(states), function(to)
    theta[[sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from)]],
    numeric(1)), numeric(states)))
  labels <- unlist(lapply(items, function(item) paste0(item, ":", sort(unique(data[[item]])))))
  response <- t(vapply(seq_len(states), function(j) vapply(labels, function(label)
    theta[[sprintf("measurement.response.profile_%d.%s", j, label)]], numeric(1)),
    numeric(length(labels))))
  list(fit = fit, log_likelihood = as.numeric(stats::logLik(fit)), n_parameters = fit$n_parameters,
       aic = stats::AIC(fit), bic = stats::BIC(fit), initial = initial, transition = transition,
       response = response)
}

#' LMest on the market data, best of 20 starts; Psi is categories x states x
#' items, padded with NA for brand's missing fifth category
fit_lmest_market <- function(data, states) {
  fits <- lapply(seq_len(20L), function(s) {
    printed <- utils::capture.output(fitted <- LMest::lmest(
      responsesFormula = brand + price ~ NULL, index = c("id", "time"), data = data, k = states,
      modBasic = 1, start = if (s == 1L) 0 else 1, seed = s, tol = 1e-12, maxit = 100000L,
      out_se = FALSE))
    fitted
  })
  best <- fits[[which.max(vapply(fits, function(f) f$lk, numeric(1)))]]
  psi <- best$Psi
  response <- cbind(t(matrix(psi[1:4, , 1L], 4L)), t(matrix(psi[1:5, , 2L], 5L)))
  list(log_likelihood = best$lk, initial = as.numeric(best$piv),
       transition = unname(best$Pi[, , 2L]), response = unname(response))
}

compare_published_market <- function(fits) {
  source <- "Pennoni, Pandolfi and Bartolucci (2024), R Journal RJ-2024-036, section 3.3"
  printed <- data.frame(k = 2:4, ll = c("-2520.131", "-2445.538", "-2434.262"), np = c(17, 29, 43),
                        aic = c("5074.263", "4949.075", "4954.524"),
                        bic = c("5130.334", "5044.727", "5096.352"))
  do.call(rbind, lapply(seq_len(nrow(printed)), function(i) {
    fit <- fits[[sprintf("k%d", printed$k[[i]])]]
    case <- sprintf("market_k%d", printed$k[[i]])
    rbind(
      .printed(case, source, "log likelihood", printed$ll[[i]], fit$log_likelihood),
      data.frame(case = case, reference = source, quantity = "number of parameters",
                 reference_value = printed$np[[i]], latents = fit$n_parameters,
                 difference = abs(fit$n_parameters - printed$np[[i]]), tolerance = 0, n_values = 1L,
                 kind = "published", agrees = fit$n_parameters == printed$np[[i]]),
      .printed(case, source, "AIC", printed$aic[[i]], fit$aic),
      .printed(case, source, "BIC", printed$bic[[i]], fit$bic))
  }))
}

.permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(.permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

compare_lta_solution <- function(case, package, ours, reference) {
  k <- length(ours$initial)
  perms <- .permutations(k)
  gap <- vapply(perms, function(perm) max(abs(reference$response[perm, , drop = FALSE] - ours$response)),
                numeric(1))
  perm <- perms[[which.min(gap)]]
  rows <- data.frame(
    case = case, reference = package,
    quantity = c("maximized log likelihood", sprintf("initial probabilities (%d values)", k),
                 sprintf("transition probabilities (%d values)", k * k),
                 sprintf("response probabilities (%d values)", length(ours$response))),
    reference_value = c(reference$log_likelihood, NA, NA, NA),
    latents = c(ours$log_likelihood, NA, NA, NA),
    difference = c(abs(reference$log_likelihood - ours$log_likelihood),
                   max(abs(reference$initial[perm] - ours$initial)),
                   max(abs(reference$transition[perm, perm, drop = FALSE] - ours$transition)), min(gap)),
    tolerance = c(1e-6 + 1e-10 * abs(ours$log_likelihood), 1e-4, 1e-4, 1e-4),
    n_values = c(1L, k, k * k, length(ours$response)), kind = "R package")
  rows$agrees <- rows$difference <= rows$tolerance
  rows
}

.printed <- function(case, source, quantity, printed, obtained) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  tolerance <- 0.5 * 10^(-digits) + 1e-9
  value <- as.numeric(printed)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, kind = "published", agrees = abs(obtained - value) <= tolerance)
}

compare_published_lta <- function(fits) {
  source <- "Bartolucci, Pandolfi and Pennoni (2017), JSS 81(4), section 3.1"
  three <- fits$k3
  four <- fits$k4
  initial <- sort(three$initial)
  rbind(
    .printed("rlms_k3", source, "log likelihood", "-13557.21", three$log_likelihood),
    data.frame(case = "rlms_k3", reference = source, quantity = "number of parameters",
               reference_value = 20, latents = three$n_parameters,
               difference = abs(three$n_parameters - 20), tolerance = 0, n_values = 1L,
               kind = "published", agrees = three$n_parameters == 20),
    .printed("rlms_k3", source, "BIC", "27263.39", three$bic),
    .printed("rlms_k3", source, "initial probability (smallest)", "0.1739", initial[[1L]]),
    .printed("rlms_k3", source, "initial probability (middle)", "0.3638", initial[[2L]]),
    .printed("rlms_k3", source, "initial probability (largest)", "0.4623", initial[[3L]]),
    .printed("rlms_k4", source, "log likelihood", "-13392.93", four$log_likelihood),
    .printed("rlms_k4", source, "BIC", "27016.78", four$bic))
}

compare_all_lta_packages <- function() {
  data <- rlms_long()
  fits <- list(k2 = fit_latents_lta(data, 2L), k3 = fit_latents_lta(data, 3L),
               k4 = fit_latents_lta(data, 4L))
  software <- do.call(rbind, lapply(names(fits), function(name) {
    k <- as.integer(sub("k", "", name))
    case <- paste0("rlms_", name)
    rbind(compare_lta_solution(case, "LMest", fits[[name]], fit_lmest(data, k)),
          if (k <= 3L) compare_lta_solution(case, "seqHMM", fits[[name]], fit_seqhmm(data, k)))
  }))
  market <- market_long()
  market_fits <- list(k2 = fit_latents_market(market, 2L), k3 = fit_latents_market(market, 3L),
                      k4 = fit_latents_market(market, 4L))
  market_software <- do.call(rbind, lapply(names(market_fits), function(name) {
    compare_lta_solution(paste0("market_", name), "LMest", market_fits[[name]],
                         fit_lmest_market(market, as.integer(sub("k", "", name))))
  }))
  out <- rbind(software, compare_published_lta(fits), market_software,
               compare_published_market(market_fits))
  out$fn <- "lta"
  out
}
