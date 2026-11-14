# Functions for compare.R: recovery of known parameters from data simulated
# with simulab, for every public function of latents.
#
# Each design draws R datasets from simulab's generators with fixed seeds,
# fits latents to each, matches the fitted classes to the generating ones and
# compares the estimates with the parameters that generated the data:
#
# - "mean estimate": the mean of the R estimates against the true value. It
#   agrees when the difference is within 3.5 Monte Carlo standard errors
#   (3.5 sd / sqrt(R)), a two-sided test at about 0.0005 per parameter, so
#   the whole kit expects no false alarm.
# - "95% interval coverage": the share of the R intervals that contain the
#   true value. It agrees when the share lies in the central 99.9% binomial
#   range for R draws at 0.95.
# - counts that must be exact: replicates fitted, parameters compared.
#
# Where simulab has no generator for a whole design, simulab draws the latent
# classes and their indicators and the kit adds what remains with base R from
# the replicate's own seed: a covariate and a distal outcome that depend on
# the class (r3step(), three_step()), the outcome of a regression mixture
# (mixture_regression()), the cluster classes of a two-level latent class
# model (multilca()), and the class sizes of a growth mixture.
#
# All tolerances are fixed above, before any comparison was run.

# permutations() and delta_terms() from the wider Latent GOLD kit.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

.recovery_reference <- "simulab (known parameters)"

#' One comparison row
#' @noRd
.recovery_row <- function(case, fn, quantity, truth, obtained, tolerance, agrees, n_values) {
  data.frame(case = case, reference = .recovery_reference, quantity = quantity,
             reference_value = truth, latents = obtained, difference = abs(obtained - truth),
             tolerance = tolerance, n_values = as.integer(n_values),
             kind = "simulation recovery", agrees = agrees, fn = fn)
}

#' The rows of one design: bias and coverage per parameter, replicates fitted
#'
#' @param truth Named numeric vector: the generating value of each parameter.
#' @param replicates List, one element per replicate: either a list with
#'   `estimate`, and optionally `lower` and `upper` (named numeric vectors
#'   over `names(truth)`) and its own `truth` where the parameter a label
#'   names depends on the fit (r3step()'s reference class), or a
#'   `recovery_failure`. Bias is then the mean of estimate minus that truth.
#' @noRd
.recovery_rows <- function(case, fn, truth, replicates) {
  fitted <- replicates[!vapply(replicates, inherits, logical(1), "recovery_failure")]
  r_total <- length(replicates)
  r <- length(fitted)
  stopifnot("at least two replicates were fitted" = r >= 2L)
  matrix_of <- function(part) {
    do.call(rbind, lapply(fitted, function(x) {
      values <- x[[part]]
      if (is.null(values)) return(stats::setNames(rep(NA_real_, length(truth)), names(truth)))
      stopifnot("every replicate names every parameter" = all(names(truth) %in% names(values)))
      values[names(truth)]
    }))
  }
  estimates <- matrix_of("estimate")
  lower <- matrix_of("lower")
  upper <- matrix_of("upper")
  truths <- do.call(rbind, lapply(fitted, function(x) (x$truth %||% truth)[names(truth)]))
  bounds <- stats::qbinom(c(0.0005, 0.9995), r, 0.95) / r
  rows <- do.call(rbind, lapply(names(truth), function(name) {
    errors <- estimates[, name] - truths[, name]
    mcse <- stats::sd(errors) / sqrt(r)
    bias <- .recovery_row(case, fn, sprintf("%s: mean estimate", name), mean(truths[, name]),
                          mean(estimates[, name]), 3.5 * mcse, abs(mean(errors)) <= 3.5 * mcse, r)
    if (all(is.na(lower[, name]))) return(bias)
    covered <- mean(lower[, name] <= truths[, name] & truths[, name] <= upper[, name])
    rbind(bias,
          .recovery_row(case, fn, sprintf("%s: 95%% interval coverage", name), 0.95, covered,
                        max(0.95 - bounds[[1L]], bounds[[2L]] - 0.95),
                        covered >= bounds[[1L]] && covered <= bounds[[2L]], r))
  }))
  rbind(.recovery_row(case, fn, "replicates fitted", r_total, r, 0, r == r_total, 1L), rows)
}

#' Run R replicates of a design, each from its own seed
#'
#' A replicate whose fit fails is kept as a `recovery_failure` holding the
#' seed and the condition message: the "replicates fitted" row of its design
#' then disagrees, so a failure is reported, never dropped. (latents' errors
#' share no base class, so every error is caught here and surfaced there.)
#' @param one Function of a seed returning one replicate.
#' @noRd
.recovery_replicates <- function(r, seed, one) {
  lapply(seed + seq_len(r), function(s) {
    tryCatch(one(s), error = function(e) {
      structure(list(seed = s, message = conditionMessage(e)), class = "recovery_failure")
    })
  })
}

#' The order of fitted classes that best matches the generating classes
#' @param fitted,truth Matrices, one row per class (fitted in latents' order,
#'   truth in simulab's), the same class-specific parameters in the columns.
#' @return `perm`: generating class j is fitted class perm[j].
#' @noRd
.recovery_order <- function(fitted, truth) {
  orders <- permutations(nrow(truth))
  gaps <- vapply(orders, function(perm) sum((fitted[perm, , drop = FALSE] - truth)^2),
                 numeric(1))
  orders[[which.min(gaps)]]
}

#' Estimates and interval bounds from parameter_inference() rows
#' @param rows Rows of parameter_inference(); @param labels Their labels.
#' @noRd
.recovery_estimates <- function(rows, labels) {
  stopifnot("one row per label" = nrow(rows) == length(labels))
  list(estimate = stats::setNames(rows$estimate, labels),
       lower = stats::setNames(rows$conf_low, labels),
       upper = stats::setNames(rows$conf_high, labels))
}

#' Join several estimate lists into one
#' @noRd
.recovery_join <- function(...) {
  parts <- list(...)
  list(estimate = unlist(lapply(parts, `[[`, "estimate")),
       lower = unlist(lapply(parts, function(p) p$lower %||%
                               stats::setNames(rep(NA_real_, length(p$estimate)), names(p$estimate)))),
       upper = unlist(lapply(parts, function(p) p$upper %||%
                               stats::setNames(rep(NA_real_, length(p$estimate)), names(p$estimate)))))
}

# ---------------------------------------------------------------------------
# Latent class design: lca(), r3step(), three_step(), enumerate_lca()

.lca_design <- list(
  n = 1000L,
  proportions = c(0.5, 0.3, 0.2),
  # P(category 1) for six binary items in three classes.
  yes = rbind(rep(0.85, 6L), c(0.85, 0.85, 0.85, 0.15, 0.15, 0.15), rep(0.15, 6L)),
  covariate_means = c(0, 1, -1),
  distal_means = c(0, 1.5, 3))

#' One latent class dataset: simulab's classes and items, plus a covariate
#' and a distal outcome normal within each class (variance 1)
#' @noRd
.lca_data <- function(seed, design = .lca_design) {
  k <- nrow(design$yes)
  probabilities <- array(0, c(k, ncol(design$yes), 2L))
  probabilities[, , 1L] <- design$yes
  probabilities[, , 2L] <- 1 - design$yes
  simulated <- simulab::simulate_lca(design$n, probabilities = probabilities,
                                     proportions = design$proportions, seed = seed)
  data <- as.data.frame(simulated)
  items <- grep("^item_", names(data), value = TRUE)
  data[items] <- lapply(data[items], factor, levels = c("1", "2"))
  class <- match(data$latent_class, sprintf("Class %d", seq_len(k)))
  stopifnot("every class label is simulab's" = !anyNA(class))
  set.seed(seed)
  data$x <- stats::rnorm(design$n, design$covariate_means[class])
  data$y <- stats::rnorm(design$n, design$distal_means[class])
  list(data = data, items = items)
}

.lca_truth <- function(design = .lca_design) {
  k <- nrow(design$yes)
  shares <- stats::setNames(design$proportions, sprintf("class %d share", seq_len(k)))
  yes <- stats::setNames(as.vector(t(design$yes)),
                         sprintf("class %d item_%d P(1)", rep(seq_len(k), each = ncol(design$yes)),
                                 rep(seq_len(ncol(design$yes)), k)))
  c(shares, yes)
}

#' Fit lca() to one dataset and align its classes with the generating ones
#' @noRd
.lca_fit <- function(seed, design = .lca_design) {
  simulated <- .lca_data(seed, design)
  fit <- lca(simulated$data, simulated$items, n_classes = nrow(design$yes), n_starts = 10L,
             seed = 1L)
  inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed", boundary = "fix"))
  yes <- inference[inference$level == "measurement" & grepl(":1$", inference$term), ]
  k <- nrow(design$yes)
  fitted <- t(vapply(seq_len(k), function(c)
    yes$estimate[match(sprintf("profile_%d %s:1", c, simulated$items),
                       paste(yes$outcome, yes$term))], numeric(length(simulated$items))))
  perm <- .recovery_order(fitted, design$yes)
  list(fit = fit, data = simulated$data, items = simulated$items, inference = inference,
       perm = perm)
}

.lca_estimates <- function(fitted, design = .lca_design) {
  k <- nrow(design$yes)
  inference <- fitted$inference
  shares <- inference[inference$level == "profile", ]
  share_rows <- shares[match(sprintf("profile_%d", fitted$perm), shares$outcome), ]
  yes <- inference[inference$level == "measurement", ]
  yes_rows <- yes[match(sprintf("profile_%d %s:1", rep(fitted$perm, each = length(fitted$items)),
                                rep(fitted$items, k)), paste(yes$outcome, yes$term)), ]
  truth <- .lca_truth(design)
  .recovery_estimates(rbind(share_rows, yes_rows), names(truth))
}

#' r3step()'s truth: x normal within class c with mean mu_c and variance 1
#' makes the multinomial logit of class c against the reference class r
#' exactly linear, with slope mu_c - mu_r and intercept
#' log(pi_c / pi_r) - (mu_c^2 - mu_r^2) / 2.
#' @noRd
.r3step_truth <- function(reference, design = .lca_design) {
  others <- setdiff(seq_len(nrow(design$yes)), reference)
  mu <- design$covariate_means
  pi <- design$proportions
  stats::setNames(
    c(log(pi[others] / pi[[reference]]) - (mu[others]^2 - mu[[reference]]^2) / 2,
      mu[others] - mu[[reference]]),
    c(sprintf("class %d vs class %d intercept", others, reference),
      sprintf("class %d vs class %d slope", others, reference)))
}

#' Every latent class design function: lca(), r3step(), three_step() with
#' BCH and ML, enumerate_lca()
#' @noRd
compare_recovery_lca <- function(r = 100L, r_enumeration = 30L, seed = 1000L) {
  design <- .lca_design
  fits <- .recovery_replicates(r, seed, function(s) .lca_fit(s, design))
  usable <- fits[!vapply(fits, inherits, logical(1), "recovery_failure")]
  failed <- fits[vapply(fits, inherits, logical(1), "recovery_failure")]
  keep <- function(f) c(lapply(usable, f), failed)
  # r3step(): latents' last fitted class is the reference. Labels name the
  # two other generating classes in ascending order, and each replicate
  # carries the truth for its own reference class.
  k <- nrow(design$yes)
  r3_labels <- c("first other class intercept", "second other class intercept",
                 "first other class slope", "second other class slope")
  r3 <- keep(function(fitted) {
    reference <- match(k, fitted$perm)
    others <- setdiff(seq_len(k), reference)
    table <- as.data.frame(r3step(fitted$fit, fitted$data, covariates = "x"))
    rows <- table[match(c(sprintf("profile_%d (Intercept)", fitted$perm[others]),
                          sprintf("profile_%d x", fitted$perm[others])),
                        paste(table$outcome, table$term)), ]
    estimates <- .recovery_estimates(rows, r3_labels)
    estimates$truth <- stats::setNames(.r3step_truth(reference, design), r3_labels)
    estimates
  })
  r3_truth <- stats::setNames(.r3step_truth(k, design), r3_labels)
  distal_truth <- stats::setNames(design$distal_means, sprintf("class %d mean", seq_len(k)))
  distal <- function(method) keep(function(fitted) {
    table <- as.data.frame(three_step(fitted$fit, fitted$data, outcome = "y", method = method))
    rows <- table[match(fitted$perm, table$class), ]
    .recovery_estimates(rows, names(distal_truth))
  })
  enumeration <- .recovery_replicates(r_enumeration, seed + 5000L, function(s) {
    simulated <- .lca_data(s, design)
    grid <- enumerate_lca(simulated$data[simulated$items], simulated$items, n_classes = 1:5,
                          n_starts = 10L, seed = 1L)
    best <- as.data.frame(get_results(grid, "best"))
    list(estimate = c(`BIC selects 3 classes` = as.numeric(
      best$n_profiles[best$criterion == "bic" & best$convention == "individuals"] == k)))
  })
  rbind(
    .recovery_rows("lca_k3", "lca", .lca_truth(design), keep(.lca_estimates)),
    .recovery_rows("lca_r3step", "three_step", r3_truth, r3),
    .recovery_rows("lca_three_step_bch", "three_step", distal_truth, distal("bch")),
    .recovery_rows("lca_three_step_ml", "three_step", distal_truth, distal("ml")),
    .recovery_selection("lca_enumeration", enumeration))
}

#' The share of replicates in which a criterion selects the generating
#' number of classes: it agrees at 0.9 or more (designs with well separated
#' classes, chosen so; fixed before the run)
#' @noRd
.recovery_selection <- function(case, replicates) {
  fitted <- replicates[!vapply(replicates, inherits, logical(1), "recovery_failure")]
  share <- mean(vapply(fitted, function(x) x$estimate[[1L]], numeric(1)))
  quantity <- names(fitted[[1L]]$estimate)[[1L]]
  rbind(.recovery_row(case, "enumeration", "replicates fitted", length(replicates),
                      length(fitted), 0, length(fitted) == length(replicates), 1L),
        .recovery_row(case, "enumeration", sprintf("%s (share of replicates)", quantity), 1,
                      share, 0.1, share >= 0.9, length(fitted)))
}

# ---------------------------------------------------------------------------
# Latent profile design: lpa(), enumerate_lpa()

.lpa_design <- list(
  n = 800L,
  proportions = c(0.4, 0.35, 0.25),
  means = rbind(c(0, 0, 0, 0), c(2, 2, 0, 0), c(0, 2, 2, 2)),
  sd = 1)

.lpa_data <- function(seed, design = .lpa_design) {
  simulated <- simulab::simulate_lpa(design$n, means = design$means, sds = design$sd,
                                     proportions = design$proportions, seed = seed)
  data <- as.data.frame(simulated)
  list(data = data, vars = grep("^indicator_", names(data), value = TRUE))
}

.lpa_truth <- function(design = .lpa_design) {
  k <- nrow(design$means)
  p <- ncol(design$means)
  c(stats::setNames(design$proportions, sprintf("profile %d share", seq_len(k))),
    stats::setNames(as.vector(t(design$means)),
                    sprintf("profile %d indicator_%d mean", rep(seq_len(k), each = p),
                            rep(seq_len(p), k))),
    stats::setNames(rep(design$sd^2, p), sprintf("indicator_%d variance", seq_len(p))))
}

compare_recovery_lpa <- function(r = 100L, r_enumeration = 30L, seed = 2000L) {
  design <- .lpa_design
  k <- nrow(design$means)
  truth <- .lpa_truth(design)
  replicates <- .recovery_replicates(r, seed, function(s) {
    simulated <- .lpa_data(s, design)
    fit <- lpa(simulated$data, simulated$vars, n_profiles = k, model = "EEI", n_starts = 10L,
               seed = 1L)
    inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed"))
    means <- inference[inference$level == "measurement" & inference$parameter == "mean", ]
    fitted <- t(vapply(seq_len(k), function(c)
      means$estimate[match(sprintf("profile_%d %s", c, simulated$vars),
                           paste(means$outcome, means$term))], numeric(length(simulated$vars))))
    perm <- .recovery_order(fitted, design$means)
    shares <- inference[inference$level == "profile", ]
    variances <- inference[inference$level == "measurement" & inference$parameter == "variance", ]
    rows <- rbind(
      shares[match(sprintf("profile_%d", perm), shares$outcome), ],
      means[match(sprintf("profile_%d %s", rep(perm, each = length(simulated$vars)),
                          rep(simulated$vars, k)), paste(means$outcome, means$term)), ],
      variances[match(simulated$vars, variances$term), ])
    .recovery_estimates(rows, names(truth))
  })
  enumeration <- .recovery_replicates(r_enumeration, seed + 5000L, function(s) {
    simulated <- .lpa_data(s, design)
    grid <- enumerate_lpa(simulated$data, simulated$vars, n_profiles = 1:5, model = "EEI",
                          n_starts = 10L, seed = 1L)
    best <- as.data.frame(get_results(grid, "best"))
    list(estimate = c(`BIC selects 3 profiles` = as.numeric(
      best$n_profiles[best$criterion == "bic" & best$convention == "individuals"] == k)))
  })
  rbind(.recovery_rows("lpa_eei_k3", "lpa", truth, replicates),
        .recovery_selection("lpa_enumeration", enumeration))
}

# ---------------------------------------------------------------------------
# Latent transition design: lta() on simulab's hidden Markov chains

# Every probability at 0.15 or more, and 1000 chains: a maximum with a
# response probability on its bound is then rare. (Transition fits accept
# parameter_inference(boundary = "fix") but do not apply it, so such a
# replicate is refused and counted as not fitted.)
.lta_design <- list(
  n = 1000L, occasions = 5L,
  initial = c(0.6, 0.4),
  transition = rbind(c(0.85, 0.15), c(0.2, 0.8)),
  emission = rbind(c(0.65, 0.2, 0.15), c(0.15, 0.2, 0.65)))

.lta_truth <- function(design = .lta_design) {
  k <- nrow(design$transition)
  m <- ncol(design$emission)
  c(stats::setNames(design$initial, sprintf("state %d initial", seq_len(k))),
    stats::setNames(as.vector(t(design$transition)),
                    sprintf("state %d to state %d", rep(seq_len(k), each = k), rep(seq_len(k), k))),
    stats::setNames(as.vector(t(design$emission)),
                    sprintf("state %d observation %d", rep(seq_len(k), each = m),
                            rep(seq_len(m), k))))
}

compare_recovery_lta <- function(r = 100L, seed = 3000L) {
  design <- .lta_design
  k <- nrow(design$transition)
  m <- ncol(design$emission)
  truth <- .lta_truth(design)
  replicates <- .recovery_replicates(r, seed, function(s) {
    simulated <- simulab::simulate_hmm(design$n, transition = design$transition,
                                       chain_length = design$occasions, emission = design$emission,
                                       initial = design$initial, seed = s)
    data <- as.data.frame(simulated)
    data$observation <- factor(data$observation, levels = sprintf("Observation %d", seq_len(m)))
    fit <- lta(data, "observation", id = "id", n_profiles = k, time = "occasion",
               categorical = "observation", n_starts = 10L, seed = 1L)
    inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed", boundary = "fix"))
    responses <- inference[inference$level == "measurement", ]
    term <- sprintf("observation:Observation %d", seq_len(m))
    fitted <- t(vapply(seq_len(k), function(c)
      responses$estimate[match(sprintf("profile_%d %s", c, term),
                               paste(responses$outcome, responses$term))], numeric(m)))
    perm <- .recovery_order(fitted, design$emission)
    initial <- inference[inference$parameter == "initial_probability", ]
    transitions <- inference[inference$parameter == "transition_probability", ]
    rows <- rbind(
      initial[match(sprintf("profile_%d", perm), initial$outcome), ],
      # outcome: the state moved to; term: group class and the state moved from.
      transitions[match(sprintf("profile_%d group_class_1:profile_%d", rep(perm, k),
                                rep(perm, each = k)),
                        paste(transitions$outcome, transitions$term)), ],
      responses[match(sprintf("profile_%d %s", rep(perm, each = m), rep(term, k)),
                      paste(responses$outcome, responses$term)), ])
    .recovery_estimates(rows, names(truth))
  })
  .recovery_rows("lta_hmm_k2", "lta", truth, replicates)
}

# ---------------------------------------------------------------------------
# Regression mixture design: mixture_regression()

.regression_design <- list(
  n = 600L, proportions = c(0.6, 0.4),
  intercepts = c(0, 2), slopes = c(1, -1), sigma = 0.5)

compare_recovery_regression <- function(r = 100L, seed = 4000L) {
  design <- .regression_design
  truth <- c(`class 1 (Intercept)` = design$intercepts[[1L]], `class 1 x` = design$slopes[[1L]],
             `class 2 (Intercept)` = design$intercepts[[2L]], `class 2 x` = design$slopes[[2L]],
             `residual sd` = design$sigma,
             `class 2 membership logit` = log(design$proportions[[2L]] / design$proportions[[1L]]))
  replicates <- .recovery_replicates(r, seed, function(s) {
    # simulab draws each unit's class and a covariate independent of it; the
    # outcome is the class's regression line plus normal error.
    simulated <- as.data.frame(simulab::simulate_lpa(design$n, means = rbind(0, 0), sds = 1,
                                                     proportions = design$proportions,
                                                     seed = s))
    class <- match(simulated$profile, c("Profile 1", "Profile 2"))
    data <- data.frame(x = simulated[[grep("^indicator_", names(simulated))[[1L]]]])
    set.seed(s)
    data$y <- design$intercepts[class] + design$slopes[class] * data$x +
      stats::rnorm(design$n, sd = design$sigma)
    fit <- mixture_regression(y ~ x, data, n_classes = 2L, variance = "equal", n_starts = 10L,
                              seed = 1L)
    theta <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    perm <- .recovery_order(
      rbind(theta[c("coefficient.class_1.(Intercept)", "coefficient.class_1.x")],
            theta[c("coefficient.class_2.(Intercept)", "coefficient.class_2.x")]),
      cbind(design$intercepts, design$slopes))
    names_of <- function(c) sprintf(c("coefficient.class_%d.(Intercept)", "coefficient.class_%d.x"), c)
    map <- function(values) {
      membership <- values[["membership.class_2.(Intercept)"]]
      c(values[names_of(perm[[1L]])], values[names_of(perm[[2L]])],
        exp(values[["log_sigma.all"]]),
        if (perm[[2L]] == 2L) membership else -membership)
    }
    terms <- delta_terms(function(values) stats::setNames(map(values), names(truth)), theta,
                         covariance)
    z <- stats::qnorm(0.975)
    list(estimate = stats::setNames(terms$estimate, names(truth)),
         lower = stats::setNames(terms$estimate - z * terms$standard_error, names(truth)),
         upper = stats::setNames(terms$estimate + z * terms$standard_error, names(truth)))
  })
  .recovery_rows("mixture_regression_k2", "mixture_regression", truth, replicates)
}

# ---------------------------------------------------------------------------
# Growth mixture design: mixture_regression(random =)

.growth_design <- list(
  n = 500L, share_1 = 0.6, times = 0:4,
  intercepts = c(0, 3), slopes = c(1, -0.5),
  random_sd = c(0.6, 0.2), residual_sd = 0.5)

compare_recovery_growth <- function(r = 50L, seed = 5000L) {
  design <- .growth_design
  truth <- c(`class 1 (Intercept)` = design$intercepts[[1L]], `class 1 time` = design$slopes[[1L]],
             `class 2 (Intercept)` = design$intercepts[[2L]], `class 2 time` = design$slopes[[2L]],
             `residual sd` = design$residual_sd,
             `random intercept sd` = design$random_sd[[1L]],
             `random slope sd` = design$random_sd[[2L]],
             `class 2 membership logit` = log((1 - design$share_1) / design$share_1))
  replicates <- .recovery_replicates(r, seed, function(s) {
    # The class sizes are binomial; simulab draws each class's trajectories.
    set.seed(s)
    n_1 <- stats::rbinom(1L, design$n, design$share_1)
    draw <- function(n, c, seed) as.data.frame(simulab::simulate_growth(
      n, times = design$times, intercept = design$intercepts[[c]], slope = design$slopes[[c]],
      random_sd = design$random_sd, residual_sd = design$residual_sd, seed = seed))
    first <- draw(n_1, 1L, s)
    second <- draw(design$n - n_1, 2L, s + 100000L)
    second$id <- second$id + n_1
    data <- rbind(first, second)
    data <- data[order(data$id, data$time), ]
    fit <- mixture_regression(outcome ~ time, data, n_classes = 2L, id = "id",
                              class_level = "group", random = ~ 1 + time,
                              random_covariance = "equal", variance = "equal", n_starts = 10L,
                              seed = 1L)
    theta <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    names_of <- function(c) sprintf(c("coefficient.class_%d.(Intercept)", "coefficient.class_%d.time"), c)
    perm <- .recovery_order(rbind(theta[names_of(1L)], theta[names_of(2L)]),
                            cbind(design$intercepts, design$slopes))
    map <- function(values) {
      l11 <- exp(values[["random.log_cholesky.shared.(Intercept):(Intercept)"]])
      l22 <- exp(values[["random.log_cholesky.shared.time:time"]])
      l21 <- values[["random.cholesky.shared.time:(Intercept)"]]
      membership <- values[["membership.class_2.(Intercept)"]]
      c(values[names_of(perm[[1L]])], values[names_of(perm[[2L]])],
        exp(values[["log_sigma.shared"]]), l11, sqrt(l21^2 + l22^2),
        if (perm[[2L]] == 2L) membership else -membership)
    }
    terms <- delta_terms(function(values) stats::setNames(map(values), names(truth)), theta,
                         covariance)
    z <- stats::qnorm(0.975)
    list(estimate = stats::setNames(terms$estimate, names(truth)),
         lower = stats::setNames(terms$estimate - z * terms$standard_error, names(truth)),
         upper = stats::setNames(terms$estimate + z * terms$standard_error, names(truth)))
  })
  .recovery_rows("growth_k2", "growth", truth, replicates)
}

# ---------------------------------------------------------------------------
# Two-level designs: multilpa(), multilca()

.multilpa_design <- list(
  clusters = 60L, cluster_size = 20L,
  means = rbind(c(0, 0, 0), c(2, 2, 0), c(0, 2, 2)),
  profile_probabilities = rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6)),
  cluster_class_proportions = c(0.5, 0.5))

#' The two-level truth: shares within each group class, the group class
#' shares, and the measurement parameters
#' @noRd
.two_level_truth <- function(profile_probabilities, group_shares, measurement) {
  g <- nrow(profile_probabilities)
  k <- ncol(profile_probabilities)
  c(stats::setNames(group_shares, sprintf("group class %d share", seq_len(g))),
    stats::setNames(as.vector(t(profile_probabilities)),
                    sprintf("group class %d profile %d share", rep(seq_len(g), each = k),
                            rep(seq_len(k), g))),
    measurement)
}

#' Rows for group class and within-group-class shares, given the profile and
#' group class orders
#' @noRd
.two_level_shares <- function(inference, perm, group_perm) {
  groups <- inference[inference$level == "group", ]
  profiles <- inference[inference$level == "profile", ]
  k <- length(perm)
  rbind(groups[match(sprintf("group_class_%d", group_perm), groups$outcome), ],
        profiles[match(sprintf("profile_%d group_class_%d", rep(perm, length(group_perm)),
                               rep(group_perm, each = k)),
                       paste(profiles$outcome, profiles$term)), ])
}

#' The group class order: matched on the within-group-class profile shares
#' (already in generating profile order)
#' @noRd
.group_order <- function(inference, perm, truth) {
  profiles <- inference[inference$level == "profile", ]
  g <- nrow(truth)
  fitted <- t(vapply(seq_len(g), function(gc)
    profiles$estimate[match(sprintf("profile_%d group_class_%d", perm, gc),
                            paste(profiles$outcome, profiles$term))], numeric(length(perm))))
  .recovery_order(fitted, truth)
}

compare_recovery_multilpa <- function(r = 40L, seed = 6000L) {
  design <- .multilpa_design
  k <- nrow(design$means)
  p <- ncol(design$means)
  truth <- .two_level_truth(
    design$profile_probabilities, design$cluster_class_proportions,
    c(stats::setNames(as.vector(t(design$means)),
                      sprintf("profile %d indicator_%d mean", rep(seq_len(k), each = p),
                              rep(seq_len(p), k))),
      stats::setNames(rep(1, p), sprintf("indicator_%d variance", seq_len(p)))))
  replicates <- .recovery_replicates(r, seed, function(s) {
    data <- as.data.frame(simulab::simulate_ml_lpa(
      clusters = design$clusters, cluster_size = design$cluster_size, means = design$means,
      profile_probabilities = design$profile_probabilities,
      cluster_class_proportions = design$cluster_class_proportions, seed = s))
    vars <- grep("^indicator_", names(data), value = TRUE)
    fit <- multilpa(data, vars, id = "cluster", n_profiles = k, n_group_classes = 2L,
                    variance_model = "equal", n_starts = 10L, seed = 1L)
    inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed"))
    means <- inference[inference$level == "measurement" & inference$parameter == "mean", ]
    fitted <- t(vapply(seq_len(k), function(c)
      means$estimate[match(sprintf("profile_%d %s", c, vars), paste(means$outcome, means$term))],
      numeric(p)))
    perm <- .recovery_order(fitted, design$means)
    group_perm <- .group_order(inference, perm, design$profile_probabilities)
    variances <- inference[inference$level == "measurement" & inference$parameter == "variance", ]
    rows <- rbind(
      .two_level_shares(inference, perm, group_perm),
      means[match(sprintf("profile_%d %s", rep(perm, each = p), rep(vars, k)),
                  paste(means$outcome, means$term)), ],
      variances[match(vars, variances$term), ])
    .recovery_estimates(rows, names(truth))
  })
  .recovery_rows("multilpa_k3_g2", "multilpa", truth, replicates)
}

.multilca_design <- list(
  clusters = 60L, cluster_size = 25L,
  yes = rbind(rep(0.85, 6L), c(0.85, 0.85, 0.85, 0.15, 0.15, 0.15), rep(0.15, 6L)),
  profile_probabilities = rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6)),
  cluster_class_proportions = c(0.5, 0.5))

compare_recovery_multilca <- function(r = 40L, seed = 7000L) {
  design <- .multilca_design
  k <- nrow(design$yes)
  items_n <- ncol(design$yes)
  truth <- .two_level_truth(
    design$profile_probabilities, design$cluster_class_proportions,
    stats::setNames(as.vector(t(design$yes)),
                    sprintf("class %d item_%d P(1)", rep(seq_len(k), each = items_n),
                            rep(seq_len(items_n), k))))
  probabilities <- array(0, c(k, items_n, 2L))
  probabilities[, , 1L] <- design$yes
  probabilities[, , 2L] <- 1 - design$yes
  replicates <- .recovery_replicates(r, seed, function(s) {
    # Each cluster's class is drawn with the replicate's seed; simulab draws
    # the members' classes and items from that cluster class's shares.
    set.seed(s)
    cluster_class <- sample.int(2L, design$clusters, replace = TRUE,
                                prob = design$cluster_class_proportions)
    data <- do.call(rbind, lapply(seq_len(design$clusters), function(j) {
      members <- as.data.frame(simulab::simulate_lca(
        design$cluster_size, probabilities = probabilities,
        proportions = design$profile_probabilities[cluster_class[[j]], ], seed = s * 1000L + j))
      members$cluster <- j
      members
    }))
    items <- grep("^item_", names(data), value = TRUE)
    data[items] <- lapply(data[items], factor, levels = c("1", "2"))
    fit <- multilca(data, items, id = "cluster", n_profiles = k, n_group_classes = 2L,
                    n_starts = 10L, seed = 1L)
    inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed", boundary = "fix"))
    yes <- inference[inference$level == "measurement" & grepl(":1$", inference$term), ]
    fitted <- t(vapply(seq_len(k), function(c)
      yes$estimate[match(sprintf("profile_%d %s:1", c, items), paste(yes$outcome, yes$term))],
      numeric(items_n)))
    perm <- .recovery_order(fitted, design$yes)
    group_perm <- .group_order(inference, perm, design$profile_probabilities)
    rows <- rbind(
      .two_level_shares(inference, perm, group_perm),
      yes[match(sprintf("profile_%d %s:1", rep(perm, each = items_n), rep(items, k)),
                paste(yes$outcome, yes$term)), ])
    .recovery_estimates(rows, names(truth))
  })
  .recovery_rows("multilca_k3_g2", "multilca", truth, replicates)
}

# ---------------------------------------------------------------------------
# Bootstrap likelihood-ratio test: size and power

#' The share of replicates in which bootstrap_lrt() rejects 2 against 3
#' classes at p <= 0.05 (19 bootstrap draws: p = 1/20 when the observed
#' statistic exceeds every draw). Under two generating classes it agrees when
#' the share lies in the central 99.9% binomial range at 0.05; under three
#' well separated classes when it is 0.9 or more.
#' @noRd
compare_recovery_bootstrap <- function(r_size = 20L, r_power = 10L, seed = 8000L) {
  design <- .lca_design
  test <- function(s, k) {
    yes <- design$yes[seq_len(k), , drop = FALSE]
    probabilities <- array(0, c(k, ncol(yes), 2L))
    probabilities[, , 1L] <- yes
    probabilities[, , 2L] <- 1 - yes
    data <- as.data.frame(simulab::simulate_lca(500L, probabilities = probabilities,
                                                proportions = design$proportions[seq_len(k)],
                                                seed = s))
    items <- grep("^item_", names(data), value = TRUE)
    data[items] <- lapply(data[items], factor, levels = c("1", "2"))
    null <- lca(data, items, n_classes = 2L, n_starts = 10L, seed = 1L)
    alternative <- lca(data, items, n_classes = 3L, n_starts = 10L, seed = 1L)
    result <- as.data.frame(bootstrap_lrt(null, alternative, iter = 19L, seed = 1L))
    list(estimate = c(reject = as.numeric(result$p_value <= 0.05)))
  }
  share <- function(replicates) {
    fitted <- replicates[!vapply(replicates, inherits, logical(1), "recovery_failure")]
    list(n = length(fitted), share = mean(vapply(fitted, function(x) x$estimate[["reject"]],
                                                 numeric(1))))
  }
  size <- share(.recovery_replicates(r_size, seed, function(s) test(s, 2L)))
  power <- share(.recovery_replicates(r_power, seed + 500L, function(s) test(s, 3L)))
  bounds <- stats::qbinom(c(0.0005, 0.9995), size$n, 0.05) / size$n
  rbind(
    .recovery_row("bootstrap_lrt_size", "bootstrap_lrt", "replicates fitted", r_size, size$n, 0,
                  size$n == r_size, 1L),
    .recovery_row("bootstrap_lrt_size", "bootstrap_lrt",
                  "rejection rate at 0.05, two generating classes", 0.05, size$share,
                  max(0.05 - bounds[[1L]], bounds[[2L]] - 0.05),
                  size$share >= bounds[[1L]] && size$share <= bounds[[2L]], size$n),
    .recovery_row("bootstrap_lrt_power", "bootstrap_lrt", "replicates fitted", r_power, power$n, 0,
                  power$n == r_power, 1L),
    .recovery_row("bootstrap_lrt_power", "bootstrap_lrt",
                  "rejection rate at 0.05, three generating classes", 1, power$share, 0.1,
                  power$share >= 0.9, power$n))
}

# ---------------------------------------------------------------------------
# Multiple imputation: pool_imputations()

.pool_design <- list(n = 500L, proportions = c(0.6, 0.4),
                     means = rbind(c(0, 0, 0), c(2.5, 2, 1.5)), missing = 0.15, m = 5L)

compare_recovery_pool <- function(r = 30L, seed = 9000L) {
  stopifnot("mice is needed for the imputation design" = requireNamespace("mice", quietly = TRUE))
  design <- .pool_design
  k <- nrow(design$means)
  p <- ncol(design$means)
  truth <- c(stats::setNames(design$proportions, sprintf("profile %d share", seq_len(k))),
             stats::setNames(as.vector(t(design$means)),
                             sprintf("profile %d indicator_%d mean", rep(seq_len(k), each = p),
                                     rep(seq_len(p), k))))
  replicates <- .recovery_replicates(r, seed, function(s) {
    simulated <- as.data.frame(simulab::simulate_lpa(design$n, means = design$means, sds = 1,
                                                     proportions = design$proportions, seed = s))
    vars <- grep("^indicator_", names(simulated), value = TRUE)
    incomplete <- as.data.frame(simulab::inject_missingness(
      simulated[c("id", vars)], mechanism = "MCAR", proportion = design$missing, variables = vars,
      seed = s))
    # The imputation model includes the generating profile, so the
    # imputations are proper and the comparison isolates the pooling.
    incomplete$profile <- factor(simulated$profile)
    predictors <- mice::make.predictorMatrix(incomplete)
    predictors[, "id"] <- 0
    imputed <- mice::mice(incomplete, m = design$m, method = "norm", predictorMatrix = predictors,
                          seed = s, printFlag = FALSE)
    completed <- lapply(seq_len(design$m), function(i) mice::complete(imputed, i)[c("id", vars)])
    pooled <- as.data.frame(pool_imputations(completed, vars, id = "id", n_profiles = k,
                                             n_group_classes = 1L, variance_model = "equal",
                                             n_starts = 10L, seed = 1L))
    means <- pooled[pooled$level == "measurement" & pooled$parameter == "mean", ]
    fitted <- t(vapply(seq_len(k), function(c)
      means$estimate[match(sprintf("profile_%d %s", c, vars), paste(means$outcome, means$term))],
      numeric(p)))
    perm <- .recovery_order(fitted, design$means)
    shares <- pooled[pooled$level == "profile", ]
    rows <- rbind(shares[match(sprintf("profile_%d", perm), shares$outcome), ],
                  means[match(sprintf("profile_%d %s", rep(perm, each = p), rep(vars, k)),
                              paste(means$outcome, means$term)), ])
    .recovery_estimates(rows, names(truth))
  })
  .recovery_rows("pool_imputations_lpa_k2", "pool_imputations", truth, replicates)
}

compare_all_recovery <- function() {
  stopifnot("simulab is needed to draw the datasets" = requireNamespace("simulab", quietly = TRUE))
  rbind(compare_recovery_lca(), compare_recovery_lpa(), compare_recovery_lta(),
        compare_recovery_regression(), compare_recovery_growth(), compare_recovery_multilpa(),
        compare_recovery_multilca(), compare_recovery_bootstrap(), compare_recovery_pool())
}
