# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Accessors and helpers outside the fitting functions: pool_imputations() and the
# latents_pooled methods; diagnostics() with its tables, get_replication(),
# get_checks(); starting_values() and the multilpa_start methods; sensitivity();
# align_labels(); fit_staged(). References: Rubin's rules from mice::pool.scalar()
# (per parameter) and mitools::MIcombine() (the whole covariance matrix); mclust's
# E-step at latents' estimates and its EM continued to the maximum; the
# classification summaries, the replication tolerance and a dense mixture
# likelihood from their definitions; a two-level likelihood with the measurement
# held, maximised by stats::optim().

# --- pool_imputations() ---------------------------------------------------------

# Two well separated profiles; three stand-in imputations that redraw the first
# fifteen values of y1, which is all the pooling step can see of an imputation model.
pool_frame <- local({
  set.seed(5)
  z <- stats::rnorm(200)
  state <- 2L - stats::rbinom(200, 1L, stats::plogis(1.2 * z))
  data.frame(y1 = stats::rnorm(200, c(0, 2.5)[state]), y2 = stats::rnorm(200, c(0, 2)[state]))
})
pool_imputed <- local({
  set.seed(1)
  lapply(1:3, \(index) transform(pool_frame, y1 = c(stats::rnorm(15, 1), y1[-(1:15)])))
})
pooled <- pool_imputations(pool_imputed, c("y1", "y2"), id = NULL, n_profiles = 2L,
                           n_starts = 2, seed = 1)
pool_keys <- c("level", "outcome", "term", "parameter")

# mice::pool.scalar() on each parameter's per-imputation estimate and variance. The
# group-class probability of a single-level model is fixed at one (variance zero),
# where Rubin's degrees of freedom are 0 / 0: it is left out of the comparisons.
pool_mice <- function() {
  per_imputation <- subset(get_results(pooled, "imputations"), level != "group")
  do.call(rbind, lapply(split(per_imputation, per_imputation[pool_keys], drop = TRUE), \(rows) {
    scalar <- mice::pool.scalar(rows$estimate, rows$standard_error^2, n = Inf)
    data.frame(rows[1L, pool_keys], estimate = scalar$qbar, standard_error = sqrt(scalar$t),
               df = scalar$df, fmi = scalar$fmi, row.names = NULL)
  }))
}

# coef() names a parameter level.parameter.outcome.term.
pool_names <- function(table) paste(table$level, table$parameter, table$outcome, table$term, sep = ".")

test_that("the pooled table, its summary and its intervals follow Rubin's rules (mice)", {
  skip_without_reference("mice")
  reference <- pool_mice()
  # Rubin's rules are a few sums of the same doubles: rounding only.
  expect_reference(subset(get_results(pooled), level != "group"), reference, reference = "mice",
                   covers = c("pool_imputations()", "get_results.latents_pooled()"),
                   quantity = "pooled estimates, standard errors, df and fmi", tolerance = 1e-10,
                   measure = "relative", key = pool_keys)
  # Wald intervals with Rubin's degrees of freedom: estimate -/+ t(0.975, df) times the total standard error.
  free <- pool_names(reference)
  half_width <- stats::qt(0.975, reference$df) * reference$standard_error
  expect_reference(confint(pooled)[free, ],
                   cbind(`2.5%` = reference$estimate - half_width, `97.5%` = reference$estimate + half_width),
                   reference = "mice", covers = c("pool_imputations()", "confint.latents_pooled()"),
                   quantity = "95% intervals with Rubin's degrees of freedom", tolerance = 1e-10,
                   measure = "relative", key = c(paste(free, "low"), paste(free, "high")))
  # The fixed group-class row has a fraction of missing information of zero in
  # latents and undefined in mice; zero does not move the maximum.
  expect_reference(cbind(summary(pooled), key = "pooled"),
                   data.frame(key = "pooled", fmi_max = max(reference$fmi)), reference = "mice",
                   covers = c("pool_imputations()", "summary.latents_pooled()"),
                   quantity = "largest fraction of missing information", tolerance = 1e-10,
                   measure = "relative", key = "key")
})

test_that("the pooled estimates and covariance matrix equal mitools::MIcombine()", {
  skip_without_reference("mitools")
  # The imputations' own fits, numbered like the first, as pool_imputations() aligns them.
  fits <- lapply(pool_imputed, \(data) lpa(data, c("y1", "y2"), n_profiles = 2, n_starts = 2, seed = 1))
  fits <- c(fits[1L], lapply(fits[-1L], \(fit) align_labels(fit, fits[[1L]])))
  # Unit check: these are the fits pool_imputations() pooled.
  expect_equal(unname(unlist(lapply(fits, coef))), get_results(pooled, "imputations")$estimate,
               tolerance = 1e-12)
  combined <- mitools::MIcombine(lapply(fits, coef), lapply(fits, vcov))
  # Means of three doubles: rounding only.
  expect_reference(coef(pooled), stats::coef(combined), reference = "mitools",
                   covers = c("pool_imputations()", "coef.latents_pooled()"),
                   quantity = "pooled estimates", tolerance = 1e-12, measure = "relative")
  # U + (1 + 1/m) B: entries of order 1e-2 and below, the fixed row and column
  # exactly zero, so the difference is absolute.
  expect_reference(vcov(pooled), stats::vcov(combined), reference = "mitools",
                   covers = c("pool_imputations()", "vcov.latents_pooled()"),
                   quantity = "pooled covariance matrix", tolerance = 1e-12, measure = "absolute")
})

test_that("the summary's likelihood range is a dense mixture likelihood at each imputation's estimates", {
  per_imputation <- get_results(pooled, "imputations")
  # Two diagonal normal profiles, from stats::dnorm(), at one imputation's estimates.
  dense <- function(index) {
    rows <- subset(per_imputation, imputation == index)
    value <- \(parameter, k, term) rows$estimate[rows$parameter == parameter &
                                                    rows$outcome == sprintf("profile_%d", k) & rows$term == term]
    data <- pool_imputed[[index]]
    density <- vapply(1:2, \(k) value("probability", k, "group_class_1") *
                        stats::dnorm(data$y1, value("mean", k, "y1"), sqrt(value("variance", k, "y1"))) *
                        stats::dnorm(data$y2, value("mean", k, "y2"), sqrt(value("variance", k, "y2"))),
                      numeric(nrow(data)))
    sum(log(rowSums(density)))
  }
  likelihoods <- vapply(seq_along(pool_imputed), dense, numeric(1))
  # Sums of 200 logs, computed twice: rounding only.
  expect_reference(cbind(summary(pooled), key = "pooled"),
                   data.frame(key = "pooled", m = length(pool_imputed),
                              log_likelihood_min = min(likelihoods), log_likelihood_max = max(likelihoods)),
                   reference = "independent: two-profile normal mixture likelihood from stats::dnorm()",
                   covers = c("pool_imputations()", "summary.latents_pooled()"),
                   quantity = "number of imputations and range of their log likelihoods", tolerance = 1e-10,
                   measure = "relative", key = "key")
})

# --- diagnostics(), get_replication(), get_checks() -----------------------------

faithful_x <- as.matrix(datasets::faithful)
faithful_n <- nrow(faithful_x)
diagnosed_fit <- lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 2, n_starts = 5,
                     seed = 1, tol = 1e-12)
diagnosed <- diagnostics(diagnosed_fit, plots = FALSE)

# mclust's E-step (its VVV model with diagonal covariances) at latents' estimates.
faithful_estep <- function(fit) {
  profiles <- get_results(fit, "profiles")
  means <- matrix(profiles$mean, 2L, 2L)
  sigma <- array(unlist(lapply(1:2, \(k) diag(subset(profiles, profile == k)$variance))), c(2L, 2L, 2L))
  shares <- get_results(fit, "profile_probabilities")$probability
  parameters <- list(pro = shares, mean = means,
                     variance = list(modelName = "VVV", d = 2L, G = 2L, sigma = sigma,
                                     cholsigma = array(apply(sigma, 3L, chol), c(2L, 2L, 2L))))
  mclust::estepVVV(data = faithful_x, parameters = parameters)
}

# mclust's EM continued from those posteriors to its maximum, as a profiles table.
faithful_maximum <- function(fit) {
  continued <- mclust::meVVI(data = faithful_x, z = faithful_estep(fit)$z,
                             control = mclust::emControl(tol = c(1e-12, 1e-12)))
  variances <- as.vector(apply(continued$parameters$variance$sigma, 3L, diag))
  list(log_likelihood = continued$loglik,
       profiles = data.frame(profile = rep(1:2, each = 2), indicator = rep(colnames(faithful_x), 2),
                             mean = as.vector(continued$parameters$mean), variance = variances,
                             standard_deviation = sqrt(variances)))
}

test_that("diagnostics() tables follow their definitions from mclust's posteriors", {
  skip_without_reference("mclust")
  z <- faithful_estep(diagnosed_fit)$z
  entropy_sum <- -sum(ifelse(z > 0, z * log(z), 0))
  modal <- max.col(z, ties.method = "first")
  shares <- colMeans(z)
  diagonal <- vapply(1:2, \(k) mean(z[modal == k, k]), numeric(1))
  definitions <- "independent: classification summaries from their definitions, mclust's E-step posteriors"
  # Ramaswamy, DeSarbo, Reibstein and Robinson (1993): 1 - EN / (n log K). Sums of
  # posteriors that agree with latents' to rounding: 1e-10 relative.
  entropy <- data.frame(level = "individuals", entropy_sum = entropy_sum,
                        relative_entropy = 1 - entropy_sum / (faithful_n * log(2)))
  expect_reference(subset(as.data.frame(diagnosed), level == "individuals"), entropy,
                   reference = definitions, covers = c("diagnostics()", "as.data.frame.multilpa_diagnostics()"),
                   quantity = "entropy", tolerance = 1e-10, measure = "relative", key = "level")
  expect_reference(subset(as.data.frame(summary(diagnosed)), level == "individuals"), entropy,
                   reference = definitions,
                   covers = c("summary.multilpa_diagnostics()", "as.data.frame.multilpa_diagnostics()"),
                   quantity = "entropy, summary()", tolerance = 1e-10, measure = "relative", key = "level")
  # Nagin (2005): modal counts, posterior-weighted sizes, average posterior of the
  # assigned class and the odds of correct classification.
  classification <- data.frame(class = 1:2, n_modal = tabulate(modal, 2L),
                               proportion_modal = tabulate(modal, 2L) / faithful_n,
                               estimated_n = colSums(z), estimated_proportion = shares,
                               average_posterior = diagonal,
                               odds_correct_classification = (diagonal / (1 - diagonal)) / (shares / (1 - shares)))
  expect_reference(get_results(diagnosed, "classification"), classification, reference = definitions,
                   covers = c("diagnostics()", "get_results.multilpa_diagnostics()"),
                   quantity = "classification summary", tolerance = 1e-10, measure = "relative", key = "class")
  # The average posterior of each class among the units assigned to each class;
  # the off-diagonal ones run down to 1e-5 and are probabilities: absolute.
  average <- expand.grid(assigned_class = 1:2, class = 1:2)
  average$n_assigned <- tabulate(modal, 2L)[average$assigned_class]
  average$average_posterior <- mapply(\(j, k) mean(z[modal == j, k]), average$assigned_class, average$class)
  expect_reference(get_results(diagnosed, "average_posteriors"), average, reference = definitions,
                   covers = c("diagnostics()", "get_results.multilpa_diagnostics()"),
                   quantity = "average posterior probabilities", tolerance = 1e-12,
                   key = c("assigned_class", "class"))
  # The numeric checks that have a definition: the smallest estimated class share,
  # the largest gap between modal and estimated shares, the relative entropy, the
  # lowest average posterior and the lowest odds of correct classification.
  odds <- (diagonal / (1 - diagonal)) / (shares / (1 - shares))
  checks <- data.frame(check = c("smallest class", "modal vs estimated proportions", "relative entropy",
                                 "lowest average posterior", "lowest odds of correct classification"),
                       value = c(min(shares), max(abs(tabulate(modal, 2L) / faithful_n - shares)),
                                 1 - entropy_sum / (faithful_n * log(2)), min(diagonal), min(odds)))
  expect_reference(subset(get_checks(diagnosed_fit), check %in% checks$check, select = c(check, value)),
                   checks, reference = definitions, covers = "get_checks()",
                   quantity = "checks with a definition", tolerance = 1e-10, measure = "relative", key = "check")
  expect_reference(subset(get_results(diagnosed, "checks"), check %in% checks$check, select = c(check, value)),
                   checks, reference = definitions,
                   covers = c("diagnostics()", "get_results.multilpa_diagnostics()"),
                   quantity = "checks with a definition, diagnostics()", tolerance = 1e-10,
                   measure = "relative", key = "check")
})

test_that("the replicated log likelihood is mclust's maximum and the tolerance its definition", {
  skip_without_reference("mclust")
  maximum <- faithful_maximum(diagnosed_fit)
  # Both EMs run to a relative change of 1e-12: the maxima agree to about 1e-10 relative.
  expect_reference(cbind(get_replication(diagnosed_fit), key = "fit"),
                   data.frame(key = "fit", log_likelihood = maximum$log_likelihood), reference = "mclust",
                   covers = c("lpa()", "get_replication()"), quantity = "best log likelihood of the starts",
                   tolerance = 1e-10, measure = "relative", key = "key")
  # The replication tolerance is 1e-6 per observed value (272 rows, two indicators).
  expect_reference(cbind(get_replication(diagnosed_fit), key = "fit"),
                   data.frame(key = "fit", tolerance = 1e-6 * faithful_n * 2),
                   reference = "independent: replication tolerance from its definition, 1e-6 per observed value",
                   covers = "get_replication()", quantity = "replication tolerance", tolerance = 1e-12,
                   measure = "relative", key = "key")
})

# --- starting_values() ------------------------------------------------------------

test_that("starting values read from a fit are mclust's maximum", {
  skip_without_reference("mclust")
  start <- starting_values(diagnosed_fit)
  # starting_values() numbers the indicators; mclust's table names them.
  expected <- transform(faithful_maximum(diagnosed_fit)$profiles,
                        indicator = match(indicator, colnames(faithful_x)))
  # Both EMs stop at a relative change of 1e-12 in the log likelihood; the
  # estimates agree to about the square root of that: 1e-6 relative.
  expect_reference(as.data.frame(start), expected, reference = "mclust",
                   covers = c("starting_values()", "as.data.frame.multilpa_start()"),
                   quantity = "starting means, variances and standard deviations", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
  expect_reference(get_results(start), expected, reference = "mclust",
                   covers = c("starting_values()", "get_results.multilpa_start()"),
                   quantity = "starting values, get_results()", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
  expect_reference(as.data.frame(summary(start)), expected, reference = "mclust",
                   covers = c("summary.multilpa_start()", "as.data.frame.multilpa_start()"),
                   quantity = "starting values, summary()", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
})

# --- sensitivity() ----------------------------------------------------------------

test_that("every seed's refit reaches mclust's maximum", {
  skip_without_reference("mclust")
  # The seeds are the refits' draws; the deterministic part checked is the log
  # likelihood each refit attains, which on this well separated two-profile
  # model is the one maximum mclust's EM reaches from latents' posteriors.
  refits <- sensitivity(diagnosed_fit, seeds = 1:4, tol = 1e-12)
  maximum <- faithful_maximum(diagnosed_fit)$log_likelihood
  # Every refit and mclust run to a relative change of 1e-12: about 1e-10 relative.
  expect_reference(as.data.frame(refits),
                   data.frame(seed = 1:4, log_likelihood = maximum), reference = "mclust",
                   covers = c("sensitivity()", "as.data.frame.latents_sensitivity()"),
                   quantity = "log likelihood reached by each seed", tolerance = 1e-10,
                   measure = "relative", key = "seed")
})

# --- align_labels() ---------------------------------------------------------------

test_that("aligned profiles are the reference's nearest profiles, with mclust's estimates", {
  skip_without_reference("mclust")
  feelings <- c("happy", "exhausted")
  first_days <- subset(student_esm, day <= 3)
  last_days <- subset(student_esm, day > 3)
  early <- lpa(first_days, feelings, n_profiles = 2, n_starts = 2, seed = 1, tol = 1e-12, max_iter = 100000L)
  late <- lpa(last_days, feelings, n_profiles = 2, n_starts = 2, seed = 1, tol = 1e-12, max_iter = 100000L)
  # mclust's EM continued from each fit's posteriors: estimates in that fit's own numbering.
  maximum <- function(fit, data) {
    posteriors <- get_results(fit, "posteriors")
    z <- unclass(tapply(posteriors$posterior, list(posteriors$row, posteriors$profile), sum))
    continued <- mclust::meVVI(data = as.matrix(data[feelings]), z = unname(z),
                               control = mclust::emControl(tol = c(1e-12, 1e-12)))
    list(means = continued$parameters$mean, variances = apply(continued$parameters$variance$sigma, 3L, diag))
  }
  reference <- maximum(early, first_days)
  other <- maximum(late, last_days)
  # The permutation of the late profiles nearest the early ones: squared mean
  # differences in the early profiles' standard deviations, every permutation tried.
  permutations <- list(c(1L, 2L), c(2L, 1L))
  distance <- vapply(permutations, \(permutation) {
    sum(((other$means[, permutation] - reference$means) / sqrt(reference$variances))^2)
  }, numeric(1))
  nearest <- permutations[[which.min(distance)]]
  expected <- data.frame(profile = rep(1:2, each = 2), indicator = rep(feelings, 2),
                         mean = as.vector(other$means[, nearest]), variance = as.vector(other$variances[, nearest]))
  # Two EM runs stopped at a relative change of 1e-12: the estimates agree to about 1e-6 relative.
  expect_reference(get_results(align_labels(late, early), "profiles"), expected,
                   reference = "mclust", covers = c("align_labels()", "lpa()"),
                   quantity = "late profiles numbered like the early ones", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"),
                   note = "Reference: mclust::meVVI() continued from each fit's posteriors; the numbering by the nearest permutation of the means, computed in the test.")
})

# --- fit_staged() -----------------------------------------------------------------

staged_vars <- c("happy", "exhausted")
staged <- fit_staged(student_esm, staged_vars, "student", n_profiles = 2, n_group_classes = 2,
                     n_starts = 3, seed = 1, tol = 1e-12, max_iter = 100000L)

test_that("the first stage of fit_staged() is mclust's single-level maximum", {
  skip_without_reference("mclust")
  x <- as.matrix(student_esm[staged_vars])
  profiles <- get_results(staged, "profiles")
  means <- matrix(profiles$mean, 2L, 2L)
  variances <- matrix(profiles$variance, 2L, 2L)
  # A start for mclust: Bayes' rule with equal shares at the held measurement.
  log_density <- vapply(1:2, \(k) stats::dnorm(x[, 1L], means[1L, k], sqrt(variances[1L, k]), log = TRUE) +
                          stats::dnorm(x[, 2L], means[2L, k], sqrt(variances[2L, k]), log = TRUE), numeric(nrow(x)))
  z <- exp(log_density - apply(log_density, 1L, max))
  # From this start mclust's EM creeps (it stops short at a relative change of
  # 1e-15): it runs 20000 iterations without a stopping rule instead.
  continued <- mclust::meVVI(data = x, z = z / rowSums(z),
                             control = mclust::emControl(tol = c(0, 0), itmax = c(20000, 20000)))
  expected <- data.frame(profile = rep(1:2, each = 2), indicator = rep(staged_vars, 2),
                         mean = as.vector(continued$parameters$mean),
                         variance = as.vector(apply(continued$parameters$variance$sigma, 3L, diag)))
  # fit_staged() stops at a relative change of 1e-12 in the log likelihood:
  # its estimates then agree with the maximum to about the square root, 1e-6 relative.
  # Recorded as a known difference after the first ci run (2026-11-20): one variance differed by
  # 1.3e-6 relative while the two log likelihoods differ by 2.2e-10 (2.3e-14 relative), below the
  # rounding of a sum over 2,582 rows. On ci, mclust's EM settles at a slightly different point of a
  # ridge that is flat to rounding (on macOS it settles within 3e-9 of latents); the tolerance is unchanged.
  expect_known_difference(profiles, expected, reference = "mclust",
                          covers = c("fit_staged()", "get_results.multilpa()"),
                          quantity = "held measurement: means and variances", tolerance = 1e-6,
                          measure = "relative", key = c("profile", "indicator"),
                          reason = "The first-stage likelihood is flat to rounding along one direction: latents' and mclust's points differ by 2.2e-10 in log likelihood (checked below to 1e-12 relative), which leaves the profile-2 variance of exhausted determined only to about 1e-6. Both are the maximum within double precision.")
  # The likelihood is flat along this ridge: the two maxima agree to rounding.
  expect_reference(subset(get_results(staged, "stages"), stage == "measurement", select = c(stage, log_likelihood)),
                   data.frame(stage = "measurement", log_likelihood = continued$loglik), reference = "mclust",
                   covers = "fit_staged()", quantity = "first-stage log likelihood", tolerance = 1e-12,
                   measure = "relative", key = "stage")
})

test_that("the second stage of fit_staged() maximises the two-level likelihood with the measurement held", {
  x <- as.matrix(student_esm[staged_vars])
  group <- student_esm$student
  profiles <- get_results(staged, "profiles")
  means <- matrix(profiles$mean, 2L, 2L)
  variances <- matrix(profiles$variance, 2L, 2L)
  log_density <- vapply(1:2, \(k) stats::dnorm(x[, 1L], means[1L, k], sqrt(variances[1L, k]), log = TRUE) +
                          stats::dnorm(x[, 2L], means[2L, k], sqrt(variances[2L, k]), log = TRUE), numeric(nrow(x)))
  log_sum_exp <- function(a) {
    top <- apply(a, 1L, max)
    top + log(rowSums(exp(a - top)))
  }
  # Students are the groups: sum over the two group classes of the product over
  # a student's beeps of the profile mixture, with the measurement held.
  log_likelihood <- function(theta) {
    profile_first <- stats::plogis(theta[1:2])
    class_first <- stats::plogis(theta[[3L]])
    per_class <- vapply(1:2, \(h) {
      rowsum(log_sum_exp(sweep(log_density, 2L, log(c(profile_first[h], 1 - profile_first[h])), "+")),
             group, reorder = FALSE)[, 1L]
    }, numeric(length(unique(group))))
    sum(log_sum_exp(sweep(per_class, 2L, log(c(class_first, 1 - class_first)), "+")))
  }
  # A neutral start that names the first group class the one richer in profile 1.
  optimum <- stats::optim(stats::qlogis(c(0.7, 0.3, 0.5)), \(theta) -log_likelihood(theta),
                          method = "BFGS", control = list(reltol = 1e-15, maxit = 10000L))
  expect_equal(optimum$convergence, 0L)
  maximiser <- "independent: two-level likelihood with the measurement held, maximised by stats::optim()"
  # Both at the maximum, where the likelihood is flat to second order: 1e-10 relative.
  expect_reference(as.numeric(logLik(staged)), -optimum$value, reference = maximiser,
                   covers = c("fit_staged()", "logLik.multilpa()"),
                   quantity = "second-stage log likelihood", tolerance = 1e-10, measure = "relative")
  shares <- stats::plogis(optimum$par)
  expected <- data.frame(group_class = rep(1:2, each = 2), profile = rep(1:2, 2),
                         probability = c(shares[1L], 1 - shares[1L], shares[2L], 1 - shares[2L]),
                         group_class_probability = rep(c(shares[3L], 1 - shares[3L]), each = 2))
  # BFGS stops on a relative change of 1e-15 in the objective with the logits'
  # gradient near 1e-5, which leaves the probabilities about 1e-7 from the
  # maximum: 1e-6 absolute.
  expect_reference(get_results(staged, "profile_probabilities"), expected, reference = maximiser,
                   covers = c("fit_staged()", "get_results.multilpa()"),
                   quantity = "second-stage profile and group-class probabilities", tolerance = 1e-6,
                   measure = "absolute", key = c("group_class", "profile"))
})
