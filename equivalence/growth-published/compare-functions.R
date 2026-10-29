# Functions for compare.R: the growth mixtures of Proust-Lima, Philipps and
# Liquet (2017, JSS 78(2), section 6) on paquid, fitted by mixture_regression()
# and compared with the published summary table and posterior classification,
# and with the current lcmm (hlme) fitting the same models.
#
# The model: normMMSE on a quadratic in age65, class-specific; CEP shared by
# the classes; correlated random intercept, slope and quadratic with one
# covariance for every class; one residual variance. In latents:
# common = ~CEP, random = ~ age65 + age65sq, random_covariance = "equal",
# variance = "equal". (age65sq is age65^2, the paper's poly(raw = TRUE).)
#
# Tolerances, fixed before the comparison: published values to half a unit of
# their last printed digit, counts exactly; lcmm's likelihood to 1e-4 (its
# convergence criteria are on the parameters, the score and 1e-4 of the
# likelihood) and its estimates to 1e-3 of max(1, |value|).

paquid_data <- function(directory) {
  data <- utils::read.csv(file.path(directory, "paquid-normmmse.csv"))
  data$age65sq <- data$age65^2
  data[!is.na(data$normMMSE), ]
}

fit_paquid <- function(data, classes) {
  mixture_regression(normMMSE ~ age65 + age65sq + CEP, data, classes, id = "ID",
                     class_level = "group", common = ~CEP, random = ~ age65 + age65sq,
                     random_covariance = "equal", variance = "equal",
                     n_starts = if (classes == 1L) 1L else 50L, seed = 1L, tol = 1e-10,
                     max_iter = 100000L)
}

.printed <- function(case, source, quantity, printed, obtained, kind = "published") {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  tolerance <- 0.5 * 10^(-digits) + 1e-9
  value <- as.numeric(printed)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, agrees = abs(obtained - value) <= tolerance, kind = kind)
}

#' Modal-class proportions (percent) and the mean posterior of each modal class
#' @noRd
.classification <- function(fit) {
  tau <- fit$expectation$tau
  modal <- max.col(tau, ties.method = "first")
  k <- ncol(tau)
  counts <- tabulate(modal, k)
  mean_posterior <- t(vapply(seq_len(k), function(j) colMeans(tau[modal == j, , drop = FALSE]),
                             numeric(k)))
  list(counts = counts, percent = 100 * counts / length(modal), mean_posterior = mean_posterior)
}

#' The published summary table and posterior classification
#' @noRd
compare_published_growth <- function(fits) {
  source <- "Proust-Lima, Philipps and Liquet (2017), JSS 78(2), section 6"
  published <- list(
    k1 = list(loglik = "-8920.623", npm = 11L, bic = "17909.61", percent = "100.0"),
    k2 = list(loglik = "-8899.228", npm = 15L, bic = "17891.67", percent = c("12.4", "87.6")),
    k3 = list(loglik = "-8890.148", npm = 19L, bic = "17898.37", percent = c("0.8", "11.8", "87.4")))
  rows <- do.call(rbind, lapply(names(published), function(name) {
    fit <- fits[[name]]
    p <- published[[name]]
    case <- paste0("paquid_", name)
    percent <- sort(.classification(fit)$percent)
    rbind(
      .printed(case, source, "log likelihood", p$loglik, as.numeric(stats::logLik(fit))),
      data.frame(case = case, reference = source, quantity = "number of parameters",
                 reference_value = p$npm, latents = fit$n_parameters,
                 difference = abs(fit$n_parameters - p$npm), tolerance = 0, n_values = 1L,
                 agrees = fit$n_parameters == p$npm, kind = "published"),
      .printed(case, source, "BIC", p$bic, stats::BIC(fit)),
      do.call(rbind, lapply(seq_along(p$percent), function(j)
        .printed(case, source, sprintf("modal class percentage %d (sorted)", j),
                 sort(p$percent)[[j]], percent[[j]]))))
  }))
  # postprob(m2): 62 and 438 subjects; mean posteriors 0.8054 / 0.1946 and
  # 0.1270 / 0.8730 for the small and the large class.
  two <- .classification(fits$k2)
  small <- which.min(two$counts)
  large <- 3L - small
  rbind(rows,
        data.frame(case = "paquid_k2", reference = source,
                   quantity = c("subjects in the small class", "subjects in the large class"),
                   reference_value = c(62, 438), latents = two$counts[c(small, large)],
                   difference = abs(two$counts[c(small, large)] - c(62, 438)), tolerance = 0,
                   n_values = 1L, agrees = two$counts[c(small, large)] == c(62, 438),
                   kind = "published"),
        .printed("paquid_k2", source, "mean posterior of the small class, small class", "0.8054",
                 two$mean_posterior[small, small]),
        .printed("paquid_k2", source, "mean posterior of the large class, large class", "0.8730",
                 two$mean_posterior[large, large]))
}

#' The same models by the current lcmm (hlme), as the paper fits them
#' @noRd
fit_lcmm <- function(data) {
  fixed <- normMMSE ~ age65 + age65sq + CEP
  m1 <- lcmm::hlme(fixed, random = ~ age65 + age65sq, subject = "ID", ng = 1, data = data,
                   verbose = FALSE)
  m2 <- lcmm::hlme(fixed, random = ~ age65 + age65sq, mixture = ~ age65 + age65sq,
                   subject = "ID", ng = 2, data = data, B = m1, verbose = FALSE)
  # gridsearch() re-evaluates its model call by the function's name, so hlme
  # must be visible here under that name.
  hlme <- lcmm::hlme
  m3 <- withr::with_seed(1L, lcmm::gridsearch(
    hlme(fixed, random = ~ age65 + age65sq, mixture = ~ age65 + age65sq, subject = "ID",
         ng = 3, data = data, verbose = FALSE), rep = 30, maxiter = 15, minit = m1))
  list(k1 = m1, k2 = m2, k3 = m3)
}

#' Fixed effects of a latents fit, sorted per class by intercept, with CEP
#' @noRd
.latents_fixed <- function(fit) {
  theta <- stats::coef(fit)
  k <- ncol(fit$expectation$tau)
  per_class <- t(vapply(seq_len(k), function(j) c(
    theta[[sprintf("coefficient.class_%d.(Intercept)", j)]],
    theta[[sprintf("coefficient.class_%d.age65", j)]],
    theta[[sprintf("coefficient.class_%d.age65sq", j)]]), numeric(3)))
  list(per_class = per_class[order(per_class[, 1L]), , drop = FALSE],
       cep = theta[["coefficient.common.CEP"]])
}

#' lcmm's fixed effects in the same shape (its class 1 intercept is absolute)
#' @noRd
.lcmm_fixed <- function(model) {
  k <- model$ng
  # hlme lists the class-membership logits first, also named "intercept
  # class1", ...; the fixed effects follow.
  best <- if (k > 1L) model$best[-seq_len(k - 1L)] else model$best
  pick <- function(pattern) best[grep(pattern, names(best))]
  intercept <- if (k == 1L) best[["intercept"]] else pick("^intercept class")
  age <- if (k == 1L) best[["age65"]] else pick("^age65 class")
  age2 <- if (k == 1L) best[["age65sq"]] else pick("^age65sq class")
  per_class <- cbind(intercept, age, age2)
  list(per_class = unname(per_class[order(per_class[, 1L]), , drop = FALSE]), cep = best[["CEP"]])
}

compare_lcmm <- function(fits, lcmm_fits) {
  do.call(rbind, lapply(names(fits), function(name) {
    ours <- fits[[name]]
    reference <- lcmm_fits[[name]]
    a <- .latents_fixed(ours)
    b <- .lcmm_fixed(reference)
    fixed_gap <- max(abs(c(a$per_class, a$cep) - c(b$per_class, b$cep)) /
                       pmax(1, abs(c(b$per_class, b$cep))))
    rows <- data.frame(
      case = paste0("paquid_", name), reference = "lcmm",
      quantity = c("maximized log likelihood", "number of parameters",
                   sprintf("fixed effects (%d values)", length(a$per_class) + 1L)),
      reference_value = c(reference$loglik, length(reference$best), NA),
      latents = c(as.numeric(stats::logLik(ours)), ours$n_parameters, NA),
      difference = c(abs(reference$loglik - as.numeric(stats::logLik(ours))),
                     abs(length(reference$best) - ours$n_parameters), fixed_gap),
      tolerance = c(1e-4, 0, 1e-3), n_values = c(1L, 1L, length(a$per_class) + 1L),
      kind = "R package")
    rows$agrees <- rows$difference <= rows$tolerance
    rows
  }))
}

# -- lcmm's data_hlme ------------------------------------------------------------
# lcmm's simulated example data (100 subjects, 326 rows) with the model of its
# hlme() help page: Y ~ Time * X1, intercept and Time class-specific, X1 and
# Time:X1 shared, a correlated random intercept and slope with one covariance
# for every class, one residual variance, class membership on X2 + X3.

fit_hlme_data <- function(data, classes) {
  mixture_regression(Y ~ Time + X1 + Time:X1, data, classes, id = "ID", class_level = "group",
                     common = ~ X1 + Time:X1, membership = if (classes == 1L) ~1 else ~ X2 + X3,
                     random = ~Time, random_covariance = "equal", variance = "equal",
                     n_starts = if (classes == 1L) 1L else 50L, seed = 1L, tol = 1e-10,
                     max_iter = 100000L)
}

#' lcmm's hlme() fits, by its grid search from the one-class fit
#' @noRd
fit_lcmm_hlme_data <- function(data) {
  m1 <- lcmm::hlme(Y ~ Time * X1, random = ~Time, subject = "ID", ng = 1, data = data,
                   verbose = FALSE)
  # gridsearch() re-evaluates its model call by the function's name.
  hlme <- lcmm::hlme
  search <- function(classes) withr::with_seed(1L, lcmm::gridsearch(
    hlme(Y ~ Time * X1, mixture = ~Time, random = ~Time, classmb = ~ X2 + X3, subject = "ID",
         ng = classes, data = data, verbose = FALSE), rep = 30, maxiter = 15, minit = m1))
  list(k1 = m1, k2 = search(2L), k3 = search(3L))
}

#' Every estimate of a data_hlme model in one shape: classes ordered by their
#' intercept; class-specific intercept and Time; shared X1 and Time:X1; the
#' membership logits (Intercept, X2, X3) of each class against the first
#' ordered class; the random covariance (var I, cov, var Time); residual SD
#' @noRd
.hlme_shape <- function(intercept, time, x1, time_x1, membership, covariance, sigma) {
  order <- order(intercept)
  membership <- membership[order, , drop = FALSE]
  relative <- sweep(membership, 2L, membership[1L, ], "-")[-1L, , drop = FALSE]
  c(intercept[order], time[order], x1, time_x1, as.vector(t(relative)), covariance, sigma)
}

.hlme_latents <- function(fit) {
  theta <- stats::coef(fit)
  k <- ncol(fit$expectation$tau)
  value <- function(name) theta[[name]]
  per_class <- function(term) vapply(seq_len(k), function(j)
    value(sprintf("coefficient.class_%d.%s", j, term)), numeric(1))
  membership <- t(vapply(seq_len(k), function(j) if (j == 1L) c(0, 0, 0) else vapply(
    c("(Intercept)", "X2", "X3"), function(term) value(sprintf("membership.class_%d.%s", j, term)),
    numeric(1)), numeric(3)))
  l11 <- exp(value("random.log_cholesky.shared.(Intercept):(Intercept)"))
  l21 <- value("random.cholesky.shared.Time:(Intercept)")
  l22 <- exp(value("random.log_cholesky.shared.Time:Time"))
  .hlme_shape(per_class("(Intercept)"), per_class("Time"), value("coefficient.common.X1"),
              value("coefficient.common.Time:X1"), if (k > 1L) membership else matrix(0, 1L, 3L),
              c(l11^2, l11 * l21, l21^2 + l22^2), exp(value("log_sigma.shared")))
}

.hlme_lcmm <- function(model) {
  best <- model$best
  k <- model$ng
  pick <- function(name) best[[name]]
  if (k == 1L) {
    intercept <- pick("intercept")
    time <- pick("Time")
    membership <- matrix(0, 1L, 3L)
  } else {
    # "intercept classj" names the membership logit (j < k) and, after it, the
    # class's fixed intercept.
    named <- function(j) best[names(best) == sprintf("intercept class%d", j)]
    intercept <- vapply(seq_len(k), function(j) rev(named(j))[[1L]], numeric(1))
    time <- vapply(seq_len(k), function(j) pick(sprintf("Time class%d", j)), numeric(1))
    # Membership logits of classes 1..k-1 against class k, listed first.
    membership <- rbind(t(vapply(seq_len(k - 1L), function(j) c(
      named(j)[[1L]],
      pick(sprintf("X2 class%d", j)), pick(sprintf("X3 class%d", j))), numeric(3))), c(0, 0, 0))
  }
  .hlme_shape(intercept, time, pick("X1"), pick("Time:X1"), membership,
              unname(best[grep("^varcov", names(best))]), pick("stderr"))
}

compare_lcmm_hlme_data <- function(fits, lcmm_fits) {
  do.call(rbind, lapply(names(fits), function(name) {
    ours <- fits[[name]]
    reference <- lcmm_fits[[name]]
    a <- .hlme_latents(ours)
    b <- .hlme_lcmm(reference)
    gap <- max(abs(a - b) / pmax(1, abs(b)))
    rows <- data.frame(
      case = paste0("data_hlme_", name), reference = "lcmm",
      quantity = c("maximized log likelihood", "number of parameters",
                   sprintf("estimates (%d values)", length(a))),
      reference_value = c(reference$loglik, length(reference$best), NA),
      latents = c(as.numeric(stats::logLik(ours)), ours$n_parameters, NA),
      difference = c(abs(reference$loglik - as.numeric(stats::logLik(ours))),
                     abs(length(reference$best) - ours$n_parameters), gap),
      tolerance = c(1e-4, 0, 1e-3), n_values = c(1L, 1L, length(a)), kind = "R package")
    rows$agrees <- rows$difference <= rows$tolerance
    rows
  }))
}

compare_all_growth_published <- function(directory) {
  data <- paquid_data(directory)
  fits <- list(k1 = fit_paquid(data, 1L), k2 = fit_paquid(data, 2L), k3 = fit_paquid(data, 3L))
  hlme_data <- lcmm::data_hlme
  hlme_fits <- list(k1 = fit_hlme_data(hlme_data, 1L), k2 = fit_hlme_data(hlme_data, 2L),
                    k3 = fit_hlme_data(hlme_data, 3L))
  out <- rbind(compare_published_growth(fits), compare_lcmm(fits, fit_lcmm(data)),
               compare_lcmm_hlme_data(hlme_fits, fit_lcmm_hlme_data(hlme_data)))
  out$fn <- "growth"
  out
}
