# Functions for compare.R: fit single-level latent class models with five
# independent R implementations (glca, randomLCA, BayesLCA's EM, multilevLCA,
# depmixS4) and with lca(), put every solution in one shape, and compare the
# maximized likelihood, class shares, item response probabilities and, where a
# package returns them per row, posteriors.
#
# Tolerances, fixed before any comparison was run: the log likelihood to 1e-6;
# shares, probabilities and posteriors to 1e-4 (each package's EM stops on its
# own relative criterion, set here to 1e-10 or tighter). Every package runs on
# one core.

.lca_case <- function(name, data, items, classes, references) {
  list(name = name, data = data, items = items, classes = classes, references = references)
}

.expand <- function(patterns, items) {
  patterns[rep(seq_len(nrow(patterns)), patterns$freq), items, drop = FALSE]
}

#' The datasets and class counts compared
#'
#' Binary data are coded 0/1, gss82 keeps its codes 1..C. Frequency tables
#' are expanded to one row per respondent.
lca_cases <- function() {
  binary <- c("glca", "randomLCA", "BayesLCA", "multilevLCA", "depmixS4")
  polytomous <- c("glca", "depmixS4")
  data_of <- function(name, package) {
    get(utils::data(list = name, package = package, envir = environment()), envir = environment())
  }
  zero_one <- function(frame) as.data.frame(lapply(frame, function(v) { v <- as.integer(v); v - min(v) }))
  carcinoma <- zero_one(data_of("carcinoma", "poLCA"))
  values <- zero_one(data_of("values", "poLCA"))
  cheating <- zero_one(data_of("cheating", "poLCA")[c("LIEEXAM", "LIEPAPER", "FRAUD", "COPYEXAM")])
  alzheimer <- data_of("Alzheimer", "BayesLCA")
  myocardial <- data_of("myocardial", "randomLCA")
  dentistry <- data_of("dentistry", "randomLCA")
  hivtests <- data_of("hivtests", "randomLCA")
  gss82 <- as.data.frame(lapply(data_of("gss82", "poLCA"), as.integer))
  help <- help_lca_data()
  myocardial_items <- setdiff(names(myocardial), "freq")
  dentistry_items <- setdiff(names(dentistry), "freq")
  hiv_items <- setdiff(names(hivtests), "freq")
  list(
    .lca_case("carcinoma_k2", carcinoma, names(carcinoma), 2L, binary),
    .lca_case("carcinoma_k3", carcinoma, names(carcinoma), 3L, binary),
    .lca_case("values_k2", values, names(values), 2L, binary),
    .lca_case("cheating_k2", cheating, names(cheating), 2L, binary),
    .lca_case("alzheimer_k2", alzheimer, names(alzheimer), 2L, binary),
    .lca_case("alzheimer_k3", alzheimer, names(alzheimer), 3L, binary),
    .lca_case("myocardial_k2", .expand(myocardial, myocardial_items), myocardial_items, 2L, binary),
    .lca_case("dentistry_k2", .expand(dentistry, dentistry_items), dentistry_items, 2L, binary),
    .lca_case("dentistry_k3", .expand(dentistry, dentistry_items), dentistry_items, 3L, binary),
    .lca_case("hivtests_k2", .expand(hivtests, hiv_items), hiv_items, 2L, binary),
    .lca_case("gss82_k2", gss82, names(gss82), 2L, polytomous),
    .lca_case("gss82_k3", gss82, names(gss82), 3L, polytomous),
    # randomLCA 1.1.4 refuses this model ("not identifiable": one residual
    # degree of freedom, a response probability at 1); the blog's Example
    # 8.23 fitted it with the randomLCA of 2011.
    .lca_case("help_k3", help, names(help), 3L, setdiff(binary, "randomLCA")))
}

#' The HELP study indicators of the SAS and R blog's Example 8.21
#' (PROVENANCE.txt): homeless, cesd > 20, satreat, linkstatus; 431 rows
help_lca_data <- function() {
  utils::read.csv(file.path(lca_packages_root(), "help-lca.csv"))
}

#' The kit folder, from the project root (compare.R) or from equivalence/
#' (testthat)
lca_packages_root <- function() {
  here <- "lca-packages"
  if (dir.exists(here)) here else file.path("equivalence", here)
}

.printed_row <- function(case, source, quantity, printed, obtained) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  tolerance <- 0.5 * 10^(-digits) + 1e-9
  value <- as.numeric(printed)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, agrees = abs(obtained - value) <= tolerance)
}

#' The SAS and R blog's Example 8.21 (Kleinman and Horton 2011): poLCA's
#' printed three-class fit, and SAS PROC LCA's log likelihood and item
#' response probabilities. SAS's class prevalences (.2163 / .0785 / .7052)
#' differ from poLCA's in the fourth decimal and are not compared.
compare_published_help <- function() {
  data <- help_lca_data()
  fit <- lca(data, names(data), n_classes = 3L, n_starts = 50L, seed = 1L, tol = 1e-12,
             max_iter = 100000L)
  theta <- stats::coef(fit)
  shares <- vapply(1:3, function(j) theta[[sprintf("profile.probability.profile_%d.group_class_1", j)]],
                   numeric(1))
  housed <- vapply(1:3, function(j) theta[[sprintf("measurement.response.profile_%d.homeless:0", j)]],
                   numeric(1))
  order <- order(shares)
  criteria <- get_results(fit, "information_criteria", format = "wide")
  # One published source: the post prints both programs' output.
  source <- "Kleinman and Horton (2011), SAS and R blog, Example 8.21"
  rbind(
    .printed_row("help_k3", source, "poLCA: log likelihood", "-1032.484",
                 as.numeric(stats::logLik(fit))),
    data.frame(case = "help_k3", reference = source, quantity = "poLCA: number of parameters",
               reference_value = 14, latents = fit$n_parameters,
               difference = abs(fit$n_parameters - 14), tolerance = 0, n_values = 1L,
               agrees = fit$n_parameters == 14),
    .printed_row("help_k3", source, "poLCA: AIC", "2092.967", criteria$aic),
    .printed_row("help_k3", source, "poLCA: BIC", "2149.893", criteria$bic_individual),
    do.call(rbind, Map(function(printed, j) {
      .printed_row("help_k3", source, sprintf("poLCA: class share %d (sorted)", j), printed,
                   shares[order][[j]])
    }, c("0.0785", "0.2162", "0.7053"), 1:3)),
    .printed_row("help_k3", source, "SAS PROC LCA: log likelihood", "-1032.48",
                 as.numeric(stats::logLik(fit))),
    do.call(rbind, Map(function(printed, j) {
      .printed_row("help_k3", source,
                   sprintf("SAS PROC LCA: P(housed) in the class with share %d (sorted)", j),
                   printed, housed[order][[j]])
    }, c("1.0000", "0.2703", "0.5625"), 1:3)))
}

.levels <- function(case) {
  stats::setNames(lapply(case$items, function(v) sort(unique(case$data[[v]]))), case$items)
}

#' A solution in one shape
#' @return `log_likelihood`, `shares` (K), `probabilities` (per item, K x C,
#'   columns in the items' sorted category order) and `posteriors` (n x K or
#'   `NULL`).
.solution <- function(log_likelihood, shares, probabilities, posteriors = NULL) {
  list(log_likelihood = log_likelihood, shares = as.numeric(shares),
       probabilities = probabilities, posteriors = posteriors)
}

fit_latents_lca <- function(case) {
  fit <- lca(case$data, case$items, n_classes = case$classes, n_starts = 50L, seed = 1L,
             tol = 1e-12, max_iter = 100000L)
  theta <- stats::coef(fit)
  k <- case$classes
  levels <- .levels(case)
  probabilities <- lapply(case$items, function(v) {
    t(vapply(seq_len(k), function(j) vapply(levels[[v]], function(c)
      theta[[sprintf("measurement.response.profile_%d.%s:%s", j, v, c)]], numeric(1)),
      numeric(length(levels[[v]]))))
  })
  names(probabilities) <- case$items
  .solution(as.numeric(stats::logLik(fit)),
            vapply(seq_len(k), function(j)
              theta[[sprintf("profile.probability.profile_%d.group_class_1", j)]], numeric(1)),
            probabilities, fit$subject_posteriors)
}

.binary_probabilities <- function(p_one, items) {
  stats::setNames(lapply(seq_along(items), function(i) cbind(1 - p_one[, i], p_one[, i])), items)
}

fit_glca <- function(case) {
  shifted <- as.data.frame(lapply(case$data[case$items], function(v) v - min(v) + 1L))
  # glca reads the items through its own item() in the formula.
  formula <- stats::as.formula(sprintf("item(%s) ~ 1", paste(case$items, collapse = ", ")),
                               env = asNamespace("glca"))
  fit <- glca::glca(formula, data = shifted, nclass = case$classes, n.init = 50L,
                    maxiter = 100000L, eps = 1e-12, seed = 1L, verbose = FALSE, std.err = FALSE)
  rho <- fit$param$rho$ALL
  .solution(fit$gof$loglik, fit$param$gamma,
            stats::setNames(lapply(case$items, function(v) unname(rho[[v]])), case$items),
            unname(as.matrix(fit$posterior$ALL)))
}

fit_randomlca <- function(case) {
  fit <- randomLCA::randomLCA(case$data[case$items], nclass = case$classes, penalty = 0,
                              EMtol = 1e-12, notrials = 50L, seed = 1L, cores = 1L,
                              calcSE = FALSE)
  .solution(fit$logLik, fit$classp, .binary_probabilities(fit$outcomep, case$items))
}

fit_bayeslca <- function(case) {
  # Uniform Dirichlet and Beta priors (alpha = beta = delta = 1): the posterior
  # mode is the maximum likelihood estimate. logpost adds the log density of
  # those priors, which is constant: 0 for each Beta(1, 1) and lgamma(G) for
  # the Dirichlet(1, ..., 1) on G class shares (log 2 = 0.693 at G = 3).
  fit <- withr::with_seed(1L, BayesLCA::blca.em(as.matrix(case$data[case$items]), case$classes,
                                                alpha = 1, beta = 1, delta = 1, restarts = 50L,
                                                conv = 1e-12, iter = 100000L, verbose = FALSE))
  .solution(fit$logpost - lgamma(case$classes), fit$classprob,
            .binary_probabilities(fit$itemprob, case$items))
}

fit_multilevlca <- function(case) {
  cores <- parallel::detectCores()
  fit <- multilevLCA::multiLCA(case$data[case$items], case$items, case$classes,
                               numFreeCores = cores - 1L, sequential = TRUE, tol = 1e-12,
                               maxIter = 100000L, verbose = FALSE)
  .solution(utils::tail(as.numeric(fit$LLKSeries), 1L), fit$vPi,
            .binary_probabilities(t(fit$mPhi), case$items))
}

fit_depmixs4 <- function(case) {
  frame <- as.data.frame(lapply(case$data[case$items], factor))
  responses <- lapply(case$items, function(v) stats::as.formula(paste(v, "~ 1")))
  model <- depmixS4::mix(responses, data = frame, nstates = case$classes,
                         family = rep(list(depmixS4::multinomial("identity")), length(case$items)))
  # depmixS4's EM stops on a relative criterion; at 1e-12 it stopped 3e-4 short
  # in the posteriors of a flat model (dentistry, three classes), so it runs
  # to 1e-15. There it ends when the likelihood "decreases" by floating-point
  # noise (1e-11), which depmixS4 reports as a warning: those are counted and
  # reported, any other warning propagates. Each fit's console line is not kept.
  floor_notices <- 0L
  fits <- lapply(seq_len(20L), function(s) {
    withr::with_seed(s, tryCatch(withCallingHandlers({
      printed <- utils::capture.output(fitted <- depmixS4::fit(
        model, emcontrol = depmixS4::em.control(maxit = 1000000L, tol = 1e-15, random.start = TRUE),
        verbose = FALSE))
      fitted
    }, warning = function(w) {
      if (grepl("^Log likelihood decreased on iteration", conditionMessage(w))) {
        floor_notices <<- floor_notices + 1L
        invokeRestart("muffleWarning")
      }
    }), error = function(e) {
      message(sprintf("depmixS4 start %d failed: %s", s, conditionMessage(e)))
      NULL
    }))
  })
  if (floor_notices) message(sprintf("%s: depmixS4 reached the floating-point floor in %d fits",
                                     case$name, floor_notices))
  fits <- Filter(Negate(is.null), fits)
  best <- fits[[which.max(vapply(fits, function(f) as.numeric(depmixS4::logLik(f)), numeric(1)))]]
  posteriors <- as.matrix(depmixS4::posterior(best, type = "smoothing"))
  # At an EM fixed point the shares are the mean posteriors and each item's
  # probabilities the posterior-weighted category frequencies.
  levels <- .levels(case)
  probabilities <- lapply(case$items, function(v) {
    indicator <- outer(case$data[[v]], levels[[v]], `==`) * 1
    crossprod(posteriors, indicator) / colSums(posteriors)
  })
  names(probabilities) <- case$items
  .solution(as.numeric(depmixS4::logLik(best)), colMeans(posteriors), probabilities, posteriors)
}

.fitters <- list(glca = fit_glca, randomLCA = fit_randomlca, BayesLCA = fit_bayeslca,
                 multilevLCA = fit_multilevlca, depmixS4 = fit_depmixs4)

#' All orderings of `k` classes
#' @noRd
.permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(.permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

#' The class order of `reference` that brings it closest to `ours`
#' @noRd
.match_classes <- function(ours, reference) {
  perms <- .permutations(length(ours$shares))
  gap <- vapply(perms, function(perm) {
    max(abs(reference$shares[perm] - ours$shares),
        vapply(names(ours$probabilities), function(v)
          max(abs(reference$probabilities[[v]][perm, , drop = FALSE] - ours$probabilities[[v]])),
          numeric(1)))
  }, numeric(1))
  perms[[which.min(gap)]]
}

#' Compare one package's solution with latents'
#' @noRd
compare_solution <- function(case, package, ours, reference) {
  perm <- .match_classes(ours, reference)
  probability_gap <- max(vapply(names(ours$probabilities), function(v)
    max(abs(reference$probabilities[[v]][perm, , drop = FALSE] - ours$probabilities[[v]])),
    numeric(1)))
  k <- length(ours$shares)
  n_probabilities <- sum(lengths(ours$probabilities))
  rows <- data.frame(
    case = case$name, reference = package,
    quantity = c("maximized log likelihood", sprintf("class shares (%d values)", k),
                 sprintf("item response probabilities (%d values)", n_probabilities)),
    reference_value = c(reference$log_likelihood, NA, NA),
    latents = c(ours$log_likelihood, NA, NA),
    difference = c(abs(reference$log_likelihood - ours$log_likelihood),
                   max(abs(reference$shares[perm] - ours$shares)), probability_gap),
    tolerance = c(1e-6, 1e-4, 1e-4), n_values = c(1L, k, n_probabilities))
  if (!is.null(reference$posteriors)) {
    stopifnot("one posterior row per respondent" = nrow(reference$posteriors) == nrow(ours$posteriors))
    rows <- rbind(rows, data.frame(
      case = case$name, reference = package,
      quantity = sprintf("posteriors (%d rows x %d classes)", nrow(ours$posteriors), k),
      reference_value = NA, latents = NA,
      difference = max(abs(reference$posteriors[, perm, drop = FALSE] - ours$posteriors)),
      tolerance = 1e-4, n_values = length(ours$posteriors)))
  }
  rows$agrees <- rows$difference <= rows$tolerance
  rows
}

#' Fit every case with latents and its references and compare
#' @return One row per compared quantity, with `fn = "lca"`.
compare_all_lca <- function(cases = lca_cases()) {
  out <- do.call(rbind, lapply(cases, function(case) {
    ours <- fit_latents_lca(case)
    do.call(rbind, lapply(case$references, function(package) {
      compare_solution(case, package, ours, .fitters[[package]](case))
    }))
  }))
  out$kind <- "R package"
  published <- compare_published_help()
  published$kind <- "published"
  out <- rbind(out, published)
  out$fn <- "lca"
  out
}
