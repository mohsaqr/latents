# What is computed from a fit: the classification tables, three_step() by
# every method, contrast, adjustment and confidence level, r3step(),
# parameter_inference() by every covariance type, compare_models(),
# bootstrap_lrt(), sensitivity() and the enumeration grids. References: each
# statistic recomputed here from the fit's posteriors by its published
# formula, the step-three likelihoods maximized by stats::optim, standard
# errors from numDeriv's Hessian and Jacobian of the likelihood rebuilt from
# the result tables, and stats::p.adjust.

sweep_module_postfit <- function() {
  continuous <- sweep_continuous()
  with_outcomes <- function(data, seed) {
    sweep_seeded(seed, function() {
      score <- as.numeric(scale(data[[1L]]))
      cbind(data, outcome = 0.8 * score + stats::rnorm(nrow(data)),
            covariate = 0.6 * score + stats::rnorm(nrow(data)))
    })
  }
  sets <- list(srl = with_outcomes(continuous$srl, 11L), blobs = with_outcomes(continuous$blobs, 12L),
               faithful = with_outcomes(continuous$faithful, 13L))
  indicators <- function(data) setdiff(names(data), c("outcome", "covariate"))
  fit <- function(data, k, ...) {
    lpa(data, indicators(data), n_profiles = k, seed = 111L, n_starts = 10L, tol = 1e-10,
        max_iter = 20000L, ...)
  }
  frame_of <- function(x) as.data.frame(unclass(x), stringsAsFactors = FALSE)
  settings <- sweep_grid(set = names(sets), k = 2:3)
  each_setting <- function(name, code) {
    lapply(settings, function(s) {
      sweep_define(sprintf("%s[data=%s,k=%d]", name, s$set, s$k),
                   function() code(sets[[s$set]], s$k))
    })
  }
  #' Classification error matrix P(assigned | true) from posteriors and modal classes
  error_matrix <- function(posteriors, modal) {
    K <- ncol(posteriors)
    t(vapply(seq_len(K), function(true) {
      vapply(seq_len(K), function(assigned) {
        sum(posteriors[modal == assigned, true]) / sum(posteriors[, true])
      }, numeric(1))
    }, numeric(K)))
  }

  classification_cases <- each_setting("classification", function(data, k) {
    model <- fit(data, k)
    posteriors <- sweep_posteriors(model)
    modal <- max.col(posteriors, ties.method = "first")
    n <- nrow(posteriors)
    entropy <- sweep_table(model, "entropy")
    average <- sweep_table(model, "average_posteriors")
    errors <- sweep_table(model, "classification_errors")
    classes <- sweep_table(model, "classification")
    weights <- sweep_table(model, "bch_weights")
    expected_average <- vapply(seq_len(k), function(class) {
      vapply(seq_len(k), function(assigned) mean(posteriors[modal == assigned, class]), numeric(1))
    }, numeric(k))
    expected_errors <- error_matrix(posteriors, modal)
    share <- colMeans(posteriors)
    correct <- diag(expected_average)
    inverse <- solve(expected_errors)
    rbind(
      sweep_check("relative entropy: 1 - sum(-p log p) / (n log K)",
                  entropy$relative_entropy[entropy$level == "individuals"],
                  1 - sum(-posteriors * log(pmax(posteriors, .Machine$double.xmin))) / (n * log(k)),
                  1e-10),
      sweep_check("average posterior of each class among the rows assigned to each class",
                  sweep_wide(average, "assigned_class", "class", "average_posterior",
                             seq_len(k), seq_len(k)), expected_average, 1e-10),
      sweep_check("classification errors P(assigned | true)",
                  sweep_wide(errors, "true_class", "assigned_class", "probability",
                             seq_len(k), seq_len(k)), expected_errors, 1e-10),
      sweep_check("modal counts", classes$n_modal, tabulate(modal, k), 0),
      sweep_check("estimated proportions are the mean posteriors", classes$estimated_proportion,
                  share, 1e-10),
      sweep_check("odds of correct classification",
                  classes$odds_correct_classification,
                  (correct / (1 - correct)) / (share / (1 - share)), 1e-8),
      # A = P D, so P is estimated by A D^-1: a unit assigned s carries row s of D^-1.
      sweep_check("BCH weights are the rows of the inverse error matrix",
                  sweep_wide(weights, "unit", "class", "weight", seq_len(n), seq_len(k)),
                  inverse[modal, , drop = FALSE], 1e-8))
  })

  three_step_cases <- each_setting("three-step", function(data, k) {
    model <- fit(data, k)
    posteriors <- sweep_posteriors(model)
    modal <- max.col(posteriors, ties.method = "first")
    y <- data$outcome
    bch <- solve(error_matrix(posteriors, modal))[modal, , drop = FALSE]
    indicator <- diag(k)[modal, , drop = FALSE]
    expected <- list(bch = colSums(bch * y) / colSums(bch),
                     proportional = colSums(posteriors * y) / colSums(posteriors),
                     modal = colSums(indicator * y) / colSums(indicator))
    means <- do.call(rbind, lapply(names(expected), function(method) {
      table <- frame_of(three_step(model, data, "outcome", method = method))
      pairs <- frame_of(three_step(model, data, "outcome", method = method, contrast = "pairs"))
      rbind(
        sweep_check(sprintf("%s: class means of the outcome", method),
                    table$estimate[order(table$class)], expected[[method]], 1e-8),
        sweep_check(sprintf("%s: pairwise differences are differences of the class means", method),
                    pairs$estimate,
                    expected[[method]][pairs$class] - expected[[method]][pairs$reference_class],
                    1e-8),
        sweep_check(sprintf("%s: z statistics and two-sided p-values", method), pairs$p_value,
                    2 * stats::pnorm(-abs(pairs$estimate / pairs$standard_error)), 1e-10))
    }))
    levels <- do.call(rbind, lapply(c(0.8, 0.9, 0.99), function(level) {
      table <- frame_of(three_step(model, data, "outcome", ci_level = level))
      sweep_check(sprintf("confidence limits at level %.2f", level),
                  c(table$conf_low, table$conf_high),
                  c(table$estimate - stats::qnorm((1 + level) / 2) * table$standard_error,
                    table$estimate + stats::qnorm((1 + level) / 2) * table$standard_error), 1e-10)
    }))
    adjustments <- do.call(rbind, lapply(
      c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none"), function(method) {
        table <- frame_of(three_step(model, data, "outcome", contrast = "pairs", adjust = method))
        sweep_check(sprintf("adjust = %s is stats::p.adjust", method), table$p_adjusted,
                    stats::p.adjust(table$p_value, method), 1e-12)
      }))
    overall <- frame_of(three_step(model, data, "outcome", contrast = "overall"))
    # The maximum-likelihood step three: class means and variances of the
    # outcome, and class shares, with the classification errors held fixed.
    errors <- error_matrix(posteriors, modal)
    objective <- function(theta) {
      mean <- theta[seq_len(k)]
      spread <- exp(theta[k + seq_len(k)])
      share <- exp(c(0, theta[2L * k + seq_len(k - 1L)]))
      share <- share / sum(share)
      density <- vapply(seq_len(k), function(class) {
        share[[class]] * errors[class, modal] * stats::dnorm(y, mean[[class]], spread[[class]])
      }, numeric(length(y)))
      -sum(log(rowSums(density)))
    }
    start <- c(expected$bch, rep(log(stats::sd(y)), k), rep(0, k - 1L))
    optimum <- stats::optim(start, objective, method = "BFGS",
                            control = list(reltol = 1e-14, maxit = 2000L))
    ml <- frame_of(three_step(model, data, "outcome", method = "ml"))
    rbind(means, levels, adjustments,
          sweep_check("overall test: degrees of freedom", overall$df, k - 1L, 0),
          sweep_check("overall test: p-value is the chi-square tail", overall$p_value,
                      stats::pchisq(overall$statistic, overall$df, lower.tail = FALSE), 1e-10),
          sweep_check("ml: class means against stats::optim on the step-three likelihood",
                      ml$estimate[order(ml$class)], optimum$par[seq_len(k)], 2e-3))
  })

  r3step_cases <- each_setting("r3step", function(data, k) {
    model <- fit(data, k)
    posteriors <- sweep_posteriors(model)
    modal <- max.col(posteriors, ties.method = "first")
    errors <- error_matrix(posteriors, modal)
    design <- cbind(1, data$covariate)
    table <- frame_of(r3step(model, data, covariates = "covariate", vcov_type = "observed"))
    robust <- frame_of(r3step(model, data, covariates = "covariate", vcov_type = "robust"))
    free <- sort(unique(as.integer(sub("^.*_", "", table$outcome))))
    reference <- setdiff(seq_len(k), free)
    #' Row-wise log likelihood of the logits of the free classes against the reference
    contributions <- function(theta) {
      scores <- matrix(0, nrow(design), k)
      scores[, free] <- design %*% matrix(theta, 2L)
      prior <- exp(scores - sweep_log_sum_exp(scores))
      log(rowSums(prior * t(errors[, modal, drop = FALSE])))
    }
    optimum <- stats::optim(rep(0, 2L * length(free)), function(theta) -sum(contributions(theta)),
                            method = "BFGS", control = list(reltol = 1e-14, maxit = 2000L))
    ordered <- table[order(as.integer(sub("^.*_", "", table$outcome)), table$term != "(Intercept)"), ]
    robust <- robust[order(as.integer(sub("^.*_", "", robust$outcome)),
                           robust$term != "(Intercept)"), ]
    hessian <- numDeriv::hessian(function(theta) -sum(contributions(theta)), optimum$par)
    scores <- numDeriv::jacobian(contributions, optimum$par)
    bread <- solve(hessian)
    units <- nrow(design)
    sandwich <- bread %*% crossprod(scores) %*% bread * units / (units - 1)
    rbind(
      sweep_true("one class is the reference", length(reference) == 1L),
      sweep_check("logits against stats::optim on the step-three likelihood", ordered$estimate,
                  optimum$par, 2e-3),
      sweep_check("observed-information standard errors against numDeriv's Hessian",
                  ordered$standard_error, sqrt(diag(bread)), 5e-3),
      sweep_check("robust standard errors against the sandwich with n / (n - 1)",
                  robust$standard_error, sqrt(diag(sandwich)), 5e-3),
      sweep_check("the covariance type keeps the estimates", robust$estimate, ordered$estimate,
                  1e-10),
      sweep_check("p-values are two-sided normal tails", ordered$p_value,
                  2 * stats::pnorm(-abs(ordered$estimate / ordered$standard_error)), 1e-10))
  })

  inference_cases <- lapply(
    sweep_grid(set = c("faithful", "blobs"), k = 1:2, model = c("VVI", "EEI")), function(s) {
      sweep_define(sprintf("parameter-inference[data=%s,k=%d,model=%s]", s$set, s$k, s$model),
                   function() {
        data <- sets[[s$set]][indicators(sets[[s$set]])]
        model <- lpa(data, names(data), n_profiles = s$k, model = s$model, seed = 111L,
                     n_starts = 10L, tol = 1e-12, max_iter = 20000L)
        x <- as.matrix(data)
        n <- nrow(x)
        p <- ncol(x)
        parameters <- sweep_parameters(model)
        equal <- identical(s$model, "EEI")
        variance_rows <- if (equal) 1L else s$k
        variances <- t(vapply(parameters$covariances, diag, numeric(p)))[seq_len(variance_rows), ,
                                                                       drop = FALSE]
        estimate <- c(as.numeric(parameters$means), as.numeric(variances),
                      log(parameters$mixing[1L, -1L] / parameters$mixing[1L, 1L]))
        #' Row-wise log likelihood in means, variances and share logits
        contributions <- function(theta) {
          means <- matrix(theta[seq_len(s$k * p)], s$k, p)
          spread <- matrix(theta[s$k * p + seq_len(variance_rows * p)], variance_rows, p)
          share <- exp(c(0, theta[-seq_len(s$k * p + variance_rows * p)]))
          share <- share / sum(share)
          density <- vapply(seq_len(s$k), function(class) {
            variance <- spread[if (equal) 1L else class, ]
            log(share[[class]]) + rowSums(stats::dnorm(
              x, matrix(means[class, ], n, p, byrow = TRUE),
              matrix(sqrt(variance), n, p, byrow = TRUE), log = TRUE))
          }, numeric(n))
          sweep_log_sum_exp(matrix(density, n, s$k))
        }
        hessian <- numDeriv::hessian(function(theta) -sum(contributions(theta)), estimate)
        scores <- numDeriv::jacobian(contributions, estimate)
        bread <- solve(hessian)
        expected <- list(observed = bread,
                         robust = bread %*% crossprod(scores) %*% bread * n / (n - 1),
                         opg = solve(crossprod(scores)))
        natural <- seq_len(s$k * p + variance_rows * p)
        rows <- do.call(rbind, lapply(names(expected), function(type) {
          table <- as.data.frame(parameter_inference(model, vcov_type = type))
          measured <- table[table$parameter %in% c("mean", "variance"), , drop = FALSE]
          means <- measured[measured$parameter == "mean", , drop = FALSE]
          spread <- measured[measured$parameter == "variance", , drop = FALSE]
          spread <- spread[!duplicated(if (equal) spread$term else
            paste(spread$outcome, spread$term)), , drop = FALSE]
          ordered <- rbind(
            means[order(match(means$term, colnames(x)), means$outcome), , drop = FALSE],
            spread[order(match(spread$term, colnames(x)), spread$outcome), , drop = FALSE])
          rbind(
            sweep_check(sprintf("%s: estimates are the fit's", type), ordered$estimate,
                        estimate[natural], 1e-8),
            sweep_check(sprintf("%s: standard errors of means and variances against numDeriv",
                                type), ordered$standard_error,
                        sqrt(diag(expected[[type]]))[natural], 1e-2))
        }))
        table <- as.data.frame(parameter_inference(model, level = 0.9, adjust = "holm"))
        tested <- table[!is.na(table$p_value), , drop = FALSE]
        rbind(rows,
              sweep_check("level = 0.9: confidence limits", c(table$conf_low, table$conf_high)[
                !is.na(c(table$conf_low, table$conf_high)) & rep(table$parameter == "mean", 2L)],
                with(table[table$parameter == "mean", ],
                     c(estimate - stats::qnorm(0.95) * standard_error,
                       estimate + stats::qnorm(0.95) * standard_error)), 1e-8),
              if (s$k == 1L && !equal) rbind(
                sweep_check("one profile: standard error of a mean is sqrt(variance / n)",
                            table$standard_error[table$parameter == "mean"],
                            sqrt(diag(parameters$covariances[[1L]]) / n), 1e-4),
                sweep_check("one profile: standard error of a variance is variance * sqrt(2 / n)",
                            table$standard_error[table$parameter == "variance"],
                            diag(parameters$covariances[[1L]]) * sqrt(2 / n), 1e-4)))
      })
    })

  comparison_cases <- lapply(names(sets), function(name) {
    sweep_define(sprintf("compare-models[data=%s]", name), function() {
      data <- sets[[name]]
      fits <- lapply(1:3, function(k) fit(data, k))
      table <- frame_of(do.call(compare_models, fits))
      n <- nrow(data)
      weights <- exp(-0.5 * (table$bic_individual - min(table$bic_individual)))
      grid <- enumerate_lpa(data, indicators(data), n_profiles = 1:3, model = c("VVI", "EEI"),
                            seed = 111L, n_starts = 10L, tol = 1e-10, max_iter = 20000L)
      candidates <- frame_of(as.data.frame(grid))
      direct <- vapply(seq_len(nrow(candidates)), function(i) {
        lpa(data, indicators(data), n_profiles = candidates$n_profiles[[i]],
            model = candidates$structure[[i]] %||% candidates$model[[i]], seed = 111L,
            n_starts = 10L, tol = 1e-10, max_iter = 20000L)$log_likelihood
      }, numeric(1))
      rbind(
        sweep_check("log likelihoods are the fits'", table$log_likelihood,
                    vapply(fits, function(f) f$log_likelihood, numeric(1)), 0),
        sweep_check("AIC", table$aic, -2 * table$log_likelihood + 2 * table$n_parameters, 1e-10),
        sweep_check("BIC", table$bic_individual,
                    -2 * table$log_likelihood + table$n_parameters * log(n), 1e-10),
        sweep_check("sample-size adjusted BIC", table$sabic_individual,
                    -2 * table$log_likelihood + table$n_parameters * log((n + 2) / 24), 1e-10),
        sweep_check("CAIC", table$caic_individual,
                    -2 * table$log_likelihood + table$n_parameters * (log(n) + 1), 1e-10),
        sweep_check("BIC differences", table$delta_bic,
                    table$bic_individual - min(table$bic_individual), 1e-10),
        sweep_check("BIC weights", table$bic_weight, weights / sum(weights), 1e-10),
        sweep_check("enumerate_lpa(): every candidate is the direct fit",
                    candidates$log_likelihood, direct, 1e-9),
        sweep_check("enumerate_lpa(): one candidate per number of profiles and structure",
                    nrow(candidates), 6L, 0))
    })
  })

  resampling_cases <- lapply(names(sets), function(name) {
    sweep_define(sprintf("bootstrap-and-sensitivity[data=%s]", name), function() {
      data <- sets[[name]]
      null <- fit(data, 1L)
      alternative <- fit(data, 2L)
      test <- bootstrap_lrt(null, alternative, iter = 19L, seed = 5L)
      again <- bootstrap_lrt(null, alternative, iter = 19L, seed = 5L)
      replicates <- test$replicates$statistic[test$replicates$valid]
      seeds <- sensitivity(alternative, seeds = 1:4)
      direct <- vapply(1:4, function(seed) {
        lpa(data, indicators(data), n_profiles = 2L, seed = seed, n_starts = 10L, tol = 1e-10,
            max_iter = 20000L)$log_likelihood
      }, numeric(1))
      rbind(
        sweep_check("bootstrap_lrt(): statistic is twice the likelihood difference",
                    test$statistic, 2 * (alternative$log_likelihood - null$log_likelihood), 1e-10),
        sweep_check("bootstrap_lrt(): p-value is (1 + exceedances) / (1 + replicates)",
                    test$p_value, (1 + sum(replicates >= test$statistic)) / (1 + length(replicates)),
                    1e-12),
        sweep_check("bootstrap_lrt(): one seed, one result", again$replicates$statistic,
                    test$replicates$statistic, 0),
        sweep_check("sensitivity(): each seed's log likelihood is the direct fit's",
                    as.data.frame(seeds)$log_likelihood, direct, 1e-9))
    })
  })

  refusal_cases <- list(sweep_define("refusals", function() {
    data <- sets$srl
    model <- fit(data, 2L)
    rbind(
      sweep_refuses("three_step(): an unknown outcome",
                    function() three_step(model, data, "nope")),
      sweep_refuses("three_step(): an unknown method",
                    function() three_step(model, data, "outcome", method = "lta")),
      sweep_refuses("three_step(): an unknown contrast",
                    function() three_step(model, data, "outcome", contrast = "all")),
      sweep_refuses("three_step(): an unknown adjustment",
                    function() three_step(model, data, "outcome", contrast = "pairs",
                                          adjust = "sidak")),
      sweep_refuses("three_step(): ci_level = 1",
                    function() three_step(model, data, "outcome", ci_level = 1)),
      sweep_refuses("three_step(): data of another length",
                    function() three_step(model, data[1:50, ], "outcome")),
      sweep_refuses("three_step(): a missing outcome value",
                    function() three_step(model, transform(data, outcome = replace(outcome, 1L, NA)),
                                          "outcome")),
      sweep_refuses("r3step(): no covariate", function() r3step(model, data)),
      sweep_refuses("r3step(): an unknown covariate",
                    function() r3step(model, data, covariates = "nope")),
      sweep_refuses("r3step(): an unknown vcov_type",
                    function() r3step(model, data, covariates = "covariate", vcov_type = "opg")),
      sweep_refuses("parameter_inference(): level = 0",
                    function() parameter_inference(model, level = 0)),
      sweep_refuses("parameter_inference(): an unknown vcov_type",
                    function() parameter_inference(model, vcov_type = "hc3")),
      sweep_refuses("parameter_inference(): an unknown method",
                    function() parameter_inference(model, method = "jackknife")),
      sweep_refuses("bootstrap_lrt(): models in the wrong order",
                    function() bootstrap_lrt(model, fit(data, 1L), iter = 3L, seed = 1L)),
      sweep_refuses("bootstrap_lrt(): iter = 0",
                    function() bootstrap_lrt(fit(data, 1L), model, iter = 0L)),
      sweep_refuses("compare_models(): one fit", function() compare_models(model)),
      sweep_refuses("compare_models(): not a fit", function() compare_models(model, data)),
      sweep_refuses("sensitivity(): no seeds", function() sensitivity(model, seeds = integer())),
      sweep_refuses("get_candidate(): a candidate outside the grid",
                    function() get_candidate(enumerate_lpa(data, indicators(data), n_profiles = 1:2,
                                                           seed = 1L, n_starts = 2L),
                                             n_profiles = 5L)))
  }))

  c(classification_cases, three_step_cases, r3step_cases, inference_cases, comparison_cases,
    resampling_cases, refusal_cases)
}
