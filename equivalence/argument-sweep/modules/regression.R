# mixture_regression() without random effects: the five families, classes of
# rows and classes of units, shared coefficients (`common`), membership
# covariates, equal or varying residual variances, weights, missing rows and
# predictor types, on the bundled weekly study data, flexmix's NPreg and a
# simulated nested set. References: lm/glm/MASS::glm.nb/MASS::polr at one
# class, flexmix's EM continued from the package's posteriors, and the
# mixture-of-regressions likelihood rebuilt from the result tables.

sweep_module_regression <- function() {
  suppressPackageStartupMessages(library(flexmix))
  environment <- new.env()
  utils::data("NPreg", package = "flexmix", envir = environment)
  study <- transform(latents::study_hours,
                     level = cut(score, c(-Inf, 50, 65, Inf), labels = c("low", "mid", "high"),
                                 ordered_result = TRUE))
  sets <- list(study = study, simulated = sweep_regression(), npreg = environment$NPreg)
  frame_of <- function(x) as.data.frame(unclass(x), stringsAsFactors = FALSE)
  table_of <- function(fit, what) frame_of(suppressMessages(get_results(fit, what)))
  fit <- function(formula, data, k, ..., n_starts = 10L) {
    mixture_regression(formula, data, n_classes = k, seed = 81L, n_starts = n_starts,
                       tol = 1e-10, max_iter = 5000L, ...)
  }
  ll_of <- function(x) table_of(x, "fit")$log_likelihood
  par_of <- function(x) table_of(x, "fit")$n_parameters

  #' Row-by-class log densities and prior probabilities rebuilt from the tables
  rebuild <- function(model, formula, data, family, unit = NULL, membership = NULL,
                      weights = NULL) {
    coefficients <- table_of(model, "coefficients")
    classes <- table_of(model, "classes")
    labels <- classes$class
    frame <- stats::model.frame(formula, data)
    design <- stats::model.matrix(formula, frame)
    response <- stats::model.response(frame)
    is_threshold <- startsWith(coefficients$term, "threshold:")
    log_density <- vapply(labels, function(label) {
      rows <- coefficients[coefficients$class %in% c(label, "common") & !is_threshold, ,
                           drop = FALSE]
      columns <- intersect(colnames(design), rows$term)
      eta <- as.numeric(design[, columns, drop = FALSE] %*%
                          rows$estimate[match(columns, rows$term)])
      of_class <- classes[classes$class == label, , drop = FALSE]
      switch(family,
        gaussian = stats::dnorm(response, eta, of_class$sigma, log = TRUE),
        binomial = stats::dbinom(as.integer(response), 1L, stats::plogis(eta), log = TRUE),
        poisson = stats::dpois(response, exp(eta), log = TRUE),
        negative_binomial = if (of_class$dispersion > 1e-7) {
          stats::dnbinom(response, mu = exp(eta), size = 1 / of_class$dispersion, log = TRUE)
        } else stats::dpois(response, exp(eta), log = TRUE),
        ordinal = {
          cuts <- coefficients$estimate[coefficients$class == label & is_threshold]
          upper <- stats::plogis(outer(-eta, c(cuts, Inf), `+`))
          lower <- stats::plogis(outer(-eta, c(-Inf, cuts), `+`))
          log((upper - lower)[cbind(seq_along(eta), as.integer(response))])
        })
    }, numeric(nrow(design)))
    prior <- if (is.null(membership)) {
      matrix(classes$prior, nrow(design), length(labels), byrow = TRUE)
    } else {
      logits <- table_of(model, "membership")
      if ("model" %in% names(logits)) logits <- logits[logits$model == "class", , drop = FALSE]
      covariates <- stats::model.matrix(membership, data[rownames(frame), , drop = FALSE])
      scores <- vapply(labels, function(label) {
        rows <- logits[logits$class == label, , drop = FALSE]
        if (!nrow(rows)) return(numeric(nrow(covariates)))
        as.numeric(covariates[, rows$term, drop = FALSE] %*% rows$estimate)
      }, numeric(nrow(covariates)))
      exp(scores - sweep_log_sum_exp(scores))
    }
    unit_of_row <- if (is.null(unit)) seq_len(nrow(design)) else
      match(data[rownames(frame), unit], unique(data[rownames(frame), unit]))
    joint <- rowsum(log_density, unit_of_row) + log(prior[!duplicated(unit_of_row), , drop = FALSE])
    unit_log_likelihood <- sweep_log_sum_exp(joint)
    unit_weights <- if (is.null(weights)) 1 else weights[!duplicated(unit_of_row)]
    list(log_likelihood = sum(unit_weights * unit_log_likelihood),
         posteriors = exp(joint - unit_log_likelihood))
  }

  #' Checks of one fit against the rebuilt likelihood
  rebuilt_checks <- function(model, formula, data, family, unit = NULL, ...) {
    rebuilt <- rebuild(model, formula, data, family, unit, ...)
    summary <- table_of(model, "fit")
    reported <- table_of(model, if (is.null(unit)) "assignments" else "groups")
    reported <- as.matrix(reported[grep("^probability_class_", names(reported))])
    n_bic <- if (is.null(unit)) summary$n_observations else summary$n_groups
    rbind(
      sweep_check("log likelihood rebuilt from the result tables", summary$log_likelihood,
                  rebuilt$log_likelihood, 1e-7),
      sweep_check("posteriors rebuilt from the result tables", reported, rebuilt$posteriors, 1e-5),
      sweep_check("AIC from the likelihood and parameter count", summary$aic,
                  -2 * summary$log_likelihood + 2 * summary$n_parameters, 1e-10),
      sweep_check("BIC from the likelihood and parameter count", summary$bic,
                  -2 * summary$log_likelihood + summary$n_parameters * log(n_bic), 1e-10))
  }

  #' flexmix's driver for a family
  flexmix_model <- function(family) {
    switch(family, gaussian = flexmix::FLXMRglm(family = "gaussian"),
           binomial = flexmix::FLXMRglm(family = "binomial"),
           poisson = flexmix::FLXMRglm(family = "poisson"))
  }
  flexmix_formula <- function(formula, family, unit) {
    text <- deparse(formula)
    if (identical(family, "binomial")) {
      response <- all.vars(formula)[[1L]]
      text <- sub(sprintf("^%s", response), sprintf("cbind(%s, 1 - %s)", response, response), text)
    }
    stats::as.formula(if (is.null(unit)) text else sprintf("%s | %s", text, unit))
  }

  specifications <- list(
    list(set = "study", formula = score ~ hours + sleep, family = "gaussian", unit = "student"),
    list(set = "study", formula = passed ~ hours + sleep, family = "binomial", unit = "student"),
    list(set = "study", formula = questions ~ hours, family = "poisson", unit = "student"),
    list(set = "study", formula = questions ~ hours, family = "negative_binomial",
         unit = "student"),
    list(set = "study", formula = level ~ hours, family = "ordinal", unit = "student"),
    list(set = "simulated", formula = y ~ x + arm, family = "gaussian", unit = "unit"),
    list(set = "simulated", formula = b ~ x, family = "binomial", unit = "unit"),
    list(set = "simulated", formula = n ~ x, family = "poisson", unit = "unit"),
    list(set = "simulated", formula = n ~ x, family = "negative_binomial", unit = "unit"),
    list(set = "simulated", formula = grade ~ x, family = "ordinal", unit = "unit"),
    list(set = "npreg", formula = yn ~ x + I(x^2), family = "gaussian", unit = NULL),
    list(set = "npreg", formula = yp ~ x, family = "poisson", unit = NULL))
  label_of <- function(s) sprintf("%s:%s:%s", s$set, s$family, deparse(s$formula))

  one_class_cases <- lapply(specifications, function(s) {
    sweep_define(sprintf("one-class[%s]", label_of(s)), function() {
      data <- sets[[s$set]]
      model <- fit(s$formula, data, 1L, family = s$family)
      coefficients <- table_of(model, "coefficients")
      reference <- switch(s$family,
        gaussian = stats::lm(s$formula, data),
        binomial = stats::glm(s$formula, stats::binomial(), data,
                              control = stats::glm.control(epsilon = 1e-12, maxit = 100L)),
        poisson = stats::glm(s$formula, stats::poisson(), data,
                             control = stats::glm.control(epsilon = 1e-12, maxit = 100L)),
        negative_binomial = MASS::glm.nb(s$formula, data,
                                         control = stats::glm.control(epsilon = 1e-12, maxit = 200L)),
        ordinal = MASS::polr(s$formula, data, method = "logistic"))
      slopes <- coefficients[!startsWith(coefficients$term, "threshold:"), , drop = FALSE]
      expected <- stats::coef(reference)
      rbind(
        rebuilt_checks(model, s$formula, data, s$family),
        sweep_check("log likelihood of the standard regression", ll_of(model),
                    as.numeric(stats::logLik(reference)), 1e-7),
        sweep_check("coefficients of the standard regression",
                    slopes$estimate[match(names(expected), slopes$term)], expected, 1e-4),
        if (identical(s$family, "ordinal")) {
          sweep_check("thresholds of MASS::polr",
                      coefficients$estimate[startsWith(coefficients$term, "threshold:")],
                      reference$zeta, 1e-4)
        },
        if (identical(s$family, "gaussian")) {
          sweep_check("residual standard deviation is the ML one",
                      table_of(model, "classes")$sigma,
                      sqrt(mean(stats::residuals(reference)^2)), 1e-6)
        })
    })
  })

  mixture_cases <- lapply(
    unlist(lapply(specifications, function(s) {
      # One binary outcome per class assignment is refused as not identified
      # (checked among the refusals), so binomial classes are classes of units.
      levels <- c(if (!identical(s$family, "binomial")) "observation",
                  if (!is.null(s$unit)) "group")
      lapply(sweep_grid(k = 2:3, level = levels), function(o) c(s, o))
    }), recursive = FALSE), function(s) {
      sweep_define(sprintf("classes[%s,k=%d,%s]", label_of(s), s$k, s$level), function() {
        data <- sets[[s$set]]
        nested <- identical(s$level, "group")
        model <- if (nested) {
          fit(s$formula, data, s$k, family = s$family, id = s$unit, class_level = "group")
        } else fit(s$formula, data, s$k, family = s$family)
        rows <- rebuilt_checks(model, s$formula, data, s$family, if (nested) s$unit)
        if (!s$family %in% c("gaussian", "binomial", "poisson")) return(rows)
        posteriors <- rebuild(model, s$formula, data, s$family, if (nested) s$unit)$posteriors
        start <- if (nested) {
          posteriors[match(data[[s$unit]], unique(data[[s$unit]])), , drop = FALSE]
        } else posteriors
        continued <- suppressWarnings(flexmix::flexmix(
          flexmix_formula(s$formula, s$family, if (nested) s$unit), data = data, k = s$k,
          model = flexmix_model(s$family), cluster = start,
          control = list(tolerance = 1e-12, iter.max = 20000L, minprior = 0)))
        # flexmix's Gaussian M-step divides the residual sum of squares by the
        # weighted degrees of freedom, not by the weight, so its fixed point is
        # not the maximum likelihood estimate and its likelihood is lower.
        if (identical(s$family, "gaussian")) {
          return(rbind(rows, sweep_at_least(
            "flexmix's EM from the package's posteriors does not improve the likelihood",
            ll_of(model), as.numeric(logLik(continued)), 1e-8)))
        }
        rbind(rows, sweep_check("flexmix's EM from the package's posteriors: log likelihood",
                                ll_of(model), as.numeric(logLik(continued)), 1e-6))
      })
    })

  option_cases <- list(
    sweep_define("equal-variance", function() {
      data <- sets$simulated
      equal <- fit(y ~ x + arm, data, 2L, variance = "equal")
      varying <- fit(y ~ x + arm, data, 2L)
      rbind(rebuilt_checks(equal, y ~ x + arm, data, "gaussian"),
            sweep_check("one residual standard deviation", diff(table_of(equal, "classes")$sigma),
                        0, 1e-10),
            sweep_check("one parameter fewer", par_of(equal), par_of(varying) - 1L, 0),
            sweep_at_least("varying variances at least the equal-variance likelihood",
                           ll_of(varying), ll_of(equal), 1e-7))
    }),
    sweep_define("common-coefficients", function() {
      data <- sets$simulated
      shared <- fit(y ~ x + arm, data, 2L, common = ~ arm)
      by_name <- fit(y ~ x + arm, data, 2L, common = "arm")
      free <- fit(y ~ x + arm, data, 2L)
      rbind(rebuilt_checks(shared, y ~ x + arm, data, "gaussian"),
            sweep_check("one parameter fewer", par_of(shared), par_of(free) - 1L, 0),
            sweep_check("a name against a formula", ll_of(by_name), ll_of(shared), 0),
            sweep_at_least("free coefficients at least the shared-coefficient likelihood",
                           ll_of(free), ll_of(shared), 1e-7))
    }),
    sweep_define("membership-covariates[observation]", function() {
      data <- sets$simulated
      model <- fit(y ~ x + arm, data, 2L, membership = ~ z)
      by_name <- fit(y ~ x + arm, data, 2L, membership = "z")
      continued <- suppressWarnings(flexmix::flexmix(
        y ~ x + arm, data = data, k = 2L, concomitant = flexmix::FLXPmultinom(~ z),
        cluster = rebuild(model, y ~ x + arm, data, "gaussian", membership = ~ z)$posteriors,
        control = list(tolerance = 1e-12, iter.max = 20000L, minprior = 0)))
      rbind(rebuilt_checks(model, y ~ x + arm, data, "gaussian", membership = ~ z),
            sweep_check("a name against a formula", ll_of(by_name), ll_of(model), 0),
            sweep_at_least("flexmix with a concomitant model, continued, does not improve",
                           ll_of(model), as.numeric(logLik(continued)), 1e-8),
            sweep_at_least("covariates at least the intercept-only likelihood", ll_of(model),
                           ll_of(fit(y ~ x + arm, data, 2L)), 1e-7))
    }),
    sweep_define("membership-covariates[group]", function() {
      data <- sets$simulated
      model <- fit(y ~ x + arm, data, 2L, id = "unit", class_level = "group", membership = ~ z)
      rbind(rebuilt_checks(model, y ~ x + arm, data, "gaussian", "unit", membership = ~ z),
            sweep_at_least("covariates at least the intercept-only likelihood", ll_of(model),
                           ll_of(fit(y ~ x + arm, data, 2L, id = "unit", class_level = "group")),
                           1e-7))
    }),
    sweep_define("vcov-type-keeps-the-estimates", function() {
      data <- sets$simulated
      base <- fit(y ~ x + arm, data, 2L)
      do.call(rbind, lapply(c("robust", "opg", "none"), function(type) {
        sweep_check(sprintf("vcov_type = %s: coefficients", type),
                    table_of(fit(y ~ x + arm, data, 2L, vcov_type = type), "coefficients")$estimate,
                    table_of(base, "coefficients")$estimate, 0)
      }))
    }),
    sweep_define("missing-rows", function() {
      data <- sweep_punch(sets$simulated, 0.08, columns = c("x", "y"))
      complete <- data[stats::complete.cases(data[c("x", "y", "arm")]), , drop = FALSE]
      rbind(sweep_check("missing = omit against the complete rows: log likelihood",
                        ll_of(fit(y ~ x + arm, data, 2L, missing = "omit")),
                        ll_of(fit(y ~ x + arm, complete, 2L)), 1e-9),
            sweep_refuses("missing rows without missing = omit",
                          function() fit(y ~ x + arm, data, 2L)))
    }),
    sweep_define("weights", function() {
      data <- sets$simulated
      w <- sweep_seeded(9L, function() sample(1:3, nrow(data), TRUE))
      weighted <- fit(y ~ x + arm, cbind(data, w = w), 2L, weights = "w")
      expanded <- fit(y ~ x + arm, data[rep(seq_len(nrow(data)), w), ], 2L)
      unit <- fit(y ~ x + arm, cbind(data, w = 1), 2L, weights = "w")
      rbind(
        sweep_check("weights of one: log likelihood", ll_of(unit),
                    ll_of(fit(y ~ x + arm, data, 2L)), 1e-9),
        sweep_check("integer weights against repeated rows: sorted coefficients",
                    sort(table_of(weighted, "coefficients")$estimate),
                    sort(table_of(expanded, "coefficients")$estimate), 1e-3),
        sweep_check("integer weights against repeated rows: pseudo log likelihood",
                    ll_of(weighted), nrow(data) / sum(w) * ll_of(expanded), 1e-6))
    }),
    sweep_define("predictor-types", function() {
      data <- sets$simulated
      base <- fit(y ~ x + arm, data, 2L)
      rbind(
        sweep_check("a character predictor: log likelihood",
                    ll_of(fit(y ~ x + arm, transform(data, arm = as.character(arm)), 2L)),
                    ll_of(base), 1e-10),
        sweep_check("a logical predictor: log likelihood",
                    ll_of(fit(y ~ x + arm, transform(data, arm = arm == "treated"), 2L)),
                    ll_of(base), 1e-10),
        sweep_check("a rescaled predictor: log likelihood",
                    ll_of(fit(y ~ x + arm, transform(data, x = 100 * x + 3), 2L)), ll_of(base),
                    1e-7),
        sweep_check("rows shuffled: log likelihood",
                    ll_of(fit(y ~ x + arm, sweep_seeded(5L, function() data[sample(nrow(data)), ]),
                              2L)), ll_of(base), 1e-7),
        sweep_check("two runs with one seed: coefficients",
                    table_of(fit(y ~ x + arm, data, 2L), "coefficients")$estimate,
                    table_of(base, "coefficients")$estimate, 0),
        sweep_check("select_start = converged: log likelihood",
                    ll_of(fit(y ~ x + arm, data, 2L, select_start = "converged")), ll_of(base),
                    1e-8))
    }),
    sweep_define("enumerate_regressions", function() {
      data <- sets$simulated
      grid <- enumerate_regressions(y ~ x + arm, data, n_classes = 1:3, seed = 81L,
                                    n_starts = 10L, tol = 1e-10, max_iter = 5000L)
      table <- frame_of(as.data.frame(grid))
      direct <- vapply(1:3, function(k) ll_of(fit(y ~ x + arm, data, k)), numeric(1))
      sweep_check("the grid's log likelihoods are those of the direct fits",
                  table$log_likelihood[match(1:3, table$n_classes)], direct, 1e-9)
    }))

  refusal_cases <- list(sweep_define("refusals", function() {
    data <- sets$simulated
    rbind(
      sweep_refuses("n_classes = 0", function() mixture_regression(y ~ x, data, n_classes = 0L)),
      sweep_refuses("an unknown predictor",
                    function() mixture_regression(y ~ nope, data, n_classes = 2L)),
      sweep_refuses("a formula without a response",
                    function() mixture_regression(~ x, data, n_classes = 2L)),
      sweep_refuses("an unknown family",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, family = "gamma")),
      sweep_refuses("a non-binary outcome under binomial",
                    function() mixture_regression(y ~ x, data, n_classes = 2L,
                                                  family = "binomial")),
      sweep_refuses("a negative outcome under poisson",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, family = "poisson")),
      sweep_refuses("one binary outcome per class assignment",
                    function() mixture_regression(b ~ x, data, n_classes = 2L,
                                                  family = "binomial")),
      sweep_refuses("class_level = group without id",
                    function() mixture_regression(y ~ x, data, n_classes = 2L,
                                                  class_level = "group")),
      sweep_refuses("an unknown id",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, id = "nope",
                                                  class_level = "group")),
      sweep_refuses("a shared intercept",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, common = ~ 1)),
      sweep_refuses("a shared term outside the formula",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, common = ~ arm)),
      sweep_refuses("a membership covariate varying within units",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, id = "unit",
                                                  class_level = "group", membership = ~ x)),
      sweep_refuses("variance = equal under poisson",
                    function() mixture_regression(n ~ x, data, n_classes = 2L, family = "poisson",
                                                  variance = "equal")),
      sweep_refuses("an unknown weights column",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, weights = "w")),
      sweep_refuses("residual_correlation without random",
                    function() mixture_regression(y ~ x, data, n_classes = 2L,
                                                  residual_correlation = "ar")),
      sweep_refuses("cluster without random",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, cluster = "unit")),
      sweep_refuses("an unknown vcov_type",
                    function() mixture_regression(y ~ x, data, n_classes = 2L, vcov_type = "hc3")))
  }))

  c(one_class_cases, mixture_cases, option_cases, refusal_cases)
}
