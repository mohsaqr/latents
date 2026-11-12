# mixture_regression() with random effects (growth mixture models): `random`,
# `random_covariance`, `random_diagonal`, `residual_correlation`,
# `residual_time`, `variance`, `membership` and `weights`, on the two bundled
# longitudinal datasets. References: lme4's maximum-likelihood mixed model at
# one class, the marginal Gaussian likelihood of each person's trajectory
# rebuilt from the result tables (fixed effects, random-effect covariance,
# residual variance and residual process), and the nesting of the covariance
# options.

sweep_module_growth <- function() {
  sets <- list(
    scores = list(data = latents::growth_scores, formula = score ~ wave, id = "student",
                  time = "wave", covariate = "motivation"),
    hours = list(data = latents::study_hours, formula = score ~ week, id = "student",
                 time = "week", covariate = "motivation"))
  frame_of <- function(x) as.data.frame(unclass(x), stringsAsFactors = FALSE)
  table_of <- function(fit, what) frame_of(suppressMessages(get_results(fit, what)))
  fit <- function(set, k, random, ..., data = set$data, n_starts = 6L) {
    mixture_regression(set$formula, data, n_classes = k, id = set$id, class_level = "group",
                       random = random, seed = 91L, n_starts = n_starts, ...)
  }
  ll_of <- function(x) table_of(x, "fit")$log_likelihood
  par_of <- function(x) table_of(x, "fit")$n_parameters

  #' The marginal likelihood of every person's rows, rebuilt from the tables
  rebuild <- function(model, set, random, data = set$data, membership = NULL, weights = NULL) {
    coefficients <- table_of(model, "coefficients")
    classes <- table_of(model, "classes")
    effects <- table_of(model, "random")
    labels <- classes$class
    design <- stats::model.matrix(set$formula, data)
    random_design <- stats::model.matrix(random, data)
    response <- data[[all.vars(set$formula)[[1L]]]]
    person <- match(data[[set$id]], unique(data[[set$id]]))
    time <- data[[set$time]]
    process <- effects[startsWith(effects$parameter, "process_"), , drop = FALSE]
    #' Random-effect covariance of a class, from its variances and correlations
    covariance_of <- function(label) {
      # A covariance shared by the classes is reported once, under "all".
      rows <- effects[effects$class == label, , drop = FALSE]
      if (!nrow(rows)) {
        rows <- effects[effects$class == "all" & !startsWith(effects$parameter, "process_"), ,
                        drop = FALSE]
      }
      terms <- colnames(random_design)
      spread <- sqrt(rows$estimate[rows$parameter == "variance"][match(
        terms, rows$term[rows$parameter == "variance"])])
      correlation <- diag(length(terms))
      pairs <- rows[rows$parameter == "correlation", , drop = FALSE]
      if (nrow(pairs)) {
        ends <- do.call(rbind, strsplit(sub("^(\\(Intercept\\)|[^:]+):", "\\1\r", pairs$term),
                                        "\r", fixed = TRUE))
        at <- cbind(match(ends[, 1L], terms), match(ends[, 2L], terms))
        correlation[at] <- pairs$estimate
        correlation[at[, 2:1, drop = FALSE]] <- pairs$estimate
      }
      outer(spread, spread) * correlation
    }
    log_density <- vapply(labels, function(label) {
      rows <- coefficients[coefficients$class %in% c(label, "common"), , drop = FALSE]
      mean <- as.numeric(design[, rows$term, drop = FALSE] %*% rows$estimate)
      between <- covariance_of(label)
      residual <- classes$residual_sd[classes$class == label]^2
      vapply(split(seq_along(person), person), function(rows_of) {
        z <- random_design[rows_of, , drop = FALSE]
        covariance <- z %*% between %*% t(z) + diag(residual, length(rows_of))
        if (nrow(process)) {
          spread <- process$estimate[process$parameter == "process_sd"]
          t <- time[rows_of]
          covariance <- covariance + spread^2 * if ("process_rate" %in% process$parameter) {
            exp(-process$estimate[process$parameter == "process_rate"] * abs(outer(t, t, `-`)))
          } else outer(t, t, pmin)
        }
        root <- chol(covariance)
        scaled <- backsolve(root, response[rows_of] - mean[rows_of], transpose = TRUE)
        -0.5 * sum(scaled^2) - sum(log(diag(root))) - 0.5 * length(rows_of) * log(2 * pi)
      }, numeric(1))
    }, numeric(length(unique(person))))
    first <- !duplicated(person)
    prior <- if (is.null(membership)) {
      matrix(classes$share, sum(first), length(labels), byrow = TRUE)
    } else {
      logits <- table_of(model, "membership")
      covariates <- stats::model.matrix(membership, data[first, , drop = FALSE])
      scores <- vapply(labels, function(label) {
        rows <- logits[logits$class == label, , drop = FALSE]
        if (!nrow(rows)) return(numeric(nrow(covariates)))
        as.numeric(covariates[, rows$term, drop = FALSE] %*% rows$estimate)
      }, numeric(nrow(covariates)))
      exp(scores - sweep_log_sum_exp(scores))
    }
    joint <- log_density + log(prior)
    person_log_likelihood <- sweep_log_sum_exp(joint)
    list(log_likelihood = sum((if (is.null(weights)) 1 else weights[first]) *
                                person_log_likelihood),
         posteriors = exp(joint - person_log_likelihood))
  }
  rebuilt_checks <- function(model, set, random, ...) {
    rebuilt <- rebuild(model, set, random, ...)
    summary <- table_of(model, "fit")
    rbind(
      sweep_check("marginal log likelihood rebuilt from the result tables",
                  summary$log_likelihood, rebuilt$log_likelihood, 1e-6),
      sweep_check("AIC from the likelihood and parameter count", summary$aic,
                  -2 * summary$log_likelihood + 2 * summary$n_parameters, 1e-10),
      sweep_check("BIC over persons from the likelihood and parameter count", summary$bic,
                  -2 * summary$log_likelihood + summary$n_parameters * log(summary$n_persons),
                  1e-10))
  }

  one_class_cases <- lapply(
    sweep_grid(set = names(sets), random = c("intercept", "slope", "diagonal")), function(option) {
      sweep_define(sprintf("one-class[data=%s,random=%s]", option$set, option$random), function() {
        set <- sets[[option$set]]
        random <- if (identical(option$random, "intercept")) ~ 1 else
          stats::as.formula(sprintf("~ %s", set$time))
        diagonal <- identical(option$random, "diagonal")
        model <- fit(set, 1L, random, random_diagonal = diagonal)
        outcome <- all.vars(set$formula)[[1L]]
        reference <- lme4::lmer(stats::as.formula(sprintf(
          "%s ~ %s + (%s %s %s)", outcome, set$time,
          if (identical(option$random, "intercept")) "1" else set$time,
          if (diagonal) "||" else "|", set$id)), set$data, REML = FALSE,
          control = lme4::lmerControl(optCtrl = list(xtol_abs = 1e-10, ftol_abs = 1e-12)))
        components <- as.data.frame(lme4::VarCorr(reference))
        effects <- table_of(model, "random")
        variances <- effects[effects$parameter == "variance", , drop = FALSE]
        expected <- components[is.na(components$var2) & components$grp != "Residual", , drop = FALSE]
        rbind(
          rebuilt_checks(model, set, random),
          sweep_check("lme4's maximum log likelihood", ll_of(model),
                      as.numeric(stats::logLik(reference)), 1e-6),
          sweep_check("lme4's fixed effects", table_of(model, "coefficients")$estimate,
                      lme4::fixef(reference), 1e-4),
          sweep_check("lme4's random-effect variances",
                      variances$estimate[match(expected$var1, variances$term)], expected$vcov,
                      2e-3),
          sweep_check("lme4's residual standard deviation", table_of(model, "classes")$residual_sd,
                      stats::sigma(reference), 1e-4),
          sweep_check("parameter count", par_of(model),
                      2L + 1L + if (identical(option$random, "intercept")) 1L else
                        if (diagonal) 2L else 3L, 0))
      })
    })

  class_cases <- lapply(
    sweep_grid(set = names(sets), k = 2:3,
               structure = c("varying", "equal", "proportional", "diagonal", "intercept")),
    function(option) {
      sweep_define(sprintf("classes[data=%s,k=%d,%s]", option$set, option$k, option$structure),
                   function() {
        set <- sets[[option$set]]
        slope <- stats::as.formula(sprintf("~ %s", set$time))
        random <- if (identical(option$structure, "intercept")) ~ 1 else slope
        arguments <- switch(option$structure,
          varying = list(), intercept = list(),
          equal = list(random_covariance = "equal"),
          proportional = list(random_covariance = "proportional"),
          diagonal = list(random_diagonal = TRUE))
        model <- do.call(fit, c(list(set, option$k, random), arguments))
        rows <- rebuilt_checks(model, set, random)
        if (identical(option$structure, "varying")) return(rows)
        # Every structure here restricts the varying one fitted from the same starts' model.
        rbind(rows, sweep_at_least("varying covariances at least this structure's likelihood",
                                   ll_of(fit(set, option$k, slope)), ll_of(model), 1e-5))
      })
    })

  option_cases <- c(
    lapply(sweep_grid(set = names(sets), process = c("ar", "brownian")), function(option) {
      sweep_define(sprintf("residual-process[data=%s,%s]", option$set, option$process),
                   function() {
        set <- sets[[option$set]]
        random <- stats::as.formula(sprintf("~ %s", set$time))
        model <- fit(set, 2L, random, residual_correlation = option$process)
        named <- fit(set, 2L, random, residual_correlation = option$process,
                     residual_time = set$time)
        rbind(rebuilt_checks(model, set, random),
              sweep_check("residual_time named against its default", ll_of(named), ll_of(model), 0),
              sweep_at_least("the process at least the independent-residual likelihood",
                             ll_of(model), ll_of(fit(set, 2L, random)), 1e-6))
      })
    }),
    lapply(names(sets), function(name) {
      sweep_define(sprintf("membership[data=%s]", name), function() {
        set <- sets[[name]]
        random <- stats::as.formula(sprintf("~ %s", set$time))
        membership <- stats::as.formula(sprintf("~ %s", set$covariate))
        model <- fit(set, 2L, random, membership = membership)
        rbind(rebuilt_checks(model, set, random, membership = membership),
              sweep_at_least("covariates at least the intercept-only likelihood", ll_of(model),
                             ll_of(fit(set, 2L, random)), 1e-6))
      })
    }),
    lapply(names(sets), function(name) {
      sweep_define(sprintf("invariances[data=%s]", name), function() {
        set <- sets[[name]]
        random <- stats::as.formula(sprintf("~ %s", set$time))
        base <- fit(set, 2L, random)
        shuffled <- sweep_seeded(5L, function() set$data[sample(nrow(set$data)), , drop = FALSE])
        labelled <- set$data
        labelled[[set$id]] <- paste0("p", labelled[[set$id]])
        rbind(
          sweep_check("two runs with one seed: coefficients",
                      table_of(fit(set, 2L, random), "coefficients")$estimate,
                      table_of(base, "coefficients")$estimate, 0),
          sweep_check("rows shuffled: log likelihood", ll_of(fit(set, 2L, random, data = shuffled)),
                      ll_of(base), 1e-7),
          sweep_check("character identifier: log likelihood",
                      ll_of(fit(set, 2L, random, data = labelled)), ll_of(base), 1e-9),
          sweep_check("random effects by name against a formula: log likelihood",
                      ll_of(fit(set, 2L, set$time)), ll_of(base), 0),
          sweep_check("weights of one: log likelihood",
                      ll_of(fit(set, 2L, random, data = cbind(set$data, w = 1), weights = "w")),
                      ll_of(base), 1e-8),
          sweep_check("vcov_type = robust keeps the coefficients",
                      table_of(fit(set, 2L, random, vcov_type = "robust"), "coefficients")$estimate,
                      table_of(base, "coefficients")$estimate, 0))
      })
    }))

  refusal_cases <- list(sweep_define("refusals", function() {
    set <- sets$scores
    call <- function(...) function() {
      mixture_regression(score ~ wave, set$data, n_classes = 2L, id = "student",
                         class_level = "group", ...)
    }
    rbind(
      sweep_refuses("random without id",
                    function() mixture_regression(score ~ wave, set$data, n_classes = 2L,
                                                  random = ~ wave)),
      sweep_refuses("an unknown random term", call(random = ~ nope)),
      sweep_refuses("an unknown random_covariance", call(random = ~ wave,
                                                         random_covariance = "free")),
      sweep_refuses("random_diagonal = NA", call(random = ~ wave, random_diagonal = NA)),
      sweep_refuses("an unknown residual_correlation",
                    call(random = ~ wave, residual_correlation = "arma")),
      sweep_refuses("an unknown residual_time",
                    call(random = ~ wave, residual_correlation = "ar", residual_time = "nope")),
      sweep_refuses("random effects with a binomial outcome",
                    function() mixture_regression(passed ~ week, latents::study_hours,
                                                  n_classes = 2L, family = "binomial",
                                                  id = "student", class_level = "group",
                                                  random = ~ 1)),
      sweep_refuses("an unknown cluster", call(random = ~ wave, cluster = "nope")))
  }))

  c(one_class_cases, class_cases, option_cases, refusal_cases)
}
