# Membership covariates of multilpa(), lpa() and lca() (`profile_covariates`,
# `group_covariates`, `profile_slopes`), `direct_effects`, and the staged
# estimators (`start` with `fixed`, fit_staged()). References: poLCA's latent
# class regression continued from the package's estimates, the likelihood
# with multinomial-logit shares rebuilt from the result tables, lm() for a
# direct effect at one profile, and nesting relations.

sweep_module_covariates <- function() {
  two_level <- sweep_seeded(909L, function() {
    data <- sweep_two_level()
    unit <- stats::rnorm(length(unique(data$id)))
    cbind(data, x = 0.5 * data$y1 + stats::rnorm(nrow(data)), z = unit[data$id],
          arm = factor(sample(c("a", "b"), nrow(data), TRUE)))
  })
  vars <- c("y1", "y2", "y3")
  items <- sweep_seeded(910L, function() {
    data <- sweep_categorical()$nominal
    cbind(data, x = as.integer(data$item1 == "high") + stats::rnorm(nrow(data)))
  })
  item_names <- paste0("item", 1:5)
  controls <- list(seed = 121L, n_starts = 10L, tol = 1e-10, max_iter = 20000L)
  ll_of <- function(x) sweep_table(x, "information_criteria")$log_likelihood
  par_of <- function(x) sweep_table(x, "information_criteria")$n_parameters
  index_of <- function(label) as.integer(sub("^.*_", "", label))

  #' The likelihood of a covariate fit rebuilt from its tables
  #' @param model The fit; `data` its data; `group` the unit of each row or `NULL`.
  rebuild <- function(model, data, group = NULL) {
    parameters <- sweep_parameters(model)
    log_density <- sweep_log_density(data, parameters)
    coefficients <- sweep_table(model, "coefficients")
    K <- model$n_profiles
    G <- model$n_group_classes
    n <- nrow(data)
    profile_rows <- coefficients[coefficients$level == "profile", , drop = FALSE]
    #' Linear predictor of each free profile in group class g
    scores <- function(g) {
      out <- matrix(0, n, K)
      free <- unique(profile_rows$outcome)
      out[, index_of(free)] <- vapply(free, function(outcome) {
        rows <- profile_rows[profile_rows$outcome == outcome, , drop = FALSE]
        intercept <- rows$estimate[rows$term %in% c("(Intercept)", sprintf("group_class_%d", g))]
        slopes <- rows[!rows$term %in% c("(Intercept)", sprintf("group_class_%d", seq_len(G))), ,
                       drop = FALSE]
        design <- stats::model.matrix(stats::reformulate(unique(slopes$term), intercept = FALSE),
                                      data)
        intercept + as.numeric(design[, slopes$term, drop = FALSE] %*% slopes$estimate)
      }, numeric(n))
      out
    }
    row_given_class <- vapply(seq_len(G), function(g) {
      linear <- scores(g)
      sweep_log_sum_exp(linear - sweep_log_sum_exp(linear) + log_density)
    }, numeric(n))
    row_given_class <- matrix(row_given_class, n, G)
    if (is.null(group)) return(sum(row_given_class))
    unit <- match(group, unique(group))
    group_rows <- coefficients[coefficients$level == "group", , drop = FALSE]
    first <- data[!duplicated(unit), , drop = FALSE]
    class_scores <- matrix(0, nrow(first), G)
    free <- unique(group_rows$outcome)
    class_scores[, index_of(free)] <- vapply(free, function(outcome) {
      rows <- group_rows[group_rows$outcome == outcome, , drop = FALSE]
      slopes <- rows[rows$term != "(Intercept)", , drop = FALSE]
      rows$estimate[rows$term == "(Intercept)"] +
        if (nrow(slopes)) as.numeric(as.matrix(first[slopes$term]) %*% slopes$estimate) else
          numeric(nrow(first))
    }, numeric(nrow(first)))
    joint <- rowsum(row_given_class, unit) + class_scores - sweep_log_sum_exp(class_scores)
    sum(sweep_log_sum_exp(joint))
  }
  fit_lpa <- function(...) do.call(lpa, c(list(two_level, vars), list(...), controls))
  fit_two <- function(...) do.call(multilpa, c(list(two_level, vars, id = "id"), list(...),
                                               controls))

  single_cases <- lapply(sweep_grid(k = 2:3, covariates = c("x", "x+z")), function(option) {
    sweep_define(sprintf("profile-covariates[k=%d,%s]", option$k, option$covariates), function() {
      covariates <- strsplit(option$covariates, "+", fixed = TRUE)[[1L]]
      model <- fit_lpa(n_profiles = option$k, profile_covariates = covariates)
      plain <- fit_lpa(n_profiles = option$k)
      rescaled <- do.call(lpa, c(list(transform(two_level, x = 100 * x + 5), vars,
                                      n_profiles = option$k, profile_covariates = covariates),
                                 controls))
      rbind(
        sweep_check("log likelihood rebuilt from the result tables", ll_of(model),
                    rebuild(model, two_level), 1e-8),
        sweep_check("parameter count", par_of(model),
                    par_of(plain) + (option$k - 1L) * length(covariates), 0),
        sweep_at_least("covariates at least the likelihood without them", ll_of(model),
                       ll_of(plain), 1e-7),
        sweep_check("a rescaled covariate keeps the likelihood", ll_of(rescaled), ll_of(model),
                    1e-7))
    })
  })

  two_level_cases <- lapply(
    sweep_grid(k = 2L, g = 2L, design = c("profile", "group", "both", "slopes-by-class")),
    function(option) {
      sweep_define(sprintf("two-level-covariates[%s]", option$design), function() {
        arguments <- switch(option$design,
          profile = list(profile_covariates = "x"),
          group = list(group_covariates = "z"),
          both = list(profile_covariates = "x", group_covariates = "z"),
          "slopes-by-class" = list(profile_covariates = "x", profile_slopes = "group_class"))
        model <- do.call(fit_two, c(list(n_profiles = option$k, n_group_classes = option$g),
                                    arguments))
        plain <- fit_two(n_profiles = option$k, n_group_classes = option$g)
        added <- switch(option$design, profile = 1L, group = 1L, both = 2L, "slopes-by-class" = 2L)
        rows <- rbind(
          sweep_check("parameter count", par_of(model), par_of(plain) + added, 0),
          sweep_at_least("covariates at least the likelihood without them", ll_of(model),
                         ll_of(plain), 1e-6))
        if (identical(option$design, "slopes-by-class")) {
          shared <- fit_two(n_profiles = option$k, n_group_classes = option$g,
                            profile_covariates = "x")
          return(rbind(rows, sweep_at_least("class-specific slopes at least the shared slope",
                                            ll_of(model), ll_of(shared), 1e-6)))
        }
        rbind(rows, sweep_check("log likelihood rebuilt from the result tables", ll_of(model),
                                rebuild(model, two_level, two_level$id), 1e-8))
      })
    })

  polca_cases <- lapply(2:3, function(k) {
    sweep_define(sprintf("class-regression-against-poLCA[k=%d]", k), function() {
      model <- do.call(lca, c(list(items, item_names, n_classes = k, profile_covariates = "x"),
                              controls))
      codes <- cbind(as.data.frame(lapply(items[item_names], as.integer)), x = items$x)
      formula <- stats::as.formula(sprintf("cbind(%s) ~ x", paste(item_names, collapse = ", ")))
      responses <- sweep_table(model, "responses")
      start <- lapply(item_names, function(v) {
        unname(sweep_wide(responses[responses$indicator == v, , drop = FALSE], "profile",
                          "category", "probability", seq_len(k), levels(items[[v]])))
      })
      continued <- suppressWarnings(poLCA::poLCA(formula, codes, nclass = k, probs.start = start,
                                                 maxiter = 50000L, tol = 1e-12, verbose = FALSE,
                                                 calc.se = FALSE))
      searched <- sweep_seeded(3L, function() {
        suppressWarnings(poLCA::poLCA(formula, codes, nclass = k, nrep = 15L, maxiter = 50000L,
                                      tol = 1e-12, verbose = FALSE, calc.se = FALSE))
      })
      rbind(
        sweep_check("log likelihood rebuilt from the result tables", ll_of(model),
                    rebuild(model, items), 1e-8),
        sweep_check("poLCA's regression continued from the package's probabilities",
                    ll_of(model), continued$llik, 1e-6),
        sweep_at_least("optimum at least that of poLCA's search", ll_of(model), searched$llik,
                       1e-6),
        sweep_check("parameter count is poLCA's", par_of(model), searched$npar, 0))
    })
  })

  staged_cases <- list(
    sweep_define("two-step-holds-the-measurement", function() {
      first <- fit_lpa(n_profiles = 2L)
      second <- fit_lpa(n_profiles = 2L, profile_covariates = "x", start = starting_values(first),
                        fixed = "measurement")
      joint <- fit_lpa(n_profiles = 2L, profile_covariates = "x")
      rbind(
        sweep_check("held means are the first step's", sweep_table(second, "profiles")$mean,
                    sweep_table(first, "profiles")$mean, 1e-12),
        sweep_check("held variances are the first step's", sweep_table(second, "profiles")$variance,
                    sweep_table(first, "profiles")$variance, 1e-12),
        sweep_check("only the logits are counted", par_of(second), 2L, 0),
        sweep_check("log likelihood rebuilt from the result tables", ll_of(second),
                    rebuild(second, two_level), 1e-8),
        sweep_at_least("the joint fit at least the two-step likelihood", ll_of(joint),
                       ll_of(second), 1e-7))
    }),
    sweep_define("fit_staged-is-the-two-stage-fit", function() {
      staged <- do.call(fit_staged, c(list(two_level, vars, id = "id", n_profiles = 2L,
                                           n_group_classes = 2L), controls))
      first <- fit_two(n_profiles = 2L, n_group_classes = 1L)
      manual <- fit_two(n_profiles = 2L, n_group_classes = 2L, start = starting_values(first, what = "measurement"),
                        fixed = "measurement")
      rbind(
        sweep_check("log likelihood of the manual two-stage fit", ll_of(staged), ll_of(manual),
                    1e-7),
        sweep_check("the measurement is the one-class fit's", sweep_table(staged, "profiles")$mean,
                    sweep_table(first, "profiles")$mean, 1e-10),
        sweep_at_least("the joint fit at least the staged likelihood",
                       ll_of(fit_two(n_profiles = 2L, n_group_classes = 2L)), ll_of(staged), 1e-7))
    }),
    sweep_define("direct-effect-at-one-profile-is-lm", function() {
      model <- fit_lpa(n_profiles = 1L, direct_effects = "y1 ~ x")
      reference <- stats::lm(y1 ~ x, two_level)
      effect <- sweep_table(model, "direct_effects")
      profile <- sweep_table(model, "profiles")
      rbind(
        sweep_check("the direct effect is lm's slope", effect$estimate, stats::coef(reference)[["x"]],
                    1e-7),
        sweep_check("the mean at covariate zero is lm's intercept",
                    profile$mean[profile$indicator == "y1"], stats::coef(reference)[[1L]], 1e-7),
        sweep_check("the variance is lm's ML residual variance",
                    profile$variance[profile$indicator == "y1"],
                    mean(stats::residuals(reference)^2), 1e-7),
        sweep_check("one more parameter than without the effect", par_of(model),
                    par_of(fit_lpa(n_profiles = 1L)) + 1L, 0))
    }),
    sweep_define("direct-effect-nesting", function() {
      model <- fit_lpa(n_profiles = 2L, direct_effects = c("y1 ~ x", "y2 ~ x"))
      one <- fit_lpa(n_profiles = 2L, direct_effects = "y1 ~ x")
      rbind(sweep_at_least("two effects at least one", ll_of(model), ll_of(one), 1e-6),
            sweep_at_least("one effect at least none", ll_of(one), ll_of(fit_lpa(n_profiles = 2L)),
                           1e-6),
            sweep_check("parameter count", par_of(model), par_of(one) + 1L, 0))
    }))

  refusal_cases <- list(sweep_define("refusals", function() {
    holes <- transform(two_level, x = replace(x, 5L, NA))
    rbind(
      sweep_refuses("a missing covariate value",
                    function() lpa(holes, vars, 2L, profile_covariates = "x")),
      sweep_refuses("a covariate that is an indicator",
                    function() lpa(two_level, vars, 2L, profile_covariates = "y1")),
      sweep_refuses("a constant covariate",
                    function() lpa(transform(two_level, x = 1), vars, 2L,
                                   profile_covariates = "x")),
      sweep_refuses("a group covariate varying within groups",
                    function() multilpa(two_level, vars, id = "id", 2L, 2L,
                                        group_covariates = "x")),
      sweep_refuses("group covariates in a single-level model",
                    function() lpa(two_level, vars, 2L, group_covariates = "z")),
      sweep_refuses("an unknown profile_slopes option",
                    function() multilpa(two_level, vars, id = "id", 2L, 2L,
                                        profile_covariates = "x", profile_slopes = "free")),
      sweep_refuses("covariates with a prior",
                    function() lpa(two_level, vars, 2L, profile_covariates = "x",
                                   prior = prior_control())),
      sweep_refuses("fixed without start", function() lpa(two_level, vars, 2L, fixed = "means")),
      sweep_refuses("an unknown fixed block",
                    function() lpa(two_level, vars, 2L, fixed = "shares",
                                   start = starting_values(lpa(two_level, vars, 2L, seed = 1L)))),
      sweep_refuses("a direct effect on an unknown indicator",
                    function() lpa(two_level, vars, 2L, direct_effects = "nope ~ x")),
      sweep_refuses("a direct effect with full covariance",
                    function() lpa(two_level, vars, 2L, direct_effects = "y1 ~ x",
                                   covariance_model = "full")),
      sweep_refuses("a malformed direct effect",
                    function() lpa(two_level, vars, 2L, direct_effects = "y1 x")))
  }))

  c(single_cases, two_level_cases, polca_cases, staged_cases, refusal_cases)
}
