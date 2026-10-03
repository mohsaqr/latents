# Section: mixture regression and growth mixture models
#
# mixture_regression() in its single-level, group-level and two-level forms,
# the three families, shared and membership terms, weights, missing data and
# every stored covariance type; mixture_regression(random = ) growth mixture
# models with each random-effect covariance; enumerate_regressions() and
# compare_models().

#' study_hours with the derived columns the regression cases need
#'
#' `weight_row` is a deterministic per-row sampling weight, `weight_student`
#' one constant within each student, `successes`/`failures` turn `questions`
#' into binomial counts out of eight, `exposure` is a Poisson exposure, and
#' `score_gappy` has every 17th score missing.
#' @return A data frame.
golden_reg_hours <- function() {
  data <- golden_hours()
  data$weight_row <- 0.5 + (seq_len(nrow(data)) %% 4L) / 4
  data$weight_student <- 0.5 + (data$student %% 3L) / 2
  data$successes <- pmin(data$questions, 8L)
  data$failures <- 8L - data$successes
  data$exposure <- data$week + 1
  data$score_gappy <- data$score
  data$score_gappy[seq(5L, nrow(data), by = 17L)] <- NA_real_
  data
}

#' growth_scores with a deterministic per-student weight
#' @return A data frame.
golden_reg_growth <- function() {
  data <- golden_growth()
  data$weight_student <- 0.5 + (data$student %% 3L) / 2
  data
}

golden_cases_regression <- function() {
  hours <- golden_reg_hours()
  growth <- golden_reg_growth()
  # Ordinal outcomes for the ordinal cases only, so the other cases' data are
  # untouched: `grade` (an ordered factor from score) and `level` (questions
  # capped at five, whole-number categories 0..5).
  ordinal_hours <- golden_reg_hours()
  ordinal_hours$grade <- cut(ordinal_hours$score, c(-Inf, 50, 65, 80, Inf),
                             labels = c("D", "C", "B", "A"),
                             ordered_result = TRUE)
  ordinal_hours$level <- pmin(ordinal_hours$questions, 5L)
  regression <- function(...) {
    mixture_regression(score ~ hours, data = hours, n_classes = 2,
                       n_starts = 2, seed = 1, ...)
  }
  growth_fit <- function(...) {
    mixture_regression(score ~ wave, data = growth, n_classes = 2,
                       id = "student", class_level = "group", n_starts = 2,
                       seed = 1, ...)
  }
  list(
    # -- mixture_regression(), gaussian -------------------------------------
    golden_case(
      "reg_gaussian_single", "mixture_regression",
      "gaussian, single level, 2 classes, varying variance, observed vcov; predict",
      fit = function() regression(),
      data = hours, services = "predict"),
    golden_case(
      "reg_gaussian_equal_variance", "mixture_regression",
      "gaussian, variance = equal",
      fit = function() regression(variance = "equal"),
      data = hours),
    golden_case(
      "reg_gaussian_common", "mixture_regression",
      "gaussian, common = ~ sleep (coefficient shared across classes)",
      fit = function() {
        mixture_regression(score ~ hours + sleep, data = hours, n_classes = 2,
                           common = ~ sleep, n_starts = 2, seed = 1)
      },
      data = hours, services = "predict"),
    golden_case(
      "reg_gaussian_membership", "mixture_regression",
      "gaussian, row-level membership = ~ sleep (concomitant logit)",
      fit = function() regression(membership = ~ sleep),
      data = hours, services = "predict"),
    golden_case(
      "reg_gaussian_select_converged", "mixture_regression",
      "gaussian, select_start = converged",
      fit = function() regression(select_start = "converged"),
      data = hours),
    golden_case(
      "reg_gaussian_vcov_robust_clustered", "mixture_regression",
      "gaussian, observation classes with id and one group class: robust vcov clustered on id",
      fit = function() regression(id = "student", vcov_type = "robust"),
      data = hours),
    golden_case(
      "reg_gaussian_vcov_opg", "mixture_regression",
      "gaussian, vcov_type = opg",
      fit = function() regression(vcov_type = "opg"),
      data = hours),
    golden_case(
      "reg_gaussian_vcov_none", "mixture_regression",
      "gaussian, vcov_type = none (no inference stored)",
      fit = function() regression(vcov_type = "none"),
      data = hours),
    golden_case(
      "reg_gaussian_weights_rows", "mixture_regression",
      "gaussian, single-level sampling weights per row (robust default)",
      fit = function() regression(weights = "weight_row"),
      data = hours),
    golden_case(
      "reg_gaussian_weights_observed_refused", "mixture_regression",
      "weights with vcov_type = observed is refused",
      fit = function() regression(weights = "weight_row", vcov_type = "observed"),
      data = hours),
    golden_case(
      "reg_gaussian_missing_omit", "mixture_regression",
      "missing = omit drops rows with NA outcome (latents_rows_dropped warning)",
      fit = function() {
        mixture_regression(score_gappy ~ hours, data = hours, n_classes = 2,
                           missing = "omit", n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_refusal_missing_error", "mixture_regression",
      "missing = error with NA outcome is refused (latents_missing_data)",
      fit = function() {
        mixture_regression(score_gappy ~ hours, data = hours, n_classes = 2,
                           n_starts = 2, seed = 1)
      },
      data = hours),
    # -- class_level = "group" and two-level --------------------------------
    golden_case(
      "reg_group_level", "mixture_regression",
      "gaussian, class_level = group (one class per student)",
      fit = function() regression(id = "student", class_level = "group"),
      data = hours, services = "predict"),
    golden_case(
      "reg_group_membership", "mixture_regression",
      "class_level = group with group-constant membership = ~ motivation",
      fit = function() {
        regression(id = "student", class_level = "group",
                   membership = ~ motivation)
      },
      data = hours),
    golden_case(
      "reg_group_weights", "mixture_regression",
      "class_level = group with per-student sampling weights",
      fit = function() {
        regression(id = "student", class_level = "group",
                   weights = "weight_student")
      },
      data = hours),
    golden_case(
      "reg_two_level", "mixture_regression",
      "two-level: observation classes, id, n_group_classes = 2",
      fit = function() regression(id = "student", n_group_classes = 2),
      data = hours, services = "predict"),
    golden_case(
      "reg_two_level_group_membership", "mixture_regression",
      "two-level with group_membership = ~ motivation and variance = equal",
      fit = function() {
        regression(id = "student", n_group_classes = 2,
                   group_membership = ~ motivation, variance = "equal")
      },
      data = hours),
    golden_case(
      "reg_two_level_robust", "mixture_regression",
      "two-level, vcov_type = robust",
      fit = function() {
        regression(id = "student", n_group_classes = 2, vcov_type = "robust")
      },
      data = hours),
    # -- binomial and poisson -----------------------------------------------
    golden_case(
      "reg_binomial_group", "mixture_regression",
      "binomial (0/1 outcome), class_level = group",
      fit = function() {
        mixture_regression(passed ~ hours, data = hours, n_classes = 2,
                           family = "binomial", id = "student",
                           class_level = "group", n_starts = 2, seed = 1)
      },
      data = hours, services = "predict"),
    golden_case(
      "reg_binomial_cbind", "mixture_regression",
      "binomial with cbind(successes, failures), single level",
      fit = function() {
        mixture_regression(cbind(successes, failures) ~ hours, data = hours,
                           n_classes = 2, family = "binomial", n_starts = 2,
                           seed = 1)
      },
      data = hours, services = "predict"),
    golden_case(
      "reg_poisson", "mixture_regression",
      "poisson, single level",
      fit = function() {
        mixture_regression(questions ~ hours, data = hours, n_classes = 2,
                           family = "poisson", n_starts = 2, seed = 1)
      },
      data = hours, services = "predict"),
    golden_case(
      "reg_poisson_offset", "mixture_regression",
      "poisson with offset(log(exposure)) and class_level = group",
      fit = function() {
        mixture_regression(questions ~ hours + offset(log(exposure)),
                           data = hours, n_classes = 2, family = "poisson",
                           id = "student", class_level = "group",
                           n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_negative_binomial", "mixture_regression",
      "negative_binomial, single level, dispersion varying (phase E7)",
      fit = function() {
        mixture_regression(questions ~ hours, data = hours, n_classes = 2,
                           family = "negative_binomial", n_starts = 2, seed = 1)
      },
      data = hours, services = c("predict", "simulate")),
    golden_case(
      "reg_negative_binomial_group", "mixture_regression",
      "negative_binomial with offset, class_level = group, dispersion equal (phase E7)",
      fit = function() {
        mixture_regression(questions ~ hours + offset(log(exposure)),
                           data = hours, n_classes = 2,
                           family = "negative_binomial", variance = "equal",
                           id = "student", class_level = "group",
                           n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_ordinal_group", "mixture_regression",
      "ordinal (ordered factor grade), class_level = group (phase E7)",
      fit = function() {
        mixture_regression(grade ~ hours, data = ordinal_hours, n_classes = 2,
                           family = "ordinal", id = "student",
                           class_level = "group", n_starts = 2, seed = 1)
      },
      data = ordinal_hours, services = c("predict", "simulate")),
    golden_case(
      "reg_ordinal_common", "mixture_regression",
      "ordinal (integer categories), single level, common = ~ sleep (phase E7)",
      fit = function() {
        mixture_regression(level ~ hours + sleep, data = ordinal_hours,
                           n_classes = 2, family = "ordinal", common = ~ sleep,
                           n_starts = 2, seed = 1)
      },
      data = ordinal_hours, services = c("predict", "simulate")),
    golden_case(
      "reg_ordinal_two_level", "mixture_regression",
      "ordinal (integer categories), two-level, n_group_classes = 2 (phase E7)",
      fit = function() {
        mixture_regression(level ~ hours, data = ordinal_hours, n_classes = 2,
                           family = "ordinal", id = "student",
                           n_group_classes = 2, n_starts = 2, seed = 1)
      },
      data = ordinal_hours),
    golden_case(
      "reg_ordinal_refusal_two_categories", "mixture_regression",
      "ordinal with two categories and one outcome per row is refused (latents_not_identified)",
      fit = function() {
        mixture_regression(factor(passed) ~ hours, data = ordinal_hours,
                           n_classes = 2, family = "ordinal", n_starts = 2,
                           seed = 1)
      },
      data = ordinal_hours),
    golden_case(
      "reg_refusal_binary_rows", "mixture_regression",
      "binary outcome with one trial per class assignment is refused",
      fit = function() {
        mixture_regression(passed ~ hours, data = hours, n_classes = 2,
                           family = "binomial", n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_refusal_group_without_id", "mixture_regression",
      "class_level = group without id is refused",
      fit = function() regression(class_level = "group"),
      data = hours),
    golden_case(
      "reg_refusal_poisson_noninteger", "mixture_regression",
      "poisson on a non-count outcome is refused",
      fit = function() {
        mixture_regression(score ~ hours, data = hours, n_classes = 2,
                           family = "poisson", n_starts = 2, seed = 1)
      },
      data = hours),
    # -- enumeration and comparison -----------------------------------------
    golden_case(
      "reg_enumerate", "enumerate_regressions",
      "enumerate_regressions over 1:2 classes, no bootstrap",
      fit = function() {
        enumerate_regressions(score ~ hours, data = hours, n_classes = 1:2,
                              n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_enumerate_bootstrap", "enumerate_regressions",
      "enumerate_regressions with a 3-replicate BLRT (bootstrap_starts = 2)",
      fit = function() {
        enumerate_regressions(score ~ hours, data = hours, n_classes = 1:2,
                              bootstrap = 3L, bootstrap_starts = 2L,
                              n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "reg_compare_models", "compare_models",
      "compare_models over 1-, 2- and 3-class single-level regressions",
      fit = function() {
        compare_models(
          one = mixture_regression(score ~ hours, data = hours, n_classes = 1,
                                   n_starts = 2, seed = 1),
          two = regression(),
          three = mixture_regression(score ~ hours, data = hours,
                                     n_classes = 3, n_starts = 2, seed = 1))
      },
      data = hours),
    # -- growth mixture models ----------------------------------------------
    golden_case(
      "growth_lcga", "mixture_regression",
      "latent class growth analysis: score ~ wave, class_level = group, no random",
      fit = function() growth_fit(),
      data = growth, services = "predict"),
    golden_case(
      "growth_random_intercept", "growth_mixture",
      "growth mixture, random = intercept, covariance varying",
      fit = function() growth_fit(random = "intercept"),
      data = growth, services = "predict"),
    golden_case(
      "growth_random_slope", "growth_mixture",
      "growth mixture, random intercept + slope (random = wave), varying; compare_models vs LCGA",
      fit = function() growth_fit(random = "wave"),
      data = growth, services = "compare_models",
      other_fits = list(function() growth_fit())),
    golden_case(
      "growth_random_slope_only", "growth_mixture",
      "growth mixture, random = ~ 0 + wave (slope without random intercept)",
      fit = function() growth_fit(random = ~ 0 + wave),
      data = growth),
    golden_case(
      "growth_covariance_equal", "growth_mixture",
      "growth mixture, random = wave, random_covariance = equal",
      fit = function() growth_fit(random = "wave", random_covariance = "equal"),
      data = growth),
    golden_case(
      "growth_covariance_proportional", "growth_mixture",
      "growth mixture, random = wave, random_covariance = proportional",
      fit = function() {
        growth_fit(random = "wave", random_covariance = "proportional")
      },
      data = growth),
    golden_case(
      "growth_random_diagonal", "growth_mixture",
      "growth mixture, random = wave, random_diagonal = TRUE",
      fit = function() growth_fit(random = "wave", random_diagonal = TRUE),
      data = growth),
    golden_case(
      "growth_variance_equal", "growth_mixture",
      "growth mixture, random intercept, residual variance = equal",
      fit = function() growth_fit(random = "intercept", variance = "equal"),
      data = growth),
    golden_case(
      "growth_membership", "growth_mixture",
      "growth mixture, random intercept, membership = ~ motivation",
      fit = function() {
        growth_fit(random = "intercept", membership = ~ motivation)
      },
      data = growth),
    golden_case(
      "growth_weights", "growth_mixture",
      "growth mixture, random intercept, per-student sampling weights",
      fit = function() {
        growth_fit(random = "intercept", weights = "weight_student")
      },
      data = growth),
    golden_case(
      "growth_vcov_robust", "growth_mixture",
      "growth mixture, random = wave, vcov_type = robust",
      fit = function() growth_fit(random = "wave", vcov_type = "robust"),
      data = growth),
    golden_case(
      "growth_vcov_opg", "growth_mixture",
      "growth mixture, random = wave, vcov_type = opg",
      fit = function() growth_fit(random = "wave", vcov_type = "opg"),
      data = growth),
    golden_case(
      "growth_refusal_poisson_random", "growth_mixture",
      "random effects with a non-gaussian family are refused",
      fit = function() {
        mixture_regression(round(score) ~ wave, data = growth, n_classes = 2,
                           family = "poisson", id = "student",
                           class_level = "group", random = "intercept",
                           n_starts = 2, seed = 1)
      },
      data = growth),
    golden_case(
      "growth_refusal_observation_level", "growth_mixture",
      "random effects with class_level = observation are refused",
      fit = function() {
        mixture_regression(score ~ wave, data = growth, n_classes = 2,
                           id = "student", random = "intercept",
                           n_starts = 2, seed = 1)
      },
      data = growth),
    # -- multilevel growth mixture models (phase E7) --------------------------
    golden_case(
      "growth_multilevel", "growth_mixture",
      "multilevel growth mixture: growth_schools, cluster = school, 2 group classes, group_membership = programme",
      fit = function() {
        mixture_regression(score ~ wave, growth_schools, n_classes = 2,
                           id = "student", class_level = "group", random = "wave",
                           random_covariance = "equal", cluster = "school",
                           n_group_classes = 2, group_membership = "programme",
                           n_starts = 2, seed = 1)
      },
      data = growth_schools, services = "simulate"),
    golden_case(
      "growth_multilevel_one_class", "growth_mixture",
      "multilevel growth mixture with one group class (schools as units)",
      fit = function() {
        mixture_regression(score ~ wave, growth_schools, n_classes = 2,
                           id = "student", class_level = "group", random = "intercept",
                           cluster = "school", n_group_classes = 1,
                           n_starts = 2, seed = 1)
      },
      data = growth_schools),
    golden_case(
      "growth_refusal_cluster_weights", "growth_mixture",
      "a multilevel growth mixture refuses sampling weights",
      fit = function() {
        weighted <- transform(growth_schools, weight = 1)
        mixture_regression(score ~ wave, weighted, n_classes = 2, id = "student",
                           class_level = "group", random = "wave", cluster = "school",
                           n_group_classes = 2, weights = "weight",
                           n_starts = 2, seed = 1)
      },
      data = growth_schools)
  )
}
