# Section: membership-covariate models and the Houle group-class families
# (additive, dispersion, additive-dispersion, and the two cross-level
# families), with their comparisons and refusals.

#' Engagement data with the covariates, weights and missing cells the cases use
#'
#' Sixty students rather than the usual twenty: with twenty, the group
#' covariate separates the group classes perfectly (extreme coefficients and
#' a singular information matrix) and the full-covariance model has as many
#' parameters as groups, so no Wald path would be recorded.
#'
#' `previous_grade` is a row covariate; `grade_mean` is its student mean, a
#' group covariate constant within each student; `w` is a sampling weight
#' constant within each student; `forum_band` is an ordinal (1-3) version of
#' `forum_post` and `posts` a count derived from it.
#' @param n Number of students.
#' @return A data frame.
golden_cf_engagement <- function(n = 60L) {
  data <- golden_engagement(n)
  data$grade_mean <- stats::ave(data$previous_grade, data$student)
  data$w <- 1 + (data$student %% 3L)
  data$forum_band <- as.integer(cut(data$forum_post, c(-Inf, -0.5, 0.5, Inf)))
  data$posts <- as.integer(round(exp(data$forum_post)))
  data
}

#' The same data with a deterministic pattern of missing indicator cells
#' @return A data frame.
golden_cf_engagement_missing <- function() {
  data <- golden_cf_engagement()
  rows <- seq_len(nrow(data))
  data$browse[rows %% 7L == 0L] <- NA
  data$lectures[rows %% 11L == 3L] <- NA
  data
}

#' Simulated additive-family data with interior estimates
#'
#' On `course_engagement` the varying-between-variance families put a
#' between-group variance on zero, so every Wald probe is refused at the
#' boundary. These data (60 groups of 8, two well-separated classes, positive
#' between variances; the draw of tests/testthat/helper-additive-reference.R)
#' give interior fits, so the standard-error paths are recorded too. The
#' caller's random state is restored.
#' @return A data frame with `group`, `y1`, `y2` and a group weight `w`.
golden_cf_additive_data <- function() {
  old_seed <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  on.exit(if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old_seed, globalenv())
  }, add = TRUE)
  set.seed(21L)
  means <- rbind(c(-1, -1), c(1, 1))
  between <- rbind(c(0.25, 0.4), c(0.4, 0.25))
  sizes <- rep(8L, 60L)
  classes <- sample.int(2L, length(sizes), replace = TRUE)
  blocks <- lapply(seq_along(sizes), function(j) {
    intercept <- stats::rnorm(2L, means[classes[j], ],
                              sqrt(between[classes[j], ]))
    ratings <- sweep(matrix(stats::rnorm(sizes[j] * 2L), sizes[j], 2L), 2L,
                     intercept, "+")
    colnames(ratings) <- c("y1", "y2")
    data.frame(group = sprintf("g%02d", j), ratings, w = 1 + j %% 3L)
  })
  do.call(rbind, blocks)
}

#' Experience-sampling data with a numeric profile covariate (`day`)
#' @return A data frame.
golden_cf_esm <- function() {
  data <- golden_esm(30L)
  data$day_mean <- stats::ave(data$day, data$student)
  data
}

golden_cases_covariates_families <- function() {
  engagement <- golden_cf_engagement()
  missing_data <- golden_cf_engagement_missing()
  esm <- golden_cf_esm()
  additive_data <- golden_cf_additive_data()
  vars <- golden_engagement_vars
  covariate_services <- c("inference_opg", "diagnostics")

  covariates <- list(
    golden_case(
      "cov_profile_shared", "covariates",
      "two-level, profile_covariates only (row covariate), shared slopes",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement,
      services = c(covariate_services, "inference_bootstrap", "bootstrap_lrt"),
      null_fit = function() multilpa(engagement, vars, id = "student",
                                     n_profiles = 2, n_group_classes = 1,
                                     profile_covariates = "previous_grade",
                                     n_starts = 2, seed = 1)),
    golden_case(
      "cov_group_only", "covariates",
      "two-level, group_covariates only (student-constant covariate)",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                group_covariates = "grade_mean",
                                n_starts = 2, seed = 1),
      data = engagement, services = covariate_services),
    golden_case(
      "cov_both_shared", "covariates",
      "two-level, profile and group covariates, shared slopes",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                profile_covariates = "previous_grade",
                                group_covariates = "grade_mean",
                                n_starts = 2, seed = 1),
      data = engagement, services = covariate_services),
    golden_case(
      "cov_slopes_group_class", "covariates",
      "two-level, profile_slopes = group_class (cross-level interaction)",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                profile_covariates = "previous_grade",
                                profile_slopes = "group_class",
                                n_starts = 2, seed = 1),
      data = engagement, services = "inference_opg"),
    golden_case(
      "cov_one_group_class", "covariates",
      "two-level with one group class, profile covariate (pooled LPA logit)",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 1,
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement, services = covariate_services),
    golden_case(
      "cov_single_level", "covariates",
      "single level (lpa, id = NULL) with a profile covariate",
      fit = function() lpa(engagement, vars, n_profiles = 2,
                           profile_covariates = "previous_grade",
                           n_starts = 2, seed = 1),
      data = engagement, services = covariate_services),
    golden_case(
      "cov_three_profiles_equal", "covariates",
      "three profiles, variance_model equal, profile covariate",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 3,
                                n_group_classes = 1, variance_model = "equal",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_fiml", "covariates",
      "profile covariate with missing indicators, missing = fiml",
      fit = function() multilpa(missing_data, vars, id = "student",
                                n_profiles = 2, n_group_classes = 2,
                                profile_covariates = "previous_grade",
                                missing = "fiml", n_starts = 2, seed = 1),
      data = missing_data, services = "inference_opg"),
    golden_case(
      "cov_full_covariance", "covariates",
      "profile covariate, covariance_model = full",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 1, covariance_model = "full",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement, services = "inference_opg"),
    golden_case(
      "cov_categorical_esm", "covariates",
      "two continuous + binary categorical (ESM) with a row covariate",
      fit = function() multilpa(esm, c("happy", "relaxed", "time_with_friends"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 1,
                                categorical = "time_with_friends",
                                profile_covariates = "day",
                                n_starts = 2, seed = 1),
      data = esm, services = c("inference_opg", "inference_bootstrap")),
    golden_case(
      "cov_categorical_boundary", "covariates",
      "categorical indicator at its probability bound with a covariate (Wald refused)",
      fit = function() multilpa(engagement, c("browse", "lectures", "engagement"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 1,
                                categorical = "engagement",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement, services = "inference_opg"),
    golden_case(
      "cov_ordinal", "covariates",
      "ordinal indicator (adjacent-category logit) with a profile covariate",
      fit = function() multilpa(engagement, c("browse", "lectures", "forum_band"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 1, ordinal = "forum_band",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_count_poisson", "covariates",
      "Poisson count indicator with a profile covariate",
      fit = function() multilpa(engagement, c("browse", "lectures", "posts"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 1, count = "posts",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_esm_mixed", "covariates",
      "ESM: gaussian + ordinal + binary categorical, row covariate day, two levels",
      fit = function() multilpa(esm, c("happy", "relaxed", "worried",
                                       "time_with_friends"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 2, ordinal = "worried",
                                categorical = "time_with_friends",
                                profile_covariates = "day",
                                group_covariates = "day_mean",
                                n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "cov_weights", "covariates",
      "sampling weights (constant within student) with profile covariate",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                profile_covariates = "previous_grade",
                                weights = "w", n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_refusal_start", "covariates",
      "refusal: start = with membership covariates",
      fit = function() {
        base <- multilpa(engagement, vars, id = "student", n_profiles = 2,
                         n_group_classes = 1, n_starts = 2, seed = 1)
        multilpa(engagement, vars, id = "student", n_profiles = 2,
                 n_group_classes = 1, profile_covariates = "previous_grade",
                 start = starting_values(base), n_starts = 2, seed = 1)
      },
      data = engagement),
    golden_case(
      "cov_refusal_fixed", "covariates",
      "refusal: fixed measurement with membership covariates",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 1, fixed = "means",
                                profile_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_refusal_missing_covariate", "covariates",
      "refusal: missing values in a membership covariate",
      fit = function() {
        data <- engagement
        data$previous_grade[c(3L, 40L)] <- NA
        multilpa(data, vars, id = "student", n_profiles = 2,
                 n_group_classes = 1, profile_covariates = "previous_grade",
                 n_starts = 2, seed = 1)
      },
      data = engagement),
    golden_case(
      "cov_refusal_slopes_without_covariates", "covariates",
      "refusal: profile_slopes = group_class without profile_covariates",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                profile_slopes = "group_class",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "cov_refusal_group_covariate_varies", "covariates",
      "refusal: a group covariate that varies within a group",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                group_covariates = "previous_grade",
                                n_starts = 2, seed = 1),
      data = engagement)
  )

  family_services <- c("inference_opg", "inference_bootstrap")
  family_fit <- function(family, between_variance = "varying",
                         n_group_classes = 2L, data = additive_data,
                         id = "group", indicators = c("y1", "y2"), ...) {
    force(family)
    force(between_variance)
    force(n_group_classes)
    force(data)
    force(id)
    force(indicators)
    extra <- list(...)
    function() {
      do.call(multilpa, c(list(data, indicators, id = id,
                               n_group_classes = n_group_classes,
                               family = family,
                               between_variance = between_variance,
                               n_starts = 2, seed = 1), extra))
    }
  }
  cross_fit <- function(family, between_variance = "varying",
                        variance_model = "varying") {
    force(family)
    force(between_variance)
    force(variance_model)
    function() {
      multilpa(engagement, vars, id = "student", n_profiles = 2,
               n_group_classes = 2, family = family,
               between_variance = between_variance,
               variance_model = variance_model, n_starts = 2, seed = 1)
    }
  }

  families <- list(
    golden_case("family_additive_varying", "additive",
                "family additive, between_variance varying, 2 group classes (simulated, interior)",
                fit = family_fit("additive"), data = additive_data,
                services = c(family_services, "inference_fix")),
    golden_case("family_additive_equal", "additive",
                "family additive, between_variance equal",
                fit = family_fit("additive", "equal"), data = additive_data,
                services = "inference_opg"),
    golden_case("family_additive_three_classes", "additive",
                "family additive, 3 group classes (a between variance at zero: boundary)",
                fit = family_fit("additive", n_group_classes = 3L),
                data = additive_data),
    golden_case("family_additive_weights", "additive",
                "family additive with sampling weights (robust default)",
                fit = family_fit("additive", weights = "w"),
                data = additive_data),
    golden_case("family_additive_engagement_boundary", "additive",
                "family additive on course_engagement: between variance at zero, Wald refused",
                fit = family_fit("additive", data = engagement, id = "student",
                                 indicators = vars),
                data = engagement, services = "inference_fix"),
    golden_case("family_dispersion_equal", "additive",
                "family dispersion (between_variance equal, the only choice)",
                fit = family_fit("dispersion", "equal"), data = additive_data,
                services = "inference_opg"),
    golden_case("family_dispersion_three_classes", "additive",
                "family dispersion, 3 group classes (singular information)",
                fit = family_fit("dispersion", "equal", n_group_classes = 3L),
                data = additive_data),
    golden_case("family_dispersion_varying_refused", "additive",
                "refusal: family dispersion with between_variance varying",
                fit = family_fit("dispersion", "varying"), data = additive_data),
    golden_case("family_additive_dispersion_varying", "additive",
                "family additive_dispersion, between_variance varying",
                fit = family_fit("additive_dispersion"), data = additive_data,
                services = "inference_opg"),
    golden_case("family_additive_dispersion_equal", "additive",
                "family additive_dispersion, between_variance equal",
                fit = family_fit("additive_dispersion", "equal"),
                data = additive_data),
    golden_case("family_additive_dispersion_three_classes", "additive",
                "family additive_dispersion, varying, 3 group classes (interior)",
                fit = family_fit("additive_dispersion", n_group_classes = 3L),
                data = additive_data),
    golden_case("family_restricted_cross_level_varying", "cross_level",
                "family restricted_cross_level, between varying, profiles varying",
                fit = cross_fit("restricted_cross_level"), data = engagement),
    golden_case("family_restricted_cross_level_equal", "cross_level",
                "family restricted_cross_level, between equal, profiles equal",
                fit = cross_fit("restricted_cross_level", "equal", "equal"),
                data = engagement),
    golden_case("family_full_cross_level_varying", "cross_level",
                "family full_cross_level, between varying, profiles varying",
                fit = cross_fit("full_cross_level"), data = engagement),
    golden_case("family_full_cross_level_equal", "cross_level",
                "family full_cross_level, between equal, profiles equal",
                fit = cross_fit("full_cross_level", "equal", "equal"),
                data = engagement),
    golden_case("family_refusal_categorical", "additive",
                "refusal: additive family with a categorical indicator",
                fit = function() multilpa(engagement, c(vars, "engagement"),
                                          id = "student", n_group_classes = 2,
                                          family = "additive",
                                          categorical = "engagement",
                                          n_starts = 2, seed = 1),
                data = engagement),
    golden_case("family_refusal_fiml", "additive",
                "refusal: additive family with missing = fiml",
                fit = function() multilpa(missing_data, vars, id = "student",
                                          n_group_classes = 2,
                                          family = "additive",
                                          missing = "fiml", n_starts = 2,
                                          seed = 1),
                data = missing_data),
    golden_case("family_refusal_covariates", "additive",
                "refusal: additive family with membership covariates",
                fit = function() multilpa(engagement, vars, id = "student",
                                          n_group_classes = 2,
                                          family = "additive",
                                          group_covariates = "grade_mean",
                                          n_starts = 2, seed = 1),
                data = engagement),
    golden_case("family_refusal_cross_level_categorical", "cross_level",
                "refusal: full_cross_level with a categorical indicator",
                fit = function() multilpa(engagement, c(vars, "engagement"),
                                          id = "student", n_profiles = 2,
                                          n_group_classes = 2,
                                          family = "full_cross_level",
                                          categorical = "engagement",
                                          n_starts = 2, seed = 1),
                data = engagement),
    golden_case("family_enumeration", "additive",
                "enumerate_classes over additive/dispersion/additive_dispersion x 1:2 classes",
                fit = function() enumerate_classes(
                  additive_data, c("y1", "y2"), "group",
                  family = c("additive", "dispersion", "additive_dispersion"),
                  n_group_classes = 1:2, seed = 1, n_starts = 2),
                data = additive_data),
    golden_case("family_bootstrap_lrt", "additive",
                "bootstrap_lrt additive 1 vs 2 group classes (iter 5, seeded)",
                fit = function() {
                  null <- multilpa(additive_data, c("y1", "y2"), id = "group",
                                   n_group_classes = 1, family = "additive",
                                   n_starts = 2, seed = 1)
                  alternative <- multilpa(additive_data, c("y1", "y2"),
                                          id = "group",
                                          n_group_classes = 2,
                                          family = "additive", n_starts = 2,
                                          seed = 1)
                  bootstrap_lrt(null, alternative, iter = 5, n_starts = 2,
                                seed = 1)
                },
                data = additive_data),
    golden_case("family_bootstrap_lrt_refused", "additive",
                "refusal: bootstrap_lrt between additive and dispersion (not nested)",
                fit = function() {
                  additive <- multilpa(additive_data, c("y1", "y2"),
                                       id = "group", n_group_classes = 2,
                                       family = "additive", n_starts = 2,
                                       seed = 1)
                  dispersion <- multilpa(additive_data, c("y1", "y2"),
                                         id = "group",
                                         n_group_classes = 2,
                                         family = "dispersion",
                                         between_variance = "equal",
                                         n_starts = 2, seed = 1)
                  bootstrap_lrt(additive, dispersion, iter = 3, seed = 1)
                },
                data = additive_data)
  )
  c(covariates, families)
}
