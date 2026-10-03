# Section: profiles
#
# lpa(), lca(), multilpa() and multilca() with the profile family: the
# fourteen covariance structures, single- and two-level fits, missing data,
# categorical, ordinal and count indicators, sampling weights, the noise
# component, the conjugate prior, held measurement, centring, ordering,
# acceleration, start selection and user starts, plus refusals.

#' First 150 rows of `srl`, three indicators (single-level structure fits)
golden_profiles_srl <- function() {
  data <- srl[seq_len(150L), c("cognitive_strategies", "intrinsic_value",
                               "self_efficacy")]
  row.names(data) <- NULL
  data
}

golden_profiles_srl_vars <- c("cognitive_strategies", "intrinsic_value",
                              "self_efficacy")

#' `golden_engagement()` with missing indicator values at fixed positions
golden_profiles_engagement_missing <- function() {
  data <- golden_engagement()
  n <- nrow(data)
  data$browse[seq(5L, n, by = 17L)] <- NA
  data$lectures[seq(9L, n, by = 23L)] <- NA
  data$forum_read[c(3L, 9L)] <- NA
  data
}

#' `golden_engagement()` with a sampling weight constant within student
golden_profiles_engagement_weighted <- function() {
  data <- golden_engagement()
  data$w <- 1 + (data$student %% 3L) / 2
  data
}

#' `golden_profiles_srl()` with a row-level sampling weight
golden_profiles_srl_weighted <- function() {
  data <- golden_profiles_srl()
  data$w <- 1 + (seq_len(nrow(data)) %% 4L) / 4
  data
}

#' The first 30 students of `student_esm`
golden_profiles_esm <- function() golden_esm(30L)

golden_profiles_esm_binary <- c("time_with_friends", "on_social_media",
                                "tv_video_games", "sports")

#' The fourteen mclust structures, each fitted single-level on the srl subset
golden_profiles_structure_cases <- function() {
  data <- golden_profiles_srl()
  vars <- golden_profiles_srl_vars
  codes <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
             "VVE", "EEV", "VEV", "EVV", "VVV")
  lapply(codes, function(code) {
    force(code)
    golden_case(
      sprintf("profiles_structure_%s", code), "multilpa",
      sprintf("lpa() single-level, 2 profiles, model = \"%s\" (srl, 150 rows)",
              code),
      fit = function() lpa(data, vars, n_profiles = 2, model = code,
                           n_starts = 2, seed = 1),
      data = data)
  })
}

golden_cases_profiles <- function() {
  engagement <- golden_engagement()
  engagement40 <- golden_engagement(40L)
  vars <- golden_engagement_vars
  srl_data <- golden_profiles_srl()
  srl_vars <- golden_profiles_srl_vars
  missing_data <- golden_profiles_engagement_missing()
  weighted <- golden_profiles_engagement_weighted()
  srl_weighted <- golden_profiles_srl_weighted()
  esm <- golden_profiles_esm()
  binary <- golden_profiles_esm_binary
  hours <- golden_hours(60L)
  stage_one <- function() {
    multilpa(engagement, vars, id = "student", n_profiles = 2,
             n_group_classes = 1, n_starts = 2, seed = 1)
  }

  basic <- list(
    golden_case(
      "profiles_two_level_basic", "multilpa",
      "multilpa two-level, 2 profiles x 2 group classes, diagonal varying",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, n_starts = 2, seed = 1),
      data = engagement,
      services = c("diagnostics", "three_step_bch", "bootstrap_lrt"),
      outcome = "previous_grade",
      null_fit = function() multilpa(engagement, vars, id = "student",
                                     n_profiles = 2, n_group_classes = 1,
                                     n_starts = 2, seed = 1)),
    golden_case(
      "profiles_refusal_noise_two_level", "multilpa",
      "noise = TRUE with two group classes is refused",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, noise = TRUE, n_starts = 2,
                                seed = 1),
      data = engagement),
    golden_case(
      "profiles_lpa_single_varying", "multilpa",
      "lpa() single-level srl (all 5 indicators), 3 profiles, varying diagonal",
      fit = function() lpa(srl, names(srl), n_profiles = 3, n_starts = 2,
                           seed = 1),
      data = srl,
      services = c("descriptives", "predict", "starting_values",
                   "inference_opg", "report")),
    golden_case(
      "profiles_lpa_single_equal", "multilpa",
      "lpa() single-level, variance_model = \"equal\"",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                           variance_model = "equal", n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_lpa_single_full", "multilpa",
      "lpa() single-level, covariance_model = \"full\", varying",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                           covariance_model = "full", n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_lpa_single_full_equal", "multilpa",
      "lpa() single-level, full covariance shared across profiles",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                           covariance_model = "full", variance_model = "equal",
                           n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_lpa_engagement_id_null", "multilpa",
      "multilpa(id = NULL) on engagement rows: single-level message, 2 profiles",
      fit = function() multilpa(engagement, vars, id = NULL, n_profiles = 2,
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_lpa_one_profile", "multilpa",
      "lpa() with a single profile",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 1, n_starts = 1,
                           seed = 1),
      data = srl_data),
    golden_case(
      "profiles_two_level_equal", "multilpa",
      "multilpa two-level, variance_model = \"equal\"",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, variance_model = "equal",
                                n_starts = 2, seed = 1),
      data = engagement,
      services = c("predict", "starting_values")),
    golden_case(
      "profiles_two_level_full", "multilpa",
      "multilpa two-level, covariance_model = \"full\" (40 students, so robust/OPG apply)",
      fit = function() multilpa(engagement40, vars, id = "student",
                                n_profiles = 2, n_group_classes = 2,
                                covariance_model = "full", n_starts = 2,
                                seed = 1),
      data = engagement40,
      services = "inference_opg"),
    golden_case(
      "profiles_two_level_three_profiles", "multilpa",
      "multilpa two-level, 3 profiles x 2 group classes",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 3,
                                n_group_classes = 2, n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_two_level_one_class", "multilpa",
      "multilpa two-level data with one group class",
      fit = stage_one, data = engagement,
      services = c("starting_values_measurement", "diagnostics"))
  )

  structures <- c(
    golden_profiles_structure_cases(),
    list(
      golden_case(
        "profiles_structure_pieces_VEV", "multilpa",
        "structure by volume/shape/orientation (varying/equal/varying = VEV)",
        fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                             volume = "varying", shape = "equal",
                             orientation = "varying", n_starts = 2, seed = 1),
        data = srl_data),
      golden_case(
        "profiles_structure_pieces_EII", "multilpa",
        "structure by pieces: equal volume, spherical shape (EII)",
        fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                             volume = "equal", shape = "spherical",
                             n_starts = 2, seed = 1),
        data = srl_data),
      golden_case(
        "profiles_structure_three_profiles_EVE", "multilpa",
        "EVE (minorize-maximize orientation step) with 3 profiles",
        fit = function() lpa(srl_data, srl_vars, n_profiles = 3, model = "EVE",
                             n_starts = 2, seed = 1),
        data = srl_data),
      golden_case(
        "profiles_two_level_VEI", "multilpa",
        "two-level multilpa with model = \"VEI\"",
        fit = function() multilpa(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  model = "VEI", n_starts = 2, seed = 1),
        data = engagement,
        services = "inference_bootstrap"),
      golden_case(
        "profiles_two_level_EEV", "multilpa",
        "two-level multilpa with model = \"EEV\"",
        fit = function() multilpa(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  model = "EEV", n_starts = 2, seed = 1),
        data = engagement),
      golden_case(
        "profiles_two_level_VVE", "multilpa",
        "two-level multilpa with model = \"VVE\"",
        fit = function() multilpa(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  model = "VVE", n_starts = 2, seed = 1),
        data = engagement),
      golden_case(
        "profiles_refusal_model_and_pieces", "multilpa",
        "model = together with volume/shape/orientation is refused",
        fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "VVV",
                             volume = "equal", n_starts = 2, seed = 1),
        data = srl_data)
    )
  )

  missing_cases <- list(
    golden_case(
      "profiles_fiml_two_level", "multilpa",
      "missing = \"fiml\", two-level diagonal",
      fit = function() multilpa(missing_data, vars, id = "student",
                                n_profiles = 2, n_group_classes = 2,
                                missing = "fiml", n_starts = 2, seed = 1),
      data = missing_data,
      services = "descriptives"),
    golden_case(
      "profiles_fiml_single_full", "multilpa",
      "missing = \"fiml\", single-level full covariance",
      fit = function() multilpa(missing_data, vars, id = NULL, n_profiles = 2,
                                covariance_model = "full", missing = "fiml",
                                n_starts = 2, seed = 1),
      data = missing_data),
    golden_case(
      "profiles_refusal_missing_error", "multilpa",
      "missing values under missing = \"error\" are refused",
      fit = function() multilpa(missing_data, vars, id = "student",
                                n_profiles = 2, n_group_classes = 2,
                                n_starts = 2, seed = 1),
      data = missing_data)
  )

  indicator_cases <- list(
    golden_case(
      "profiles_lca_single", "multilca",
      "lca() single-level on 4 binary student_esm items, 2 classes",
      fit = function() lca(esm, binary, n_classes = 2, n_starts = 2, seed = 1),
      data = esm,
      services = c("predict", "starting_values", "inference_bootstrap")),
    golden_case(
      "profiles_multilca_two_level", "multilca",
      "multilca() two-level on binary items, 2 classes x 2 group classes",
      fit = function() multilca(esm, binary, id = "student", n_profiles = 2,
                                n_group_classes = 2, n_starts = 2, seed = 1),
      data = esm,
      services = "diagnostics"),
    golden_case(
      "profiles_mixed_categorical", "multilpa",
      "multilpa mixed: 2 continuous (happy, relaxed) + 2 categorical items",
      fit = function() multilpa(esm, c("happy", "relaxed", "sports", "reading"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 2,
                                categorical = c("sports", "reading"),
                                n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "profiles_ordinal_single", "multilpa",
      "lpa() with ordinal items (worried, exhausted) + continuous happy",
      fit = function() lpa(esm, c("happy", "worried", "exhausted"),
                           n_profiles = 2, ordinal = c("worried", "exhausted"),
                           n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "profiles_ordinal_two_level", "multilpa",
      "multilpa two-level with an ordinal item and a categorical item",
      fit = function() multilpa(esm, c("happy", "worried", "sports"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 2, ordinal = "worried",
                                categorical = "sports", n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "profiles_count_poisson", "multilpa",
      "lpa() with a Poisson count (questions) + continuous score",
      fit = function() lpa(hours, c("score", "questions"), n_profiles = 2,
                           count = "questions", n_starts = 2, seed = 1),
      data = hours,
      services = "predict"),
    golden_case(
      "profiles_count_nb_varying", "multilpa",
      "negative binomial count, count_dispersion = \"varying\"",
      fit = function() lpa(hours, c("score", "questions"), n_profiles = 2,
                           count = "questions",
                           count_model = "negative_binomial",
                           count_dispersion = "varying", n_starts = 2, seed = 1),
      data = hours),
    golden_case(
      "profiles_count_nb_equal", "multilpa",
      "negative binomial count, count_dispersion = \"equal\"",
      fit = function() lpa(hours, c("score", "questions"), n_profiles = 2,
                           count = "questions",
                           count_model = "negative_binomial",
                           count_dispersion = "equal", n_starts = 2, seed = 1),
      data = hours),
    golden_case(
      "profiles_count_two_level", "multilpa",
      "multilpa two-level: Poisson count + ordinal + continuous (mixed)",
      fit = function() multilpa(esm, c("happy", "worried", "exhausted"),
                                id = "student", n_profiles = 2,
                                n_group_classes = 2, ordinal = "worried",
                                count = "exhausted", n_starts = 2, seed = 1),
      data = esm,
      services = "inference_bootstrap"),
    golden_case(
      "profiles_mixed_all_types", "multilpa",
      "single-level: continuous + categorical + ordinal + count",
      fit = function() lpa(esm, c("happy", "sports", "worried", "exhausted"),
                           n_profiles = 2, categorical = "sports",
                           ordinal = "worried", count = "exhausted",
                           n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "profiles_extra_fiml", "multilpa",
      "ordinal + count indicators with missing = \"fiml\"",
      fit = function() {
        data <- esm
        data$worried[seq(4L, nrow(data), by = 19L)] <- NA
        data$exhausted[seq(7L, nrow(data), by = 29L)] <- NA
        lpa(data, c("happy", "worried", "exhausted"), n_profiles = 2,
            ordinal = "worried", count = "exhausted", missing = "fiml",
            n_starts = 2, seed = 1)
      },
      data = esm),
    golden_case(
      "profiles_refusal_ordinal_fixed", "multilpa",
      "ordinal indicators with fixed = measurement are refused",
      fit = function() {
        first <- lpa(esm, c("happy", "worried"), n_profiles = 2,
                     ordinal = "worried", n_starts = 2, seed = 1)
        lpa(esm, c("happy", "worried"), n_profiles = 2, ordinal = "worried",
            start = starting_values(first, what = "measurement"),
            fixed = "measurement", n_starts = 2, seed = 1)
      },
      data = esm)
  )

  weight_cases <- list(
    golden_case(
      "profiles_weights_single", "multilpa",
      "lpa() with row sampling weights",
      fit = function() lpa(srl_weighted, srl_vars, n_profiles = 2,
                           weights = "w", n_starts = 2, seed = 1),
      data = srl_weighted),
    golden_case(
      "profiles_weights_two_level", "multilpa",
      "multilpa two-level with group sampling weights",
      fit = function() multilpa(weighted, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, weights = "w",
                                n_starts = 2, seed = 1),
      data = weighted),
    golden_case(
      "profiles_weights_categorical", "multilpa",
      "weighted two-level fit with categorical indicators",
      fit = function() {
        data <- esm
        data$w <- 1 + (data$student %% 2L)
        multilca(data, binary, id = "student", n_profiles = 2,
                 n_group_classes = 2, weights = "w", n_starts = 2, seed = 1)
      },
      data = esm),
    golden_case(
      "profiles_refusal_weights_prior", "multilpa",
      "weights together with prior are refused",
      fit = function() lpa(srl_weighted, srl_vars, n_profiles = 2,
                           weights = "w", prior = prior_control(),
                           n_starts = 2, seed = 1),
      data = srl_weighted)
  )

  noise_prior_cases <- list(
    golden_case(
      "profiles_noise_single", "multilpa",
      "noise = TRUE single-level, VVI",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, noise = TRUE,
                           n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_noise_VVV", "multilpa",
      "noise = TRUE single-level, model = \"VVV\"",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "VVV",
                           noise = TRUE, n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_prior_VVI", "multilpa",
      "prior = prior_control(), VVI",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2,
                           prior = prior_control(), n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_prior_EEE", "multilpa",
      "prior = prior_control(), EEE",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "EEE",
                           prior = prior_control(), n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_prior_VEV_shrinkage", "multilpa",
      "prior_control(shrinkage = 0.1), VEV",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "VEV",
                           prior = prior_control(shrinkage = 0.1),
                           n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_prior_noise_EII", "multilpa",
      "prior and noise together, EII",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "EII",
                           prior = prior_control(), noise = TRUE,
                           n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_refusal_prior_VEE", "multilpa",
      "prior with VEE (no mclust prior) is refused",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "VEE",
                           prior = prior_control(), n_starts = 2, seed = 1),
      data = srl_data)
  )

  fixed_cases <- list(
    golden_case(
      "profiles_fixed_measurement", "multilpa",
      "fixed = \"measurement\" from a one-class stage, 2 group classes",
      fit = function() multilpa(
        engagement, vars, id = "student", n_profiles = 2, n_group_classes = 2,
        start = starting_values(stage_one(), what = "measurement"),
        fixed = "measurement", n_starts = 2, seed = 1),
      data = engagement,
      services = c("inference_fix", "bootstrap_lrt"),
      null_fit = function() multilpa(
        engagement, vars, id = "student", n_profiles = 2, n_group_classes = 1,
        start = starting_values(stage_one(), what = "measurement"),
        fixed = "measurement", n_starts = 2, seed = 1)),
    golden_case(
      "profiles_fixed_means", "multilpa",
      "fixed = \"means\" only",
      fit = function() multilpa(
        engagement, vars, id = "student", n_profiles = 2, n_group_classes = 2,
        start = starting_values(stage_one(), what = "measurement"),
        fixed = "means", n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_fixed_categorical", "multilpa",
      "fixed = \"response_probabilities\" on a latent class model",
      fit = function() {
        first <- multilca(esm, binary, id = "student", n_profiles = 2,
                          n_group_classes = 1, n_starts = 2, seed = 1)
        multilca(esm, binary, id = "student", n_profiles = 2,
                 n_group_classes = 2,
                 start = starting_values(first, what = "measurement"),
                 fixed = "response_probabilities", n_starts = 2, seed = 1)
      },
      data = esm),
    golden_case(
      "profiles_evaluate_only", "multilpa",
      "max_iter = 0 with start: evaluate a supplied parameter set",
      fit = function() multilpa(
        engagement, vars, id = "student", n_profiles = 2, n_group_classes = 1,
        start = starting_values(stage_one()), max_iter = 0, n_starts = 2,
        seed = 1),
      data = engagement),
    golden_case(
      "profiles_user_start", "multilpa",
      "user start list (means, variances, probabilities), two-level",
      fit = function() multilpa(
        engagement, vars, id = "student", n_profiles = 2, n_group_classes = 2,
        start = list(
          means = rbind(c(0.5, 0.5, 0.5), c(-0.8, -0.6, -0.8)),
          variances = matrix(0.5, 2, 3),
          profile_probabilities = rbind(c(0.8, 0.2), c(0.2, 0.8)),
          group_probabilities = c(0.6, 0.4)),
        n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_refusal_fixed_without_start", "multilpa",
      "fixed = \"means\" without start is refused",
      fit = function() multilpa(engagement, vars, id = "student",
                                n_profiles = 2, n_group_classes = 2,
                                fixed = "means", n_starts = 2, seed = 1),
      data = engagement)
  )

  option_cases <- list(
    golden_case(
      "profiles_centering_person", "multilpa",
      "centering = \"person\", two-level",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, centering = "person",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_centering_grand", "multilpa",
      "centering = \"grand\", two-level",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, centering = "grand",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_time_sequences", "multilpa",
      "time = \"sequence\": sequence tables, get_tna / get_group_tna",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, time = "sequence",
                                n_starts = 2, seed = 1),
      data = engagement,
      services = c("get_tna", "get_group_tna")),
    golden_case(
      "profiles_acceleration_none", "multilpa",
      "acceleration = \"none\" (plain EM), two-level",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 2,
                                n_group_classes = 2, acceleration = "none",
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_acceleration_squarem_full", "multilpa",
      "acceleration = \"squarem\" explicit, single-level full covariance",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 3,
                           covariance_model = "full",
                           acceleration = "squarem", n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_acceleration_none_full", "multilpa",
      "acceleration = \"none\", single-level full covariance",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 3,
                           covariance_model = "full", acceleration = "none",
                           n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_select_converged", "multilpa",
      "select_start = \"converged\": max_iter = 60 leaves the best start unconverged, a worse one converged",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 3,
                                n_group_classes = 2, select_start = "converged",
                                max_iter = 60, n_starts = 3, seed = 2),
      data = engagement),
    golden_case(
      "profiles_select_likelihood", "multilpa",
      "select_start = \"likelihood\" on the same short run (picks the unconverged best)",
      fit = function() multilpa(engagement, vars, id = "student", n_profiles = 3,
                                n_group_classes = 2, select_start = "likelihood",
                                max_iter = 60, n_starts = 3, seed = 2),
      data = engagement),
    golden_case(
      "profiles_min_variance_tol", "multilpa",
      "non-default tol and min_variance",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, tol = 1e-6,
                           min_variance = 1e-3, n_starts = 2, seed = 1),
      data = srl_data),
    golden_case(
      "profiles_refusal_two_classes_no_id", "multilpa",
      "n_group_classes > 1 with id = NULL is refused",
      fit = function() multilpa(engagement, vars, id = NULL, n_profiles = 2,
                                n_group_classes = 2, n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_refusal_missing_id", "multilpa",
      "omitting id is refused",
      fit = function() multilpa(engagement, vars, n_profiles = 2,
                                n_starts = 2, seed = 1),
      data = engagement),
    golden_case(
      "profiles_refusal_bad_model", "multilpa",
      "an unknown model code is refused",
      fit = function() lpa(srl_data, srl_vars, n_profiles = 2, model = "XYZ",
                           n_starts = 2, seed = 1),
      data = srl_data)
  )

  c(basic, structures, missing_cases, indicator_cases, weight_cases,
    noise_prior_cases, fixed_cases, option_cases)
}
