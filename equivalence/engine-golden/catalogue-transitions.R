# Section: transitions. lta() through both of its engines -- the homogeneous
# engine (`multilpa_transitions`) and the general engine (`multilpa_lta`) --
# plus the transition bootstrap, the bootstrap likelihood-ratio test between
# nested transition fits, the transition enumeration and the tna hand-over.

#' A short engagement panel: the first `n` students, first `occasions` courses
#' @param n Number of students.
#' @param occasions Number of leading courses (`sequence`) kept per student.
#' @return A data.frame with a per-student sampling `weight` column added.
golden_lta_panel <- function(n = 40L, occasions = 4L) {
  data <- subset(course_engagement, student <= n & sequence <= occasions)
  row.names(data) <- NULL
  data$weight <- 1 + (data$student %% 3L) / 2
  data
}

#' The engagement panel with holes: some students skip a course position, so
#' "observed" and "grid" occasions differ
#' @return A data.frame.
golden_lta_gapped <- function() {
  data <- golden_lta_panel(40L, 5L)
  drop <- data$sequence == 3L & data$student %% 4L == 0L
  data <- data[!drop, , drop = FALSE]
  row.names(data) <- NULL
  data
}

#' The engagement panel with missing indicator values (for missing = "fiml")
#' @return A data.frame.
golden_lta_missing <- function() {
  data <- golden_lta_panel(40L, 4L)
  data$browse[seq(3L, nrow(data), by = 11L)] <- NA_real_
  data$forum_read[seq(5L, nrow(data), by = 13L)] <- NA_real_
  data
}

#' Weekly study panel: `week` occasions, binary, factor, count indicators
#' @param n Number of students.
#' @return A data.frame with a per-student `weight` and an ordered `effort`.
golden_lta_hours <- function(n = 60L) {
  data <- subset(study_hours, student <= n)
  row.names(data) <- NULL
  data$weight <- 1 + (data$student %% 4L) / 3
  data$effort <- cut(data$hours, c(-Inf, 4, 9, Inf),
                     labels = c("low", "mid", "high"), ordered_result = TRUE)
  data
}

#' Experience-sampling panel with a single within-student occasion index
#' @param n Number of students.
#' @return A data.frame with `occasion = day * 5 + beep`.
golden_lta_esm <- function(n = 20L) {
  data <- subset(student_esm, student <= n)
  row.names(data) <- NULL
  data$occasion <- data$day * 5L + data$beep
  data
}

golden_cases_transitions <- function() {
  vars <- golden_engagement_vars
  panel <- golden_lta_panel()
  gapped <- golden_lta_gapped()
  missing_panel <- golden_lta_missing()
  hours <- golden_lta_hours()
  esm <- golden_lta_esm()
  engagement <- golden_engagement()
  panel_lta <- function(...) {
    lta(panel, vars, id = "student", n_profiles = 2, time = "sequence",
        n_starts = 2, seed = 1, ...)
  }
  homogeneous_services <- c("diagnostics", "descriptives", "inference_opg",
                            "get_tna", "get_group_tna", "inference_bootstrap")
  # The transition bootstrap needs `data` (each resample keeps a person's
  # occasions and covariates); the `inference_bootstrap` service omits it and
  # so records that refusal, and these cases record the bootstrap itself.
  bootstrap <- function(fit, data) {
    parameter_inference(fit, data = data, method = "bootstrap", iter = 5L,
                        n_starts = 2L, seed = 1L)
  }

  homogeneous <- list(
    golden_case(
      "lta_homogeneous_basic", "lta (homogeneous)",
      "lta() 2 profiles, 1 group class, diagonal varying, 20 students x ~14 courses",
      fit = function() {
        lta(engagement, vars, id = "student", n_profiles = 2, time = "sequence",
            n_starts = 2, seed = 1)
      },
      data = engagement, services = homogeneous_services),
    golden_case(
      "lta_homogeneous_group_classes", "lta (homogeneous)",
      "lta() 2 profiles x 2 group classes (mixture of transition matrices)",
      fit = function() panel_lta(n_group_classes = 2),
      data = panel, services = c("get_tna", "get_group_tna", "diagnostics")),
    golden_case(
      "lta_homogeneous_equal_variance", "lta (homogeneous)",
      "lta() variance_model = equal",
      fit = function() panel_lta(variance_model = "equal"),
      data = panel),
    golden_case(
      "lta_homogeneous_full_covariance", "lta (homogeneous)",
      "lta() covariance_model = full (varying)",
      fit = function() panel_lta(covariance_model = "full"),
      data = panel, services = "inference_opg"),
    golden_case(
      "lta_homogeneous_model_EEI", "lta (homogeneous)",
      "lta() model = EEI",
      fit = function() panel_lta(model = "EEI"),
      data = panel),
    golden_case(
      "lta_homogeneous_model_VVI", "lta (homogeneous)",
      "lta() model = VVI",
      fit = function() panel_lta(model = "VVI"),
      data = panel),
    golden_case(
      "lta_homogeneous_model_EEE", "lta (homogeneous)",
      "lta() model = EEE",
      fit = function() panel_lta(model = "EEE"),
      data = panel),
    golden_case(
      "lta_homogeneous_model_VVV", "lta (homogeneous)",
      "lta() model = VVV",
      fit = function() panel_lta(model = "VVV"),
      data = panel),
    golden_case(
      "lta_homogeneous_model_VEI", "lta (homogeneous)",
      "lta() model = VEI (inference refused for this structure)",
      fit = function() panel_lta(model = "VEI"),
      data = panel),
    golden_case(
      "lta_homogeneous_categorical", "lta (homogeneous)",
      "lta() categorical binary indicators (ESM activities)",
      fit = function() {
        lta(esm, c("time_with_friends", "on_social_media", "tv_video_games"),
            id = "student", n_profiles = 2, time = "occasion",
            categorical = c("time_with_friends", "on_social_media",
                            "tv_video_games"),
            n_starts = 2, seed = 1)
      },
      data = esm, services = "inference_fix"),
    golden_case(
      "lta_homogeneous_mixed", "lta (homogeneous)",
      "lta() mixed indicators: two Gaussian, one categorical factor",
      fit = function() {
        lta(hours, c("hours", "sleep", "strategy"), id = "student",
            n_profiles = 2, time = "week", categorical = "strategy",
            n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_homogeneous_fiml", "lta (homogeneous)",
      "lta() missing = fiml with missing indicator values",
      fit = function() {
        lta(missing_panel, vars, id = "student", n_profiles = 2,
            time = "sequence", missing = "fiml", n_starts = 2, seed = 1)
      },
      data = missing_panel),
    golden_case(
      "lta_homogeneous_occasions_observed", "lta (homogeneous)",
      "lta() occasions = observed on a panel with skipped positions",
      fit = function() {
        lta(gapped, vars, id = "student", n_profiles = 2, time = "sequence",
            occasions = "observed", n_starts = 2, seed = 1)
      },
      data = gapped),
    golden_case(
      "lta_homogeneous_occasions_grid", "lta (homogeneous)",
      "lta() occasions = grid on a panel with skipped positions",
      fit = function() {
        lta(gapped, vars, id = "student", n_profiles = 2, time = "sequence",
            occasions = "grid", n_starts = 2, seed = 1)
      },
      data = gapped),
    golden_case(
      "lta_homogeneous_select_converged", "lta (homogeneous)",
      "lta() select_start = converged, 3 profiles, 2 group classes",
      fit = function() {
        lta(panel, vars, id = "student", n_profiles = 3, n_group_classes = 2,
            time = "sequence", select_start = "converged", n_starts = 3,
            seed = 1)
      },
      data = panel, services = "inference_fix"),
    golden_case(
      "lta_homogeneous_bootstrap", "lta (homogeneous)",
      "parameter_inference(method = bootstrap, iter = 5) on a homogeneous fit",
      fit = function() bootstrap(panel_lta(), panel),
      data = panel)
  )

  general <- list(
    golden_case(
      "lta_general_occasion_transitions", "lta (general)",
      "lta() transitions = occasion (one matrix per move)",
      fit = function() panel_lta(transitions = "occasion"),
      data = panel, services = "inference_opg"),
    golden_case(
      "lta_general_order2", "lta (general)",
      "lta() order = 2 (second-order Markov chain)",
      fit = function() panel_lta(order = 2),
      data = panel, services = "inference_fix"),
    golden_case(
      "lta_general_mover_stayer", "lta (general)",
      "lta() mover_stayer = TRUE",
      fit = function() panel_lta(mover_stayer = TRUE),
      data = panel),
    golden_case(
      "lta_general_transition_covariates", "lta (general)",
      "lta() transition_covariates = previous_grade",
      fit = function() panel_lta(transition_covariates = "previous_grade"),
      data = panel, services = c("inference_opg", "get_tna")),
    golden_case(
      "lta_general_initial_covariates", "lta (general)",
      "lta() initial_covariates = previous_grade",
      fit = function() panel_lta(initial_covariates = "previous_grade"),
      data = panel),
    golden_case(
      "lta_general_occasion_measurement", "lta (general)",
      "lta() measurement = occasion (non-invariant profiles)",
      fit = function() panel_lta(measurement = "occasion"),
      data = panel),
    golden_case(
      "lta_general_group_classes", "lta (general)",
      "lta() n_group_classes = 2 with transition covariates",
      fit = function() {
        panel_lta(n_group_classes = 2, transition_covariates = "previous_grade")
      },
      data = panel),
    golden_case(
      "lta_general_full_covariance", "lta (general)",
      "lta() covariance_model = full with occasion transitions",
      fit = function() {
        panel_lta(covariance_model = "full", transitions = "occasion")
      },
      data = panel),
    golden_case(
      "lta_general_model_VEI", "lta (general)",
      "lta() model = VEI with transition covariates (EM, inference refused)",
      fit = function() {
        panel_lta(model = "VEI", transition_covariates = "previous_grade")
      },
      data = panel),
    golden_case(
      "lta_general_categorical", "lta (general)",
      "lta() categorical indicators with occasion transitions",
      fit = function() {
        lta(hours, c("passed", "strategy"), id = "student", n_profiles = 2,
            time = "week", categorical = c("passed", "strategy"),
            transitions = "occasion", n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_ordinal", "lta (general)",
      "lta() ordinal indicator (adjacent-category logit) with a Gaussian one",
      fit = function() {
        lta(hours, c("sleep", "effort"), id = "student", n_profiles = 2,
            time = "week", ordinal = "effort", n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_count_poisson", "lta (general)",
      "lta() Poisson count indicator with a Gaussian one",
      fit = function() {
        lta(hours, c("sleep", "questions"), id = "student", n_profiles = 2,
            time = "week", count = "questions", n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_count_nb_varying", "lta (general)",
      "lta() negative binomial count, count_dispersion = varying",
      fit = function() {
        lta(hours, c("sleep", "questions"), id = "student", n_profiles = 2,
            time = "week", count = "questions",
            count_model = "negative_binomial", count_dispersion = "varying",
            n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_count_nb_equal", "lta (general)",
      "lta() negative binomial count, count_dispersion = equal",
      fit = function() {
        lta(hours, c("sleep", "questions"), id = "student", n_profiles = 2,
            time = "week", count = "questions",
            count_model = "negative_binomial", count_dispersion = "equal",
            n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_weights", "lta (general)",
      "lta() weights (pseudo ML; robust default, observed refused)",
      fit = function() panel_lta(weights = "weight"),
      data = panel),
    golden_case(
      "lta_general_fiml_covariates", "lta (general)",
      "lta() missing = fiml with transition covariates",
      fit = function() {
        lta(missing_panel, vars, id = "student", n_profiles = 2,
            time = "sequence", missing = "fiml",
            transition_covariates = "previous_grade", n_starts = 2, seed = 1)
      },
      data = missing_panel),
    golden_case(
      "lta_general_mixed_extensions", "lta (general)",
      paste("lta() combination: mixed categorical/ordinal/count indicators,",
            "weights, fiml, occasion transitions"),
      fit = function() {
        data <- hours
        data$sleep[seq(4L, nrow(data), by = 17L)] <- NA_real_
        lta(data, c("sleep", "strategy", "effort", "questions"),
            id = "student", n_profiles = 2, time = "week",
            categorical = "strategy", ordinal = "effort", count = "questions",
            missing = "fiml", weights = "weight", transitions = "occasion",
            n_starts = 2, seed = 1)
      },
      data = hours),
    golden_case(
      "lta_general_order2_mover_stayer", "lta (general)",
      "lta() combination: order 2, mover-stayer, occasion transitions",
      fit = function() {
        panel_lta(order = 2, mover_stayer = TRUE, transitions = "occasion")
      },
      data = panel),
    golden_case(
      "lta_general_bootstrap", "lta (general)",
      "parameter_inference(method = bootstrap, iter = 5) on a covariate transition fit",
      fit = function() {
        bootstrap(panel_lta(transition_covariates = "previous_grade"), panel)
      },
      data = panel),
    golden_case(
      "lta_general_weights_bootstrap", "lta (general)",
      "parameter_inference(method = bootstrap) on a weighted transition fit",
      fit = function() bootstrap(panel_lta(weights = "weight"), panel),
      data = panel)
  )

  refusals <- list(
    golden_case(
      "lta_refusal_one_profile", "lta (refusal)",
      "lta() n_profiles = 1 is refused (latents_bad_transition)",
      fit = function() {
        lta(panel, vars, id = "student", n_profiles = 1, time = "sequence",
            n_starts = 2, seed = 1)
      },
      data = panel),
    golden_case(
      "lta_refusal_order3", "lta (refusal)",
      "lta() order = 3 is refused",
      fit = function() panel_lta(order = 3),
      data = panel),
    golden_case(
      "lta_refusal_no_time", "lta (refusal)",
      "lta() without time is refused",
      fit = function() {
        lta(panel, vars, id = "student", n_profiles = 2, time = NULL,
            n_starts = 2, seed = 1)
      },
      data = panel),
    golden_case(
      "lta_refusal_missing_error", "lta (refusal)",
      "lta() missing indicators under missing = error are refused",
      fit = function() {
        lta(missing_panel, vars, id = "student", n_profiles = 2,
            time = "sequence", n_starts = 2, seed = 1)
      },
      data = missing_panel),
    golden_case(
      "lta_refusal_bad_model", "lta (refusal)",
      "lta() an unknown covariance structure code is refused",
      fit = function() panel_lta(model = "XYZ"),
      data = panel)
  )

  comparison <- list(
    golden_case(
      "lta_bootstrap_lrt_covariates", "lta (compare)",
      "bootstrap_lrt() homogeneous null vs transition-covariate alternative",
      fit = function() {
        null <- panel_lta()
        alternative <- panel_lta(transition_covariates = "previous_grade")
        bootstrap_lrt(null, alternative, data = panel, iter = 5, n_starts = 2,
                      seed = 1)
      },
      data = panel),
    golden_case(
      "lta_bootstrap_lrt_group_classes", "lta (compare)",
      "bootstrap_lrt() 1 vs 2 group classes of homogeneous transition fits",
      fit = function() {
        null <- panel_lta()
        alternative <- panel_lta(n_group_classes = 2)
        bootstrap_lrt(null, alternative, data = panel, iter = 5, n_starts = 1,
                      seed = 1)
      },
      data = panel),
    golden_case(
      "lta_refusal_bootstrap_lrt_reversed", "lta (compare)",
      "bootstrap_lrt() with the models in the wrong order is refused",
      fit = function() {
        null <- panel_lta()
        alternative <- panel_lta(transition_covariates = "previous_grade")
        bootstrap_lrt(alternative, null, data = panel, iter = 5, seed = 1)
      },
      data = panel),
    golden_case(
      "lta_enumeration", "lta (compare)",
      "enumerate_classes(time =) over 2:3 profiles x 1:2 group classes",
      fit = function() {
        enumerate_classes(panel, vars, "student", time = "sequence",
                          n_profiles = 2:3, n_group_classes = 1:2,
                          n_starts = 2, seed = 1)
      },
      data = panel),
    golden_case(
      "lta_enumeration_structures", "lta (compare)",
      "enumerate_classes(time =) crossing structures VVI and VEI",
      fit = function() {
        enumerate_classes(panel, vars, "student", time = "sequence",
                          n_profiles = 2, model = c("VVI", "VEI"),
                          n_starts = 2, seed = 1)
      },
      data = panel),
    golden_case(
      "lta_enumeration_candidate", "lta (compare)",
      "candidate_fit() picks the 3-profile, 1-class fit out of a transition grid",
      fit = function() {
        grid <- enumerate_classes(panel, vars, "student", time = "sequence",
                                  n_profiles = 2:3, n_group_classes = 1,
                                  n_starts = 2, seed = 1)
        candidate_fit(grid, n_profiles = 3, n_group_classes = 1)
      },
      data = panel),
    golden_case(
      "lta_refusal_enumeration_one_profile", "lta (compare)",
      "enumerate_classes(time =) with a one-profile candidate is refused",
      fit = function() {
        enumerate_classes(panel, vars, "student", time = "sequence",
                          n_profiles = 1:2, n_starts = 2, seed = 1)
      },
      data = panel)
  )

  c(homogeneous, general, refusals, comparison)
}
