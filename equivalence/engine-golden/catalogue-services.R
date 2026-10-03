# Section: services -- the verbs that sit on top of a fit: enumeration,
# bootstrap likelihood-ratio tests, Wald / bootstrap / boundary inference,
# three-step and R3STEP, sensitivity, multiple-imputation pooling, staged fits,
# starting values, diagnostics, descriptives, prediction and report.
#
# A case whose `fit()` returns the service's own result (a table, an
# enumeration, a pooled object) is fingerprinted like any other object.

#' Evaluate `expr` under `seed`, restoring the caller's random state
#' @param seed Seed.
#' @param expr Expression, evaluated lazily.
#' @return The value of `expr`.
golden_svc_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE, after = FALSE)
  set.seed(seed)
  expr
}

#' Engagement data with a group-level outcome/covariate (the student mean of
#' `previous_grade`), constant within each student
golden_svc_engagement <- function() {
  data <- golden_engagement()
  data$grade_mean <- stats::ave(data$previous_grade, data$student)
  data
}

#' The first 150 `srl` respondents
golden_svc_srl <- function() {
  data <- srl[seq_len(150L), ]
  row.names(data) <- NULL
  data
}

#' Binary-item two-level design with a group covariate `z` (tests'
#' `.lca_design()`), plus an item `y6` whose second category never occurs in
#' class-1-defining rows, driving one response probability to its floor
golden_svc_lca_boundary <- function() {
  golden_svc_with_seed(2024L, {
    n_groups <- 60L
    per <- 8L
    g <- rep(seq_len(n_groups), each = per)
    cluster <- rep(sample(1:2, n_groups, TRUE), each = per)
    z <- rep(stats::rnorm(n_groups), each = per)
    class <- ifelse(stats::runif(length(g)) <
                      stats::plogis(ifelse(cluster == 1L, 1, -1) + 0.8 * z),
                    1L, 2L)
    draw <- function(p1, p2) {
      ifelse(stats::runif(length(g)) < ifelse(class == 1L, p1, p2), 1L, 2L)
    }
    data <- data.frame(g = g, z = z, y1 = draw(0.85, 0.2), y2 = draw(0.8, 0.25),
                       y3 = draw(0.75, 0.3), y4 = draw(0.2, 0.7),
                       y5 = draw(0.9, 0.35))
    data$y6 <- ifelse(data$y1 == 1L, 1L, sample(1:2, nrow(data), TRUE))
    data
  })
}

#' Two-level data for pooling (tests' `.pool_fixture()`): groups `g`, a row
#' covariate `z` predicting the profile, and two indicators
golden_svc_pool_frame <- function() {
  golden_svc_with_seed(5L, {
    n_groups <- 40L
    size <- 8L
    n <- n_groups * size
    frame <- data.frame(g = rep(seq_len(n_groups), each = size),
                        z = stats::rnorm(n))
    state <- 2L - stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z))
    frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
    frame$y2 <- stats::rnorm(n, c(0, 2)[state])
    frame
  })
}

#' Stand-in imputations: the same data with the first 20 covariate values
#' redrawn in each
golden_svc_pool_imputed <- function(frame, m = 3L) {
  golden_svc_with_seed(1L, lapply(seq_len(m), function(index) {
    copy <- frame
    copy$z[seq_len(20L)] <- stats::rnorm(20L)
    copy
  }))
}

#' Binary `student_esm` items of the first 30 students
golden_svc_esm_items <- c("time_with_friends", "on_social_media",
                          "tv_video_games", "listened_music")

golden_cases_services <- function() {
  engagement <- golden_svc_engagement()
  vars <- golden_engagement_vars
  srl_small <- golden_svc_srl()
  srl_vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy")
  esm <- golden_esm()
  lca_boundary <- golden_svc_lca_boundary()
  pool_frame <- golden_svc_pool_frame()
  srl_missing <- srl_small
  srl_missing$cognitive_strategies[seq(3L, 150L, by = 11L)] <- NA
  srl_missing$self_efficacy[seq(7L, 150L, by = 13L)] <- NA

  two_level <- function() {
    multilpa(engagement, vars, id = "student", n_profiles = 2,
             n_group_classes = 2, n_starts = 2, seed = 1)
  }
  single_level <- function(n_profiles = 2) {
    lpa(srl_small, srl_vars, n_profiles = n_profiles, n_starts = 2, seed = 1)
  }
  lca_fit <- function(n_classes) {
    lca(esm, golden_svc_esm_items, n_classes = n_classes, n_starts = 2,
        seed = 1)
  }

  list(
    # -- enumeration --------------------------------------------------------
    golden_case(
      "svc_enumerate_lpa_basic", "enumeration",
      "enumerate_lpa(): 1-3 profiles x the four basic structures (EEI/VVI/EEE/VVV)",
      fit = function() enumerate_lpa(srl_small, srl_vars, n_profiles = 1:3,
                                     n_starts = 2, seed = 1),
      data = srl_small),
    golden_case(
      "svc_enumerate_lpa_all_structures", "enumeration",
      "enumerate_lpa(model = \"all\"): all 14 structures at 2 profiles",
      fit = function() enumerate_lpa(srl_small, srl_vars, n_profiles = 2,
                                     model = "all", n_starts = 2, seed = 1),
      data = srl_small),
    golden_case(
      "svc_enumerate_lca", "enumeration",
      "enumerate_lca(): 1-3 classes on four binary ESM items",
      fit = function() enumerate_lca(esm, golden_svc_esm_items,
                                     n_classes = 1:3, n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "svc_enumerate_classes_two_level", "enumeration",
      "enumerate_classes(): two-level grid, 1-2 profiles x 1-2 group classes x EEI/VVI",
      fit = function() enumerate_classes(engagement, vars, id = "student",
                                         n_profiles = 1:2, n_group_classes = 1:2,
                                         model = c("EEI", "VVI"), n_starts = 2,
                                         seed = 1),
      data = engagement),
    golden_case(
      "svc_candidate_fit", "enumeration",
      "candidate_fit() picks the 2x2 VVI fit out of a two-level grid",
      fit = function() {
        grid <- enumerate_classes(engagement, vars, id = "student",
                                  n_profiles = 1:2, n_group_classes = 1:2,
                                  model = c("EEI", "VVI"), n_starts = 2,
                                  seed = 1)
        candidate_fit(grid, n_profiles = 2, n_group_classes = 2, model = "VVI")
      },
      data = engagement),
    golden_case(
      "svc_candidate_fit_ambiguous", "enumeration",
      "candidate_fit() by counts alone over a multi-structure grid is refused",
      fit = function() {
        grid <- enumerate_classes(engagement, vars, id = "student",
                                  n_profiles = 1:2, n_group_classes = 1,
                                  model = c("EEI", "VVI"), n_starts = 2,
                                  seed = 1)
        candidate_fit(grid, n_profiles = 2)
      },
      data = engagement),

    # -- bootstrap likelihood-ratio tests -----------------------------------
    golden_case(
      "svc_bootstrap_lrt_lpa", "bootstrap_lrt",
      "bootstrap_lrt(): single-level LPA 2 vs 3 profiles, 5 replicates",
      fit = function() single_level(3), data = srl_small,
      services = "bootstrap_lrt", null_fit = function() single_level(2)),
    golden_case(
      "svc_bootstrap_lrt_lca", "bootstrap_lrt",
      "bootstrap_lrt(): LCA 1 vs 2 classes on binary items, 5 replicates",
      fit = function() lca_fit(2), data = esm,
      services = "bootstrap_lrt", null_fit = function() lca_fit(1)),
    golden_case(
      "svc_bootstrap_lrt_fiml", "bootstrap_lrt",
      "bootstrap_lrt() with missing = \"fiml\": LPA 1 vs 2 profiles",
      fit = function() lpa(srl_missing, srl_vars, n_profiles = 2,
                           missing = "fiml", n_starts = 2, seed = 1),
      data = srl_missing, services = "bootstrap_lrt",
      null_fit = function() lpa(srl_missing, srl_vars, n_profiles = 1,
                                missing = "fiml", n_starts = 2, seed = 1)),
    golden_case(
      "svc_bootstrap_lrt_profiles", "bootstrap_lrt",
      "bootstrap_lrt(): two-level 1 vs 2 profiles (one group class)",
      fit = function() multilpa(engagement, vars, id = "student",
                                n_profiles = 2, n_group_classes = 1,
                                n_starts = 2, seed = 1),
      data = engagement, services = "bootstrap_lrt",
      null_fit = function() multilpa(engagement, vars, id = "student",
                                     n_profiles = 1, n_group_classes = 1,
                                     n_starts = 2, seed = 1)),

    golden_case(
      "svc_bootstrap_lrt_covariates", "bootstrap_lrt",
      "bootstrap_lrt(): covariate model (profile covariate z), 2 vs 3 profiles",
      fit = function() multilpa(pool_frame, c("y1", "y2"), id = "g",
                                n_profiles = 3, n_group_classes = 1,
                                profile_covariates = "z", n_starts = 2,
                                seed = 1),
      data = pool_frame, services = "bootstrap_lrt",
      null_fit = function() multilpa(pool_frame, c("y1", "y2"), id = "g",
                                     n_profiles = 2, n_group_classes = 1,
                                     profile_covariates = "z", n_starts = 2,
                                     seed = 1)),
    golden_case(
      "svc_bootstrap_lrt_mixed_families", "bootstrap_lrt",
      "bootstrap_lrt() of a plain null against a covariate alternative: refused (currently by an unclassed error)",
      fit = function() multilpa(pool_frame, c("y1", "y2"), id = "g",
                                n_profiles = 2, n_group_classes = 1,
                                profile_covariates = "z", n_starts = 2,
                                seed = 1),
      data = pool_frame, services = "bootstrap_lrt",
      null_fit = function() multilpa(pool_frame, c("y1", "y2"), id = "g",
                                     n_profiles = 1, n_group_classes = 1,
                                     n_starts = 2, seed = 1)),
    golden_case(
      "svc_bootstrap_lrt_fixed", "bootstrap_lrt",
      "bootstrap_lrt(): fixed measurement, 1 vs 2 group classes (fit_staged alternative, fixed = \"measurement\" null)",
      fit = function() fit_staged(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  n_starts = 2, seed = 1),
      data = engagement, services = "bootstrap_lrt",
      null_fit = function() {
        staged <- fit_staged(engagement, vars, id = "student", n_profiles = 2,
                             n_group_classes = 2, n_starts = 2, seed = 1)
        multilpa(engagement, vars, id = "student", n_profiles = 2,
                 n_group_classes = 1, fixed = "measurement",
                 start = starting_values(staged, what = "measurement"),
                 n_starts = 2, seed = 1)
      }),
    golden_case(
      "svc_fit_staged_one_class", "fit_staged",
      "fit_staged(n_group_classes = 1) is refused (currently by an unclassed error)",
      fit = function() fit_staged(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 1,
                                  n_starts = 2, seed = 1),
      data = engagement),

    # -- parameter inference --------------------------------------------------
    golden_case(
      "svc_inference_lpa", "inference",
      "parameter_inference(): Wald observed/robust/OPG and bootstrap (5) on a single-level LPA",
      fit = function() single_level(2), data = srl_small,
      services = c("inference_opg", "inference_bootstrap")),
    golden_case(
      "svc_inference_lpa_full", "inference",
      "parameter_inference(): Wald and bootstrap on a VVV single-level LPA",
      fit = function() lpa(srl_small, srl_vars, n_profiles = 2, model = "VVV",
                           n_starts = 2, seed = 1),
      data = srl_small, services = c("inference_opg", "inference_bootstrap")),
    golden_case(
      "svc_inference_two_level", "inference",
      "parameter_inference(): Wald (all three), bootstrap and boundary = \"fix\" on a two-level fit",
      fit = two_level, data = engagement,
      services = c("inference_opg", "inference_bootstrap", "inference_fix")),
    golden_case(
      "svc_inference_boundary_lca", "inference",
      "a covariate LCA with a response probability at its floor: Wald refused, boundary = \"fix\" holds it",
      fit = function() multilca(lca_boundary, vars = paste0("y", 1:6), id = "g",
                                n_profiles = 2, n_group_classes = 2,
                                profile_covariates = "z", n_starts = 2, seed = 1,
                                select_start = "converged"),
      data = lca_boundary, services = c("inference_fix", "inference_opg")),
    golden_case(
      "svc_inference_bootstrap_lca", "inference",
      "parameter_inference(method = \"bootstrap\") on a categorical LCA (5 replicates)",
      fit = function() lca_fit(2), data = esm,
      services = "inference_bootstrap"),
    golden_case(
      "svc_inference_bootstrap_mixed_extra", "inference",
      "bootstrap inference on a mixed ordinal + count (Poisson) + continuous LPA",
      fit = function() lpa(esm, c("happy", "worried", "exhausted"),
                           n_profiles = 2, ordinal = "happy", count = "worried",
                           n_starts = 2, seed = 1),
      data = esm, services = "inference_bootstrap"),
    golden_case(
      "svc_inference_bootstrap_vei", "inference",
      "bootstrap inference on a constrained VEI structure",
      fit = function() lpa(srl_small, srl_vars, n_profiles = 2, model = "VEI",
                           n_starts = 2, seed = 1),
      data = srl_small, services = c("inference_bootstrap", "inference_opg")),
    golden_case(
      "svc_inference_bootstrap_staged", "inference",
      "bootstrap inference on a fit_staged() fit (measurement held)",
      fit = function() fit_staged(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  n_starts = 2, seed = 1),
      data = engagement, services = "inference_bootstrap"),
    golden_case(
      "svc_inference_adjust_level", "inference",
      "parameter_inference(level = 0.9, adjust = \"holm\", step = 1e-5)",
      fit = function() {
        parameter_inference(two_level(), level = 0.9, adjust = "holm",
                            step = 1e-5)
      },
      data = engagement),

    # -- three-step and R3STEP ------------------------------------------------
    golden_case(
      "svc_three_step_two_level", "three_step",
      "three_step(): BCH, proportional, modal, pairwise contrasts and a group-level outcome",
      fit = two_level, data = engagement,
      services = c("three_step_bch", "three_step_proportional",
                   "three_step_modal", "three_step_pairs", "three_step_groups"),
      outcome = "previous_grade", group_outcome = "grade_mean"),
    golden_case(
      "svc_three_step_independent", "three_step",
      "three_step(vcov_type = \"independent\", ci_level = 0.9, adjust = \"holm\") with contrasts",
      fit = function() {
        three_step(two_level(), data = engagement, outcome = "previous_grade",
                   vcov_type = "independent", ci_level = 0.9,
                   contrast = "pairs", adjust = "holm")
      },
      data = engagement),
    golden_case(
      "svc_three_step_single_level", "three_step",
      "three_step() on a single-level 3-profile LPA: BCH, modal, pairs",
      fit = function() single_level(3), data = srl_small,
      services = c("three_step_bch", "three_step_modal", "three_step_pairs"),
      outcome = "test_anxiety"),
    golden_case(
      "svc_three_step_indicator_refused", "three_step",
      "three_step() with an outcome that is one of the indicators is refused",
      fit = function() {
        three_step(two_level(), data = engagement, outcome = "browse")
      },
      data = engagement),
    golden_case(
      "svc_r3step_two_level", "r3step",
      "r3step(): individual-level and group-level covariates",
      fit = two_level, data = engagement,
      services = c("r3step", "r3step_groups"),
      covariates = "previous_grade", group_covariates = "grade_mean"),
    golden_case(
      "svc_r3step_by_group_class", "r3step",
      "r3step(by_group_class = TRUE)",
      fit = function() {
        r3step(two_level(), data = engagement, covariates = "previous_grade",
               by_group_class = TRUE)
      },
      data = engagement),
    golden_case(
      "svc_r3step_observed", "r3step",
      "r3step(vcov_type = \"observed\", adjust = \"holm\", ci_level = 0.9)",
      fit = function() {
        r3step(two_level(), data = engagement, covariates = "previous_grade",
               vcov_type = "observed", adjust = "holm", ci_level = 0.9)
      },
      data = engagement),

    # -- sensitivity ---------------------------------------------------------
    golden_case(
      "svc_sensitivity_two_level", "sensitivity",
      "sensitivity() over seeds 1:2 on a two-level fit",
      fit = two_level, data = engagement, services = "sensitivity"),
    golden_case(
      "svc_sensitivity_options", "sensitivity",
      "sensitivity(seeds = 3:4, n_starts = 1, max_iter = 50, tol = 1e-6) on an LPA",
      fit = function() {
        sensitivity(single_level(3), seeds = 3:4, n_starts = 1, max_iter = 50,
                    tol = 1e-6)
      },
      data = srl_small),

    # -- multiple imputation -------------------------------------------------
    golden_case(
      "svc_pool_covariates", "pool_imputations",
      "pool_imputations(): three hand-made imputations of a covariate model (Rubin's rules)",
      fit = function() {
        pool_imputations(golden_svc_pool_imputed(pool_frame), c("y1", "y2"),
                         "g", 2L, n_group_classes = 1L,
                         profile_covariates = "z", n_starts = 2, seed = 1)
      },
      data = pool_frame),
    golden_case(
      "svc_pool_robust", "pool_imputations",
      "pool_imputations(vcov_type = \"robust\", level = 0.9) of a plain two-level LPA",
      fit = function() {
        pool_imputations(golden_svc_pool_imputed(pool_frame), c("y1", "y2"),
                         "g", 2L, n_group_classes = 1L, n_starts = 2, seed = 1,
                         vcov_type = "robust", level = 0.9)
      },
      data = pool_frame),
    golden_case(
      "svc_pool_mice", "pool_imputations",
      "pool_imputations() of a seeded mice::mice() `mids` (pmm, m = 3); covariate with missing values",
      fit = function() {
        if (!requireNamespace("mice", quietly = TRUE)) return("<mice not installed>")
        frame <- pool_frame
        frame$z[seq(2L, nrow(frame), by = 9L)] <- NA
        imputed <- mice::mice(frame[c("g", "z", "y1", "y2")], m = 3,
                              method = c("", "pmm", "", ""), seed = 1,
                              printFlag = FALSE)
        pool_imputations(imputed, c("y1", "y2"), "g", 2L, n_group_classes = 1L,
                         profile_covariates = "z", n_starts = 2, seed = 1)
      },
      data = pool_frame),

    # -- staged fits and starting values ---------------------------------------
    golden_case(
      "svc_fit_staged", "fit_staged",
      "fit_staged(): measurement at one group class, then two group classes",
      fit = function() fit_staged(engagement, vars, id = "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  n_starts = 2, seed = 1),
      data = engagement,
      services = c("starting_values", "starting_values_measurement")),
    golden_case(
      "svc_fit_staged_categorical", "fit_staged",
      "fit_staged() with categorical indicators and full covariance continuous parts",
      fit = function() fit_staged(esm, c(golden_svc_esm_items[1:3], "happy",
                                         "relaxed"),
                                  id = "student", n_profiles = 2,
                                  n_group_classes = 2,
                                  categorical = golden_svc_esm_items[1:3],
                                  n_starts = 2, seed = 1),
      data = esm),
    golden_case(
      "svc_starting_values_full", "starting_values",
      "starting_values(): full covariance, covariance auto/drop/keep and measurement only",
      fit = function() {
        fit <- lpa(srl_small, srl_vars, n_profiles = 2, model = "VVV",
                   n_starts = 2, seed = 1)
        list(auto = get_results(starting_values(fit), "all"),
             drop = get_results(starting_values(fit, covariance = "drop"),
                                "all"),
             keep = get_results(starting_values(fit, covariance = "keep"),
                                "all"),
             measurement = get_results(starting_values(fit,
                                                       what = "measurement"),
                                       "all"))
      },
      data = srl_small),

    # -- diagnostics, descriptives, prediction, report -------------------------
    golden_case(
      "svc_diagnostics_report", "diagnostics",
      "diagnostics(), descriptives(), report() and predict() on a two-level fit",
      fit = two_level, data = engagement,
      services = c("diagnostics", "descriptives", "report", "predict")),
    golden_case(
      "svc_diagnostics_overall", "diagnostics",
      "diagnostics(by = \"overall\") on a single-level LPA",
      fit = function() diagnostics(single_level(3), plots = FALSE,
                                   by = "overall"),
      data = srl_small),
    golden_case(
      "svc_predict_types", "predict",
      "predict(type = \"posterior\") and predict(type = \"density\") on new rows",
      fit = function() {
        fit <- two_level()
        newdata <- subset(golden_engagement(30L), student > 20)
        list(class = stats::predict(fit, newdata = newdata),
             posterior = stats::predict(fit, newdata = newdata,
                                        type = "posterior"),
             density = stats::predict(fit, newdata = newdata,
                                      type = "density"))
      },
      data = engagement),
    golden_case(
      "svc_descriptives_data", "descriptives",
      "descriptives() of a data frame: overall, with id (ICC), and by a factor",
      fit = function() {
        list(plain = descriptives(engagement, vars = vars),
             groups = descriptives(engagement, vars = vars, id = "student"),
             by = descriptives(engagement, vars = vars, by = "student_type"))
      },
      data = engagement)
  )
}
