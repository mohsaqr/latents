# The identity panel for `search = "kmeans"`: one case per fit path and
# control that the start search touches. Each case is a function of nothing
# that returns a fitted object; record.R runs every case under a library and
# keeps the result, the RNG state on exit and the warnings.
#
# A case leaves `search` unset, so the same panel records the released
# behaviour (no such argument) and, once emEM is the default, is run again
# with `options(latents.panel_search = "kmeans")` read by `with_search()`.

with_search <- function(call_args) {
  search <- getOption("latents.panel_search")
  if (is.null(search)) call_args else c(call_args, list(search = search))
}

fit <- function(verb, ...) do.call(verb, with_search(list(...)))

esm_items <- c("time_with_friends", "on_social_media", "tv_video_games",
               "listened_music", "sports", "walking", "reading", "part_time_job")
engagement_items <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
srl_items <- names(srl)

# Draws without touching the caller's stream, so building the panel data does
# not move the RNG state a case starts from.
withr_free_sample <- function(n, size, seed) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  sample.int(n, size)
}

with_holes <- function(d) {
  holes <- withr_free_sample(nrow(d), 30L, 11L)
  d$self_efficacy[holes] <- NA
  d$test_anxiety[rev(holes)[seq_len(10L)]] <- NA
  d
}

# A group-level covariate: the number of rows of each course, in tens.
with_course_size <- function(d) {
  sizes <- table(d$course)
  d$course_size <- as.numeric(sizes[as.character(d$course)]) / 10
  d
}

srl_missing <- with_holes(srl)
srl_weighted <- transform(srl, w = rep(c(0.5, 1, 1.5, 2), length.out = nrow(srl)))
engagement_groups <- with_course_size(course_engagement)

cases <- list(
  lpa_default = function() fit(lpa, srl, srl_items, 3, seed = 1),
  lpa_seed_null = function() fit(lpa, srl, srl_items, 3),
  lpa_vvv = function() fit(lpa, srl, srl_items, 3, model = "VVV", seed = 2),
  lpa_eee_k1 = function() fit(lpa, srl, srl_items, 1, model = "EEE", seed = 2),
  lpa_noise = function() fit(lpa, srl, srl_items, 2, noise = TRUE, seed = 3),
  lpa_fiml = function() fit(lpa, srl_missing, srl_items, 2, missing = "fiml", seed = 4),
  lpa_weights = function() fit(lpa, srl_weighted, srl_items, 2, weights = "w", seed = 5),
  lpa_prior = function() fit(lpa, srl, srl_items, 3, prior = prior_control(), seed = 6),
  lpa_no_acceleration = function() fit(lpa, srl, srl_items, 3, acceleration = "none", seed = 7),
  lpa_converged_exclude = function() fit(lpa, srl, srl_items, 4, model = "VVV",
                                         select_start = "converged",
                                         boundary = "exclude", seed = 8),
  lpa_tiny_budget = function() fit(lpa, srl, srl_items, 3, max_iter = 3, seed = 9),
  lpa_cores = function() fit(lpa, srl, srl_items, 3, cores = 2, seed = 10),
  lpa_volume = function() fit(lpa, srl, srl_items, 3, volume = "equal",
                              shape = "varying", orientation = "axis", seed = 11),
  lca_esm = function() fit(lca, student_esm, esm_items, 3, seed = 12),
  lca_local_dependence = function() fit(lca, student_esm, esm_items, 2,
                                        local_dependence = "sports:walking", seed = 13),
  mixed_mode = function() fit(lpa, student_esm, c("happy", "relaxed", "sports", "reading"),
                              2, categorical = c("sports", "reading"), seed = 14),
  ordinal_count = function() fit(lpa, study_hours, c("hours", "questions", "passed"), 2,
                                 ordinal = "passed", count = "questions", seed = 15),
  direct_effects = function() fit(lpa, srl, setdiff(srl_items, "self_efficacy"), 2,
                                  direct_effects = "test_anxiety ~ self_efficacy", seed = 16),
  multilpa_two_level = function() fit(multilpa, course_engagement, engagement_items,
                                      "course", 3, 2, seed = 17),
  multilpa_two_level_null = function() fit(multilpa, course_engagement, engagement_items,
                                           "course", 3, 2),
  multilca_two_level = function() fit(multilca, student_esm, esm_items, "student", 3, 2,
                                      seed = 18),
  profile_covariates = function() fit(multilpa, course_engagement, engagement_items,
                                      "course", 2, 2, profile_covariates = "previous_grade",
                                      seed = 19),
  group_covariates = function() fit(multilpa, engagement_groups, engagement_items,
                                    "course", 2, 2, group_covariates = "course_size",
                                    seed = 20),
  additive = function() fit(multilpa, course_engagement, engagement_items, "course",
                            n_group_classes = 2, family = "additive", seed = 21),
  dispersion = function() fit(multilpa, course_engagement, engagement_items, "course",
                              n_group_classes = 2, family = "dispersion", seed = 22),
  restricted_cross_level = function() fit(multilpa, course_engagement, engagement_items,
                                          "course", 2, 2,
                                          family = "restricted_cross_level", seed = 23),
  full_cross_level = function() fit(multilpa, course_engagement, engagement_items,
                                    "course", 2, 2, family = "full_cross_level",
                                    seed = 24),
  transitions = function() fit(multilpa, course_engagement, engagement_items, "student",
                               2, 1, time = "sequence", seed = 25),
  lta = function() fit(lta, course_engagement, engagement_items, "student", 2,
                       time = "sequence", seed = 26),
  lta_mover_stayer = function() fit(lta, course_engagement, engagement_items, "student",
                                    2, time = "sequence", mover_stayer = TRUE,
                                    n_starts = 3, seed = 27),
  staged = function() fit(fit_staged, course_engagement, engagement_items, "course",
                          n_profiles = 2, n_group_classes = 2, seed = 28),
  supplied_start = function() {
    first <- fit(lpa, srl, srl_items, 3, seed = 29, n_starts = 2)
    fit(lpa, srl, srl_items, 3, start = starting_values(first), seed = 30)
  },
  zero_budget = function() {
    first <- fit(lpa, srl, srl_items, 2, seed = 31, n_starts = 2)
    fit(lpa, srl, srl_items, 2, start = starting_values(first), max_iter = 0)
  },
  enumeration = function() fit(enumerate_lpa, srl, srl_items, n_profiles = 1:3,
                               seed = 32),
  bootstrap = function() {
    one <- fit(lpa, srl, srl_items, 1, model = "EEE", seed = 33)
    two <- fit(lpa, srl, srl_items, 2, model = "EEE", seed = 33)
    bootstrap_lrt(one, two, iter = 9, seed = 34)
  },
  sensitivity = function() {
    three <- fit(lpa, srl, srl_items, 3, seed = 35)
    sensitivity(three, seeds = 1:3)
  }
)
