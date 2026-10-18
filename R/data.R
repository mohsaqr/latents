#' Student engagement across a sequence of courses
#'
#' A simulated dataset of course enrolments nested in students. Each row is
#' one student's enrolment in one course. Students take 10 to 15 courses,
#' ordered by `sequence`, so the data support a cross-sectional two-level
#' model (`id = "student"`), a latent transition model (`time = "sequence"`)
#' and a membership-covariate model (`profile_covariates = "previous_grade"`).
#'
#' Enrolments follow one of two engagement patterns, and students are of two
#' kinds that differ in the mix of patterns: a committed student is engaged
#' in about four enrolments in five, a wavering student in about one in four.
#' The pooled engaged share describes neither kind; a two-level model recovers
#' the two kinds.
#'
#' The five activity indicators were simulated on the `log1p(events)` scale
#' of a logged activity export (no raw counts were generated) and then
#' standardised within course: each value is the enrolment's distance from
#' its course's mean in that course's standard deviations, so every indicator
#' has mean zero and standard deviation one within every course. A `browse`
#' of 1.2 is 1.2 standard deviations above the average enrolment in that
#' course. The simulated courses do not differ in activity volume, so the
#' standardisation sets the scale but removes no course effect. The log
#' transformation removes the right skew of raw activity counts, which
#' standardisation within course would leave unchanged.
#'
#' @section Local dependence:
#' `attendance` is generated from the four click measures, so the indicators
#' are not independent given the profile, as a latent profile model with a
#' diagonal residual covariance assumes. For the two-profile, two-group-class
#' fit of the example, `get_results(fit, "residuals")` reports a largest
#' residual of 0.23 (`lectures` with `attendance` within the first profile,
#' `p = 6.4e-12`), and `covariance_model = "full"` improves the log likelihood
#' from -8439.8 to -8313.1 for twenty more parameters, preferred by every
#' information criterion except the individual-level AWE. Under a diagonal
#' model, [enumerate_classes()] tends to favour more classes, which absorb
#' the residual association.
#'
#' @format A data frame with 1422 rows (106 students, 32 courses, up to 15
#'   enrolments each) and 11 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 106. The nesting unit;
#'     pass it as `id`.}
#'   \item{course}{Integer course identifier, 1 to 32.}
#'   \item{sequence}{Integer position of the course in the student's own
#'     order, 1 to 15. Pass it as `time`. Students who took fewer courses have
#'     shorter sequences, so the panel is unbalanced.}
#'   \item{browse}{Numeric, course page views on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{lectures}{Numeric, lecture videos viewed on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{forum_read}{Numeric, forum posts read on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{forum_post}{Numeric, forum posts written on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{attendance}{Numeric, days the student was active in the course,
#'     on the `log1p(events)` scale, standardized within course. Generated from
#'     the four click measures with added noise, so it correlates 0.51 to 0.67
#'     with them overall and 0.10 to 0.23 among the enrolments assigned to one
#'     profile of the example fit. See *Local dependence*.}
#'   \item{previous_grade}{Numeric, the grade the student earned in the
#'     previous course, standardised across the cohort; an entry grade at a
#'     student's first course. A predictor of engagement, for
#'     `multilpa(profile_covariates = )` or [r3step()].}
#'   \item{engagement}{Factor, `"disengaged"` or `"engaged"`: the pattern
#'     the row was generated from. Not for fitting; kept to check recovery.}
#'   \item{student_type}{Factor, `"committed"` or `"wavering"`: the
#'     student's generating kind, corresponding to the group classes. Not for
#'     fitting; kept to check recovery.}
#' }
#' @source Simulated by `data-raw/course-engagement.R`, seed 2026. The
#'   simulation is calibrated to aggregate constants of a real
#'   learning-analytics export (scales, spreads, correlations and the
#'   persistence of engagement across courses). No row, identifier or value
#'   of any real student is included, and the source data are not
#'   distributed.
#' @seealso `vignette("lpa", package = "latents")` for the analysis these data are used in,
#'   [multilpa()] for the two-level model, [lta()] for the
#'   sequence, and `multilpa(profile_covariates = )` for `previous_grade`.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' as.data.frame(fit)
#' get_results(fit, what = "profile_probabilities")
#' # The generating truth ships alongside, so recovery can be checked directly.
#' get_results(fit, what = "assignments", data = course_engagement,
#'          truth = "engagement")
"course_engagement"

#' Leisure activities of university students in daily life
#'
#' Experience-sampling data on the leisure activities of university students
#' between prompts at which they had not studied. Each row is one answered
#' non-study prompt, with all eight activity items and all four affect
#' ratings observed. Prompts are nested in students, so the data suit a
#' two-level latent class analysis.
#'
#' In the original study students answered six prompts a day for fourteen
#' days; the activity items were asked only at prompts where the student had
#' not studied since the previous prompt. The dataset holds 100 students,
#' drawn with a fixed seed from those with at least 15 such prompts, and every
#' such prompt they answered: 15 to 51 per student (median 24). The eight
#' activity items are categorical indicators; the four affect ratings can be
#' added as continuous indicators in a mixed measurement model. `worried` has
#' a floor (49.7% of ratings are 1), so a profile of low-worry prompts can
#' have zero variance on it and reach the `min_variance` bound.
#'
#' @format A data frame with 2582 rows (100 students) and 15 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 100. The nesting unit;
#'     pass it as `id`.}
#'   \item{day}{Integer study day, 0 to 13.}
#'   \item{beep}{Integer prompt of the day, 1 to 5 among the prompts kept.}
#'   \item{time_with_friends}{Factor, `"no"` or `"yes"`: spent time with
#'     friends since the previous prompt.}
#'   \item{on_social_media}{Factor, `"no"` or `"yes"`: used social media.}
#'   \item{tv_video_games}{Factor, `"no"` or `"yes"`: watched TV or played
#'     video games.}
#'   \item{listened_music}{Factor, `"no"` or `"yes"`: listened to music.}
#'   \item{sports}{Factor, `"no"` or `"yes"`: did sports.}
#'   \item{walking}{Factor, `"no"` or `"yes"`: went for a walk.}
#'   \item{reading}{Factor, `"no"` or `"yes"`: read.}
#'   \item{part_time_job}{Factor, `"no"` or `"yes"`: worked in a part-time
#'     job.}
#'   \item{happy}{Integer rating of feeling happy, 1 to 7.}
#'   \item{relaxed}{Integer rating of feeling relaxed, 1 to 7.}
#'   \item{worried}{Integer rating of feeling worried, 1 to 7.}
#'   \item{exhausted}{Integer rating of feeling exhausted, 1 to 7.}
#' }
#' @source openESM dataset 0062, deposited by Neubauer and Schmiedek at
#'   Zenodo, \doi{10.5281/zenodo.17347974}, under the Creative Commons
#'   Attribution 4.0 licence (CC-BY 4.0), and retrieved through the openESM
#'   database (<https://openesmdata.org>). The subset is produced by
#'   `data-raw/student-esm.R` in the source repository; identifiers are
#'   renumbered and the values are otherwise unchanged.
#' @references Neubauer, A. B., & Schmiedek, F. (2024). Approaching academic
#'   adjustment on multiple time scales. *Zeitschrift für
#'   Erziehungswissenschaft*, 27, 147--168. \doi{10.1007/s11618-023-01182-8}
#' @seealso [multilpa()] with `categorical =` for the two-level latent class
#'   model.
#' @examples
#' head(student_esm)
#' \donttest{
#' activities <- c("time_with_friends", "on_social_media", "tv_video_games",
#'                 "listened_music", "sports", "walking", "reading",
#'                 "part_time_job")
#' first_week <- subset(student_esm, day <= 6)
#' lca <- multilpa(first_week, vars = activities, id = "student",
#'                 n_profiles = 2, n_group_classes = 2,
#'                 categorical = activities, n_starts = 3, seed = 1)
#' get_results(lca, what = "responses")
#' get_results(lca, what = "profile_probabilities")
#' }
"student_esm"

#' Study hours and quiz scores under two study strategies
#'
#' A simulated dataset for [mixture_regression()]. Each row is one student's
#' study week: 150 students, six weeks each. The generating classes are
#' included for checking recovery.
#'
#' Each week a student uses one of two study strategies, which sets the
#' regression of quiz `score` on `hours`: a deep week gains about 4.5 points
#' per hour from a base of 35, a surface week under one point per hour from a
#' base of 55. Steady students study deeply in about 85% of weeks and erratic
#' students in about 25%; higher `motivation` makes a student more likely to
#' be steady, and more `sleep` makes a deep week more likely.
#'
#' `score ~ hours` is therefore a two-class mixture at the week level
#' (`class_level = "observation"`) with class shares that differ by student
#' (`id = "student"`, `n_group_classes = 2`). `passed` depends on the
#' student's kind, not the week's strategy, so it is a mixture at the student
#' level (`class_level = "group"`). `questions` is a count whose Poisson
#' regression depends on the week's strategy.
#'
#' @format A data frame with 900 rows (150 students by 6 weeks) and 10
#'   columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 150.}
#'   \item{week}{Integer study week, 1 to 6.}
#'   \item{hours}{Hours studied that week, 0 to 12.}
#'   \item{sleep}{Hours slept the night before the quiz.}
#'   \item{motivation}{Standardized motivation, constant within student.}
#'   \item{score}{Quiz score.}
#'   \item{passed}{Integer 0/1: passed the weekly check.}
#'   \item{questions}{Count of questions the student asked that week.}
#'   \item{strategy}{Factor, `"deep"` or `"surface"`: the week's generating
#'     strategy. Not for fitting; kept to check recovery.}
#'   \item{student_type}{Factor, `"steady"` or `"erratic"`: the student's
#'     generating kind. Not for fitting; kept to check recovery.}
#' }
#' @source Simulated by `data-raw/study-hours.R` in the source repository.
#' @seealso [mixture_regression()], [enumerate_regressions()].
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' get_results(fit, "coefficients")
"study_hours"

#' Self-regulated learning scores simulated by a large language model
#'
#' Scores on five constructs of the Motivated Strategies for Learning
#' Questionnaire (MSLQ) for 300 simulated respondents, one row per
#' respondent. Each score is the mean of the construct's Likert items, on a
#' scale from 1 to 7. The responses were **generated by a large language
#' model (GPT)** in a study of how language models answer psychological
#' surveys; they are not data on students.
#'
#' Four constructs correlate about 0.8 with each other; test anxiety
#' correlates weakly and negatively with them. The data illustrate how the
#' covariance structure affects a latent profile analysis.
#'
#' @format A data frame with 300 rows and 5 columns:
#' \describe{
#'   \item{`cognitive_strategies`}{cognitive strategy use}
#'   \item{`intrinsic_value`}{intrinsic value of the course}
#'   \item{`self_efficacy`}{self-efficacy}
#'   \item{`self_regulation`}{self-regulation}
#'   \item{`test_anxiety`}{test anxiety}
#' }
#' @source The `SRL_GPT` data of the psychnets package, with descriptive
#'   column names. Saqr, M. (2025). Delving into the psychology of machines:
#'   Exploring the structure of self-regulated learning via LLM-generated
#'   survey responses. *Computers in Human Behavior*, 173, 108769.
#'   \doi{10.1016/j.chb.2025.108769}
#' @examples
#' descriptives(srl, names(srl))
#' @keywords datasets
"srl"

#' Achievement growth of students with three kinds of trajectory
#'
#' A simulated longitudinal dataset for growth models. Each row is one
#' student's test score at one occasion: 300 students tested on up to seven
#' occasions (`wave` 0 to 6). Three kinds of student follow three
#' trajectories: *improving* (starts near 48, gains 3.5
#' points a wave), *stable* (starts near 62, gains 0.5) and *declining*
#' (starts near 58, loses 2.5). Within each kind, students vary around the
#' trajectory: each has their own start and rate of change (random intercept
#' SD 4, random slope SD 0.8, correlation -0.3, the same spread in every kind),
#' with residual SD 3. Motivation, measured once, makes a student more likely
#' to be improving and less likely to be declining. About 12% of tests were
#' missed at random, so the panel is unbalanced.
#'
#' Because students vary within their kind, a trajectory model without random
#' effects needs more than three classes; a growth mixture model with random
#' intercepts and slopes recovers three.
#'
#' @format A data frame with 1844 rows (300 students) and 5 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 300. Pass it as `id`.}
#'   \item{wave}{Integer test occasion, 0 to 6. The time variable.}
#'   \item{score}{Achievement score.}
#'   \item{motivation}{Standardized motivation, constant within student.}
#'   \item{trajectory}{Factor, `"improving"`, `"stable"` or `"declining"`:
#'     the student's generating kind. Not for fitting; kept to check
#'     recovery.}
#' }
#' @source Simulated by `data-raw/growth-scores.R` in the source repository.
#' @seealso [mixture_regression()], [compare_models()].
#' @examples
#' descriptives(growth_scores, c("score", "motivation"))
#' \donttest{
#' fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
#'                           id = "student", class_level = "group",
#'                           random = "wave", random_covariance = "equal",
#'                           n_starts = 4, seed = 1)
#' get_results(fit, "recovery", data = growth_scores, truth = "trajectory")
#' }
"growth_scores"

#' Reading growth of students nested in schools
#'
#' A simulated multilevel longitudinal dataset for the multilevel growth
#' mixture model. Each row is one student's reading score at one occasion:
#' 600 students in 40 schools (15 per school), tested on up to five occasions
#' (`wave` 0 to 4). Two kinds of student follow two
#' trajectories: *improving* (starts near 48, gains 3 points a wave) and
#' *plateauing* (starts near 56, gains 0.3). Students vary around their
#' trajectory (random intercept SD 3, random slope SD 0.6, correlation -0.2;
#' residual SD 2). Schools are of two types: in *supportive* schools 80% of
#' students improve, in *struggling* schools 25% do, and a school's support
#' programme (standardized, measured once) makes it more likely to be
#' supportive. About 8% of tests were missed at random.
#'
#' The school types correspond to the group classes of
#' `mixture_regression(..., cluster = "school", n_group_classes = 2)`, with
#' `group_membership = "programme"`.
#'
#' @format A data frame with 2780 rows (600 students, 40 schools) and 7
#'   columns:
#' \describe{
#'   \item{school}{Integer school identifier, 1 to 40. Pass it as `cluster`.}
#'   \item{student}{Integer student identifier, 1 to 600. Pass it as `id`.}
#'   \item{wave}{Integer test occasion, 0 to 4. The time variable.}
#'   \item{score}{Reading score.}
#'   \item{programme}{Standardized school support programme, constant within
#'     school.}
#'   \item{school_type}{Factor, `"supportive"` or `"struggling"`: the school's
#'     generating type. Not for fitting; kept to check recovery.}
#'   \item{trajectory}{Factor, `"improving"` or `"plateauing"`: the student's
#'     generating kind. Not for fitting; kept to check recovery.}
#' }
#' @source Simulated by `data-raw/growth-schools.R` in the source repository.
#' @seealso [mixture_regression()], [get_results()].
#' @examples
#' descriptives(growth_schools, c("score", "programme"))
#' \donttest{
#' fit <- mixture_regression(score ~ wave, growth_schools, n_classes = 2,
#'                           id = "student", class_level = "group",
#'                           random = "wave", random_covariance = "equal",
#'                           cluster = "school", n_group_classes = 2,
#'                           group_membership = "programme",
#'                           n_starts = 2, seed = 1)
#' get_results(fit, "group_classes")
#' get_results(fit, "recovery", data = growth_schools, truth = "trajectory")
#' }
#' @keywords datasets
"growth_schools"
