#' Student engagement across a sequence of courses
#'
#' A simulated dataset with the structure this package exists for, and all of
#' it in one place: enrolments nested in students, ordered within each student,
#' with a covariate that predicts which pattern an enrolment follows.
#'
#' Each row is one student's enrolment in one course. A student takes up to
#' fifteen courses, and `sequence` orders them, so the same rows support a
#' cross-sectional two-level model (`id = "student"`), a latent transition model
#' (`time = "sequence"`), and a membership-covariate model
#' (`profile_covariates = "previous_grade"`).
#'
#' Two kinds of student were generated. Both contain both patterns of
#' engagement; they differ in the *mix*, which is the two-level quantity this
#' package estimates. A committed student is engaged in about four enrolments in
#' five, a wavering student in about one in four. Pooled across the sample the
#' engaged share describes neither kind, which is the point: recovering the
#' split is what a two-level model does and a pooled one cannot.
#'
#' @section Local independence does not hold here, on purpose:
#' A latent profile model with a diagonal residual covariance assumes the
#' indicators are independent given the profile. `attendance` is built from the
#' click measures, so for these data that assumption is false, and the package
#' finds it: `get_results(fit, "residuals")` reports a largest residual of 0.21
#' `p = 1.7e-15`, and `covariance_model = "full"` improves the log likelihood
#' from -7618.5 to -7494.7 for twenty more parameters, which every information
#' criterion prefers.
#'
#' This is deliberate. A bundled dataset that satisfies local independence
#' exactly can only ever show those verbs agreeing that nothing is wrong; here
#' the diagnostic has a real violation to find, its cause is known, and the fix
#' can be seen to work. One consequence is worth knowing before it surprises
#' you: under a diagonal model [enumerate_classes()] keeps preferring more
#' classes, because extra classes absorb residual association the model is not
#' allowed to express directly.
#'
#' The five activity indicators are **standardized within course**: each value is
#' that enrolment's distance from its own course's mean, in that course's
#' standard deviations. Every indicator therefore has mean zero and standard
#' deviation one within every course, and a `browse` of 1.2 reads as "1.2
#' standard deviations more active than the average enrolment in that course".
#' This is the scale a learning-analytics reader meets, because activity volume
#' is not comparable across courses: a course with weekly quizzes generates more
#' clicks than a seminar whatever its students are like, so an analyst compares a
#' student with their coursemates rather than with the whole cohort.
#'
#' In this simulation the courses do not in fact differ in activity volume, so
#' the standardization sets the scale rather than removing a course effect. It is
#' stated here rather than left implicit because a reader who assumed otherwise
#' would misread what the `course` column does.
#'
#' The indicators are simulated directly on the `log1p(events)` scale, the scale
#' of a logged activity export; no raw counts are generated. They are
#' standardized on that scale, not as counts.
#' Standardizing is affine within a course, so it leaves the within-course skew
#' exactly as it found it, and raw learning-analytics counts are strongly
#' right-skewed: on an unlogged draft of these data the skew alone looked like an
#' extra class, exactly as a floor or ceiling effect does. The log removes the
#' skew; the standardization puts the courses on a common footing afterwards.
#'
#' @format A data frame with 1422 rows (106 students, 32 courses, up to 15
#'   enrolments each) and 11 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 106. The nesting unit;
#'     pass it as `id`.}
#'   \item{course}{Integer course identifier, 1 to 32. Which course the
#'     enrolment was in.}
#'   \item{sequence}{Integer, the position of this course in the student's own
#'     order, 1 upwards. Pass it as `time`. Students who took fewer courses have
#'     shorter sequences, so the panel is ragged.}
#'   \item{browse}{Numeric, course page views on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{lectures}{Numeric, lecture videos viewed on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{forum_read}{Numeric, forum posts read on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{forum_post}{Numeric, forum posts written on the `log1p(events)` scale, standardized
#'     within course.}
#'   \item{attendance}{Numeric, days the student was active in the course as
#'     on the `log1p(events)` scale, standardized within course. It is **generated from the four click measures** and then
#'     blurred, because a student is recorded present on a day precisely because
#'     they clicked something, and a day is a coarse unit -- one tick however
#'     much happened inside it. The indicators are therefore locally dependent
#'     on purpose: `attendance` correlates 0.51 to 0.67 with the click measures
#'     overall and 0.16 to 0.26 within a single profile, which is more
#'     association than a shared profile alone would produce. See the note on
#'     local dependence below.}
#'   \item{previous_grade}{Numeric, the grade the student earned in the course
#'     *before* this one, standardised across the cohort. At a student's first
#'     course an entry grade stands in its place. It is a predictor of
#'     engagement here, not an outcome of it: pass it as `profile_covariates`
#'     to `multilpa(profile_covariates = )`, or as the predictor to
#'     [r3step()].}
#'   \item{engagement}{Factor, `"disengaged"` or `"engaged"`, the pattern the
#'     row was generated from. This is the truth behind the data, not something
#'     a study would observe; it is shipped so a fitted model can be checked
#'     against what produced it. It is deliberately not called `profile`:
#'     `get_results(fit, "assignments")` adds a column of that name, and a
#'     dataset already
#'     carrying it would make the verb refuse to join its own output.}
#'   \item{student_type}{Factor, `"committed"` or `"wavering"`, the kind of
#'     student the row's student was generated as. The truth behind the group
#'     classes, shipped for the same reason.}
#' }
#' @source Simulated by `data-raw/course-engagement.R`, seed 2026. The
#'   simulation is calibrated to the shape of a real learning-analytics export
#'   -- its scales, its spreads, its correlations and how persistent engagement
#'   is across a student's courses -- but only those aggregate constants were
#'   used. No row, identifier or value of any real student is present, and the
#'   source is not distributed.
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
#' Experience-sampling data on what university students did between prompts
#' when they were not studying. Prompts are nested within students, so the
#' data suit a two-level latent class analysis: each prompt belongs to an
#' activity profile, and students differ in how often their prompts fall in
#' each profile.
#'
#' The students answered six prompts a day for fourteen days. At a prompt where
#' the student had not studied since the previous prompt, the questionnaire
#' asked which of eight leisure activities they had done. The study asked these
#' items only at non-study prompts, so the dataset contains those prompts alone.
#' Each row is one answered non-study prompt with all eight activity items and
#' all four affect ratings observed.
#'
#' The dataset is a subset of the original study: 100 of its students, drawn
#' with a fixed seed from the students with at least 15 such prompts, with every
#' such prompt they answered. Each student contributes 15 to 51 prompts (median
#' 24). The eight activity items are categorical indicators; the four affect
#' ratings are continuous-scored and can be combined with them in a mixed
#' measurement model. `worried` has a floor: 49.7% of its ratings are 1, so a
#' profile of low-worry prompts can have zero variance on it and reach the
#' `min_variance` bound.
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
#' A simulated dataset for [mixture_regression()], built so that every nesting of the
#' mixture regression has a model that generated some of its columns, and the
#' truth is kept for checking recovery.
#'
#' 150 students each report six study weeks. In every week a student studies
#' with one of two strategies, and the strategy decides how hours become a quiz
#' `score`: a deep week gains about 4.5 points per hour from a base of 35, a
#' surface week under one point per hour from a base of 55. Students are of two
#' kinds that differ in their *mix* of strategies -- steady students study
#' deeply in about 85% of weeks, erratic ones in about 25% -- and higher
#' `motivation` makes a student more likely to be steady. More `sleep` makes a
#' deep week more likely.
#'
#' So `score ~ hours` is a two-class mixture at the level of the week
#' (`class_level = "observation"`), whose class shares differ by student (the
#' two-level model, `id = "student"`, `n_group_classes = 2`). `passed` depends
#' on the student's kind, not the week's strategy, so it is a mixture at the
#' level of the student (`class_level = "group"`). `questions` is a count whose
#' Poisson regression depends on the week's strategy.
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
#' Questionnaire (MSLQ) for 300 simulated respondents. Each score is the mean
#' of the construct's Likert items, on a scale from 1 to 7. The responses were
#' **generated by a large language model (GPT)**, not collected from students:
#' they come from a study of how language models answer psychological
#' surveys, and are used here because they are a well-structured, openly
#' shareable cross-sectional example. They should not be read as evidence
#' about students.
#'
#' The five constructs are strongly correlated (about 0.8, except test
#' anxiety), which makes the data a good illustration of how the choice of
#' covariance structure shapes a latent profile analysis.
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
