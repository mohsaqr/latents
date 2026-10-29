#' Student engagement across a sequence of courses
#'
#' @description
#' `course_engagement` is a simulated dataset for illustrating student
#' engagement across successive courses. Each row represents one student's
#' enrolment in one course, with courses ordered within students. The data
#' include five activity indicators, standardised within course, and a
#' simulated previous-course grade.
#'
#' The data were generated from two engagement profiles and two student-level
#' classes that differ in the prevalence of those profiles. Engagement can
#' change across a student's course sequence. The generating profile and
#' student-class labels are retained, allowing fitted classifications to be
#' compared with the structure used to generate the data.
#'
#' The dataset supports examples of multilevel profile analysis, latent
#' transition analysis and analyses relating previous grades to profile
#' membership.
#'
#' @details
#' Students take 10 to 15 courses, ordered by `sequence`. Use `id = "student"`
#' for the nesting structure, `time = "sequence"` for occasion order, and
#' `profile_covariates = "previous_grade"` for a membership-covariate model.
#'
#' In the generating model, students in the `"committed"` class have a
#' stationary probability of an engaged enrolment of
#' approximately 0.8, compared with approximately 0.25 in the `"wavering"`
#' class.
#'
#' The activity indicators were simulated as continuous values intended to
#' represent log-transformed activity measures, then standardized within course
#' and rounded. Raw event counts were not generated. Each indicator has
#' approximately mean zero and standard deviation one within course; for example, `browse = 1.2` represents a value 1.2 standard
#' deviations above the course mean. The simulation does not include
#' between-course differences in activity volume.
#'
#' @section Local dependence:
#' `attendance` is generated from the four click measures, introducing
#' within-profile dependence among the indicators. A model with diagonal
#' residual covariance therefore does not reproduce the generating covariance
#' structure. Use `get_results(fit, "residuals")` to examine residual
#' associations and `covariance_model = "full"` to estimate within-profile
#' covariances. Class enumeration under a diagonal model may represent
#' residual association by additional classes. The two generating profiles
#' therefore do not imply that every fitted specification will favour two
#' profiles.
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
#'   \item{browse}{Numeric, simulated browsing activity, standardized within
#'     course.}
#'   \item{lectures}{Numeric, simulated lecture-viewing activity, standardized
#'     within course.}
#'   \item{forum_read}{Numeric, simulated forum-reading activity, standardized
#'     within course.}
#'   \item{forum_post}{Numeric, simulated forum-posting activity, standardized
#'     within course.}
#'   \item{attendance}{Numeric, simulated attendance activity, standardized
#'     within course. Generated from the four other activity indicators with
#'     added noise, introducing within-profile dependence; see
#'     *Local dependence*.}
#'   \item{previous_grade}{Numeric, simulated grade in the previous course,
#'     standardized across the cohort; an entry grade for the first course.
#'     Can be used as a membership predictor through the `profile_covariates`
#'     argument of [multilpa()] or in [r3step()].}
#'   \item{engagement}{Factor, `"disengaged"` or `"engaged"`: the pattern
#'     used to generate the observation; retained for assessing class recovery.}
#'   \item{student_type}{Factor, `"committed"` or `"wavering"`: the
#'     student's generating class; retained for assessing group-class recovery.}
#' }
#' @source Simulated by `data-raw/course-engagement.R` in the source
#'   repository, seed 2026. The
#'   simulation is calibrated to aggregate constants of a real
#'   learning-analytics export (scales, spreads, correlations and the
#'   persistence of engagement across courses). No row, identifier or value
#'   of any real student is included, and the source data are not
#'   distributed.
#' @seealso `vignette("lpa", package = "latents")` for the analysis these data are used in,
#'   [multilpa()] for the two-level model, [lta()] for the
#'   sequence, and [r3step()] for membership-covariate analysis.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' as.data.frame(fit)
#' get_results(fit, what = "profile_probabilities")
#' # Compare fitted assignments with the generating engagement profiles.
#' get_results(fit, what = "assignments", data = course_engagement,
#'          truth = "engagement")
"course_engagement"

#' Leisure activities of university students in daily life
#'
#' @description
#' `student_esm` contains experience-sampling reports of leisure activities
#' and affect from 100 university students. Each row represents a survey
#' response collected when the student reported no studying since the previous
#' prompt. Eight binary items describe activities undertaken during that
#' interval, and four ratings describe affect.
#'
#' Repeated reports are nested within students, providing an example for
#' multilevel latent class analysis of activity patterns and differences in
#' their prevalence between students. The dataset retains complete reports
#' from students with at least 15 eligible responses. The affect ratings can
#' also be included as indicators in a mixed measurement model.
#'
#' @details
#' In the original study students answered six prompts a day for fourteen
#' days; the activity items were asked only at prompts where the student had
#' not studied since the previous prompt. The dataset contains 100 students,
#' drawn with a fixed seed from those with at least 15 complete non-study
#' reports, retaining every such report: 15 to 51 per student (median 24). The eight
#' activity items are categorical indicators; the four affect ratings can be
#' added as continuous indicators in a mixed measurement model. `worried` has
#' a pronounced floor effect (49.7% of ratings are 1). A Gaussian component
#' representing low-worry prompts may therefore reach the `min_variance`
#' bound.
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
#' @encoding UTF-8
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
#' @description
#' `study_hours` is a simulated dataset illustrating how the relationship
#' between study time and quiz performance can differ across latent classes.
#' Each row represents one study week for a student, with six weeks recorded
#' for each of 150 students.
#'
#' The data were generated from two week-level regression classes, labelled
#' "deep" and "surface", with different intercepts and slopes relating study
#' hours to quiz scores. Two student-level classes differ in the prevalence
#' of these weekly classes, allowing class membership to vary across weeks
#' within students.
#'
#' The generating labels are retained for assessing class recovery. Additional
#' binary and count outcomes provide examples of mixture regression with other
#' response distributions.
#'
#' @details
#' The generating model contains two week-level study strategies with
#' different regressions of quiz `score` on `hours`. The `"deep"` strategy has
#' an intercept of approximately 35 and a slope of 4.5 points per hour; the
#' `"surface"` strategy has an intercept of approximately 55 and a slope
#' below one point per hour. Students in the `"steady"` class use the deep
#' strategy in approximately 85% of weeks, compared with 25% for the
#' `"erratic"` class. In the simulation, `motivation` predicts student-level
#' membership and `sleep` predicts week-level membership.
#'
#' `score ~ hours` is therefore a two-class mixture at the week level
#' (`class_level = "observation"`) with class shares that differ by student
#' (`id = "student"`, `n_group_classes = 2`). The binary outcome `passed`
#' depends on the student's generating class and illustrates a mixture at
#' the student level (`class_level = "group"`). `questions` is a count whose Poisson
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
#'     strategy; retained for assessing class recovery.}
#'   \item{student_type}{Factor, `"steady"` or `"erratic"`: the student's
#'     generating class; retained for assessing class recovery.}
#' }
#' @source Simulated by `data-raw/study-hours.R` in the source repository.
#' @seealso [mixture_regression()], [enumerate_regressions()].
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' get_results(fit, "coefficients")
"study_hours"

#' Self-regulated learning scores generated by a large language model
#'
#' @description
#' `srl` contains self-regulated learning scores derived from questionnaire
#' responses generated by a large language model (GPT). Each row represents
#' one of 300 synthetic respondents and contains five construct scores from
#' the Motivated Strategies for Learning Questionnaire. Scores are means of
#' the corresponding items, expressed on a scale from 1 to 7.
#'
#' The dataset illustrates latent profile analysis with correlated indicators.
#' Strong associations among four constructs provide an example for examining
#' how assumptions about within-profile covariance influence fitted solutions.
#' Test anxiety shows a different pattern of association with the remaining
#' constructs.
#'
#' @details
#' Correlations among cognitive strategies, intrinsic value, self-efficacy and
#' self-regulation range from approximately 0.78 to 0.86. Their correlations
#' with test anxiety range from approximately -0.30 to -0.10.
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

#' Simulated achievement scores with three growth trajectories
#'
#' @description
#' `growth_scores` is a simulated longitudinal dataset illustrating growth
#' mixture models. Each row represents an achievement score for one of 300
#' students, measured on up to seven occasions.
#'
#' The data were generated from three trajectory classes, labelled
#' "improving", "stable" and "declining". Within each class, students differ
#' in their initial achievement and rate of change through random intercepts
#' and slopes.
#'
#' The dataset illustrates how growth mixture models distinguish differences
#' between class mean trajectories from variation among individuals within the
#' same class. Generating class labels are retained for assessing recovery,
#' and randomly omitted observations produce an unbalanced panel.
#'
#' @details
#' Occasions are indexed by `wave`, from 0 to 6. The generating classes have
#' three linear mean
#' trajectories: *improving* (intercept 48, slope 3.5 points per occasion),
#' *stable* (intercept 62, slope 0.5), and *declining* (intercept 58, slope
#' -2.5). Within each class, individual trajectories vary according to random
#' intercepts and slopes with standard deviations 4 and 0.8 and correlation
#' -0.3. The residual standard deviation is 3. Motivation is a time-invariant
#' predictor of class membership, with higher values associated with greater
#' probability of the improving class and lower probability of the declining
#' class. Approximately 12% of test occasions were removed at random,
#' producing an unbalanced panel.
#'
#' A model without random effects may represent within-class variation by
#' additional trajectory classes.
#'
#' @format A data frame with 1844 rows (300 students) and 5 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 300. Pass it as `id`.}
#'   \item{wave}{Integer test occasion, 0 to 6. The time variable.}
#'   \item{score}{Achievement score.}
#'   \item{motivation}{Standardized motivation, constant within student.}
#'   \item{trajectory}{Factor, `"improving"`, `"stable"` or `"declining"`:
#'     the student's generating class; retained for assessing class recovery.}
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
#' @description
#' `growth_schools` is a simulated dataset illustrating multilevel growth
#' mixture models. It contains repeated reading scores for 600 students nested
#' within 40 schools, with up to five observations per student.
#'
#' The data were generated from two student trajectory classes, labelled
#' "improving" and "plateauing". Students vary in their initial reading
#' achievement and rate of change within each class. Two school-level classes
#' differ in the prevalence of these student trajectories.
#'
#' The dataset illustrates how multilevel growth mixture models represent
#' individual variation in development alongside differences in trajectory
#' membership between schools. Generating class labels are retained at both
#' levels for assessing recovery.
#'
#' @details
#' There are 15 students per school, with occasions indexed by `wave` from 0
#' to 4. The generating student classes have two linear mean
#' trajectories: *improving* (intercept 48, slope 3 points per occasion) and
#' *plateauing* (intercept 56, slope 0.3). Random intercepts and slopes have
#' standard deviations 3 and 0.6 and correlation -0.2; the residual standard
#' deviation is 2. The generating school classes, labelled *supportive* and
#' *struggling*, have improving-class probabilities of 0.8 and 0.25,
#' respectively. The standardized school-level `programme` variable predicts
#' school-class membership. Approximately 8% of test occasions were removed
#' at random.
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
#'     generating class; retained for assessing class recovery.}
#'   \item{trajectory}{Factor, `"improving"` or `"plateauing"`: the student's
#'     generating class; retained for assessing class recovery.}
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
