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
#' finds it: [bivariate_residuals()] reports a largest residual of 0.21 with
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
#' The activity indicators are natural-log counts, `log1p(events)`. Raw
#' learning-analytics counts are strongly right-skewed, and on an unlogged draft
#' of these data the skew alone looked like an extra class, exactly as a floor
#' or ceiling effect does. Back-transform with `expm1()` to read them as events:
#' a `browse` of 3.9 is about 48 page views.
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
#'   \item{browse}{Numeric, `log1p` of course page views.}
#'   \item{lectures}{Numeric, `log1p` of lecture videos viewed.}
#'   \item{forum_read}{Numeric, `log1p` of forum posts read.}
#'   \item{forum_post}{Numeric, `log1p` of forum posts written.}
#'   \item{attendance}{Numeric, `log1p` of days the student was active in the
#'     course. It is **generated from the four click measures** and then
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
#'     to [fit_covariates()], or as the predictor to [r3step()].}
#'   \item{engagement}{Factor, `"disengaged"` or `"engaged"`, the pattern the
#'     row was generated from. This is the truth behind the data, not something
#'     a study would observe; it is shipped so a fitted model can be checked
#'     against what produced it. It is deliberately not called `profile`:
#'     [assignments()] adds a column of that name, and a dataset already
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
#' @seealso `vignette("multilpa")` for the analysis these data are used in,
#'   [multilpa()] for the two-level model, [fit_transitions()] for the
#'   sequence, and [fit_covariates()] for `previous_grade`.
#' @examples
#' activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
#' fit <- multilpa(course_engagement, vars = activity, id = "student",
#'                 n_profiles = 2, n_group_classes = 2, n_starts = 4, seed = 1)
#' as.data.frame(fit)
#' as.data.frame(fit, what = "profile_probabilities")
#' # The generating truth ships alongside, so recovery can be checked directly.
#' xtabs(~ profile + engagement, data = assignments(fit, data = course_engagement))
"course_engagement"
