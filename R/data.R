#' Student engagement in sixty schools
#'
#' A simulated cross-sectional dataset with the structure this package exists
#' for: students nested in schools, where the schools differ in the *mix* of
#' students they serve rather than in the kinds of student that exist.
#'
#' Half the schools were generated drawing three quarters of their students
#' from an engaged kind, half drawing one quarter. Pooled across the sample
#' 48% of students are engaged, a headline number that describes no school in
#' the data: the per-school share runs from 0.08 to 1.00. Recovering that split
#' is what a two-level model does and a pooled one cannot.
#'
#' The indicators are on wide, unbounded scales on purpose. Bounding them to a
#' 1--7 response scale piles mass against the ceiling, and that looks exactly
#' like an extra class: on a bounded draft of these data the criteria
#' confidently chose three profiles where two were simulated.
#'
#' @format A data frame with 720 rows (60 schools x 12 students) and 6 columns:
#' \describe{
#'   \item{school}{Integer school identifier, 1 to 60. The nesting unit; pass it
#'     as `id`.}
#'   \item{term}{Integer 1 to 12, the student's position within the school. Not
#'     a repeated measurement -- each row is a different student.}
#'   \item{engaged}{Logical. The kind the row was generated from. This is the
#'     truth behind the data, not something a study would observe; it is
#'     shipped so a fitted model can be checked against what produced it.}
#'   \item{homework_hours}{Numeric, hours of homework per week.}
#'   \item{participation}{Numeric, a participation score.}
#'   \item{interest}{Numeric, a self-reported interest score.}
#' }
#' @source Simulated by `data-raw/engagement.R`, seed 2024.
#' @seealso [engagement_panel] for the repeated-measures counterpart,
#'   `vignette("multilpa")` for the analysis these data are used in.
#' @examples
#' fit <- multilpa(school_engagement,
#'                 vars = c("homework_hours", "participation", "interest"),
#'                 id = "school", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 4, seed = 1)
#' as.data.frame(fit)
#' as.data.frame(fit, what = "profile_probabilities")
"school_engagement"

#' Student engagement over four waves
#'
#' A simulated panel with the structure a latent transition model is for: the
#' same students measured repeatedly, where the question is not what kinds
#' there are but who moves between them.
#'
#' Students were generated to stay in their current state from one wave to the
#' next with probability 0.85. Four waves of 120 students is enough to recover
#' that persistence to within a few points, which is the honest amount of
#' evidence a panel this size carries -- close, not exact.
#'
#' @format A data frame with 480 rows (120 students x 4 waves) and 5 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 120. Pass it as `id`.}
#'   \item{wave}{Integer 1 to 4, the occasion. Pass it as `time`.}
#'   \item{state}{Integer 1 or 2, the state the row was generated from. The
#'     truth behind the data, shipped so the fitted transition matrix can be
#'     compared with the one that produced it.}
#'   \item{homework_hours}{Numeric, hours of homework per week.}
#'   \item{participation}{Numeric, a participation score.}
#' }
#' @source Simulated by `data-raw/engagement.R`, seed 11.
#' @seealso [school_engagement] for the cross-sectional counterpart,
#'   [fit_transitions()] for the model these data are for.
#' @examples
#' moves <- fit_transitions(engagement_panel,
#'                          vars = c("homework_hours", "participation"),
#'                          id = "student", time = "wave", n_profiles = 2,
#'                          n_starts = 4, seed = 1)
#' transitions(moves, stable = TRUE)
"engagement_panel"
