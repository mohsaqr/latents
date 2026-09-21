# Generates the bundled example dataset `course_engagement`. Run with:
#   Rscript data-raw/course-engagement.R
#
# The data are SIMULATED. Nothing here is a measurement of any real person.
#
# The simulation is calibrated to the shape of a real learning-analytics export
# -- students taking a sequence of courses, with activity counts recorded per
# enrolment -- so that the scales, the spread and the correlations a reader
# meets are the ones such data actually has. Only aggregate constants from that
# calibration appear below: the marginal means and standard deviations on the
# log scale, the separation between an engaged and a disengaged pattern, and how
# persistent that pattern is across a student's courses. No row, identifier or
# value of the source survives, and the source is not distributed.
#
# Because the data are simulated, the truth behind them ships beside the
# observations: `profile` is the generating enrolment pattern and
# `student_type` the generating student kind. A reader can therefore check
# whether the model recovers what produced the data, which is the one thing real
# data can never be used to check.

set.seed(2026)

n_students <- 106L
n_courses <- 32L
max_sequence <- 15L

# ---- measurement ----------------------------------------------------------
# Activity counts are strongly right-skewed, so both the calibration and the
# model live on the log scale; the shipped indicators are `log1p(count)`.
# An earlier draft shipped the raw counts and the skew alone looked like an
# extra class, exactly as a floor effect does.
indicators <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
#
# The real export's disengaged pattern piles up against zero on the two forum
# measures -- most disengaged enrolments never post at all. That pile-up is real
# but it is a floor, and a floor is indistinguishable from an extra class: an
# earlier draft reproduced it faithfully and the criteria then chose three
# profiles out of a simulated two. The disengaged means below are therefore
# lifted clear of zero, keeping the separation the real data shows without the
# artefact. Every mean sits more than three standard deviations above zero, so
# the indicators stay non-negative, as a log count must be, without clamping.
#
# The separation is deliberately smaller than the raw export's. Fitted to the
# real data the two patterns barely overlap, and a bundled example that
# classifies every case with probability one teaches nothing about
# classification uncertainty -- the entropy is 1.000 and every diagnostic in the
# package returns the same answer. These means keep the ordering and the shape
# of the real contrast at a separation a reader can actually see working.
profile_means <- rbind(
  disengaged = c(browse = 2.90, lectures = 2.70, forum_read = 2.40,
                 forum_post = 2.10, attendance = 2.30),
  engaged    = c(browse = 3.90, lectures = 3.60, forum_read = 3.90,
                 forum_post = 3.20, attendance = 3.00))
# The disengaged spreads are tighter than the engaged ones, which is what the
# real data shows too: there are many ways to be busy in a course and rather few
# ways to be absent from it. Each keeps its mean at least four standard
# deviations above zero.
profile_sds <- rbind(
  disengaged = c(0.60, 0.62, 0.58, 0.55, 0.45),
  engaged    = c(0.63, 0.75, 0.70, 0.72, 0.43))

# ---- structure ------------------------------------------------------------
# Two kinds of student. Both kinds contain both patterns; they differ in the
# MIX, which is exactly the two-level quantity this package estimates. An
# earlier draft separated the kinds by how PERSISTENT the pattern was instead,
# and the group classes then recovered nothing: persistence is a property of the
# sequence, not of the mix, so a model of the mix cannot see it.
type_share <- c(committed = 0.55, wavering = 0.45)
engaged_share <- c(committed = 0.80, wavering = 0.25)
# Persistence from one course to the next. Both kinds are sticky -- engagement
# is a habit -- and the stationary distribution of each chain is its kind's
# engaged share above, so the sequence and the mix tell the same story.
transition <- list(
  committed = rbind(c(0.60, 0.40), c(0.10, 0.90)),
  wavering  = rbind(c(0.90, 0.10), c(0.30, 0.70)))

student_type <- sample(names(type_share), n_students, TRUE, type_share)
# A few students take only part of the programme, so the sequences are ragged.
sequence_length <- pmin(max_sequence,
                        pmax(3L, stats::rbinom(n_students, max_sequence, 0.88)))

rows <- lapply(seq_len(n_students), function(student) {
  kind <- student_type[student]
  n <- sequence_length[student]
  state <- integer(n)
  state[1L] <- 1L + (stats::runif(1) < engaged_share[[kind]])
  # A first-order chain over the student's courses. Vectorising this would need
  # the whole path before its own transitions are drawn, so the loop is the
  # honest expression of the recursion.
  for (occasion in seq_len(n)[-1L]) {
    state[occasion] <- sample.int(2L, 1L,
                                  prob = transition[[kind]][state[occasion - 1L], ])
  }
  data.frame(student = student, sequence = seq_len(n),
             course = sample.int(n_courses, n, replace = FALSE),
             student_type = kind, engagement = state, stringsAsFactors = FALSE)
})
course_engagement <- do.call(rbind, rows)

draws <- vapply(seq_len(nrow(course_engagement)), function(i) {
  k <- course_engagement$engagement[i]
  stats::rnorm(length(indicators), profile_means[k, ], profile_sds[k, ])
}, numeric(length(indicators)))
activity <- as.data.frame(t(draws))
names(activity) <- indicators

# `attendance` is the count of days the student was active, so it is not an
# independent measurement of engagement: a student is recorded present on a day
# BECAUSE they clicked something. It is therefore built from the click measures
# rather than drawn beside them, and then blurred, because a day is a coarse
# unit -- one tick however much happened inside it.
#
# This makes the indicators locally dependent, which the real export shows too
# (Active_days correlates 0.68 to 0.92 with the click counts, far above what a
# shared profile alone would produce). It is deliberate: it gives
# `bivariate_residuals()` a real violation to find and `covariance_model =
# "full"` something real to model, on data whose answer is known. A dataset
# where local independence holds exactly can only ever show those verbs
# agreeing that nothing is wrong.
clicks <- rowMeans(activity[c("browse", "lectures", "forum_read", "forum_post")])
activity$attendance <- activity$attendance +
  0.62 * (clicks - mean(clicks)) + stats::rnorm(nrow(activity), 0, 0.30)
# ---- per-course standardization -------------------------------------------
# The shipped indicators are standardized within course: each enrolment's value
# is expressed in standard deviations from its own course's mean. This is the
# scale a learning-analytics reader meets, because activity volume is not
# comparable across courses -- a course with weekly quizzes generates more
# clicks than a seminar whatever its students are like -- and an analyst
# therefore compares a student to their coursemates rather than to the cohort.
#
# It is applied to the log values, not to counts: standardizing is affine
# within each course, so it leaves the within-course skew exactly as it found
# it, and the skew is what made an earlier draft's criteria read a floor as an
# extra class. Across the pooled data the map is piecewise affine, one shift
# and scale per course, so pooled skew does move -- from 3.35 on counts to 1.96
# -- but not nearly far enough on its own. The log does that work; this puts
# the courses on a common footing afterwards.
standardize_within <- function(values, by) {
  ave(values, by, FUN = function(z) {
    spread <- stats::sd(z)
    # A course with one enrolment, or no variation in it, has nothing to divide
    # by; centring still applies and the scale is left alone.
    if (!is.finite(spread) || spread <= 0) return(z - mean(z))
    (z - mean(z)) / spread
  })
}
activity[] <- lapply(activity, standardize_within, by = course_engagement$course)
course_engagement <- cbind(course_engagement,
                           as.data.frame(lapply(activity, round, 2)))

# ---- previous_grade -------------------------------------------------------
# The grade the student earned in the course BEFORE this one, standardised
# across the cohort. At a student's first course there is no previous course,
# so an entry grade stands in its place. It is a predictor of engagement here,
# not an outcome of it: a student who did well last time is more likely to be
# engaged this time, which is what `multilpa()` and `r3step()` model.
grade_now <- stats::rnorm(nrow(course_engagement),
                          ifelse(course_engagement$engagement == 2L, 0.45, -0.45), 0.85)
first <- !duplicated(course_engagement$student)
previous <- c(NA_real_, grade_now[-length(grade_now)])
previous[first] <- stats::rnorm(sum(first), 0, 1)
course_engagement$previous_grade <- round(as.numeric(scale(previous)), 2)

# The truth column is `engagement`, not `profile`: `assignments()` adds a column
# called `profile`, and a bundled dataset that already carries that name makes
# the package refuse to join its own output onto it. A reader checking recovery
# would have had to drop the column first, which is exactly the ritual the
# package exists to remove.
course_engagement$engagement <- factor(
  c("disengaged", "engaged")[course_engagement$engagement],
  levels = c("disengaged", "engaged"))
course_engagement$student_type <- factor(course_engagement$student_type,
                                         levels = c("committed", "wavering"))
course_engagement <- course_engagement[c("student", "course", "sequence",
                                         indicators, "previous_grade",
                                         "engagement", "student_type")]
row.names(course_engagement) <- NULL

stopifnot(
  "every student must have a clean 1..n sequence" =
    all(vapply(split(course_engagement$sequence, course_engagement$student),
               function(s) identical(s, seq_along(s)), logical(1))),
  "no missing values" = !anyNA(course_engagement),
  "indicators must be finite" =
    all(vapply(course_engagement[indicators], function(v)
      all(is.finite(v)), logical(1))),
  # The standardization is asserted rather than assumed: within every course
  # each indicator must centre on zero and scale to one. Rounding to two
  # decimals is what the tolerance allows for, and the assertion is what says
  # the shipped scale is the one the documentation claims.
  "each indicator must be standardized within every course" =
    all(vapply(indicators, function(v) {
      by_course <- split(course_engagement[[v]], course_engagement$course)
      all(vapply(by_course, function(z) {
        abs(mean(z)) < 0.02 && (length(z) < 2L || abs(stats::sd(z) - 1) < 0.02)
      }, logical(1)))
    }, logical(1))))

usethis::use_data(course_engagement, overwrite = TRUE)
