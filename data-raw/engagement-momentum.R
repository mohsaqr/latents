# Generates the bundled example dataset `engagement_momentum`. Run with:
#   Rscript data-raw/engagement-momentum.R
#
# The data are SIMULATED. Nothing here is a measurement of any real person.
#
# The design follows the study of engagement dynamics over a full programme by
# López-Pernas and Saqr (2024, Computers in Human Behavior, 153, 108126):
# students take a sequence of courses; in every course their online activity
# places them in one of three engagement states (active, average,
# disengaged); and over the programme they follow one of three trajectories
# (engaged, fluctuating, mostly disengaged). Engaged students show momentum:
# an active course tends to be followed by another, and achievement builds on
# the previous course. Mostly disengaged students show inertia: once
# disengaged they tend to stay so. The parameter values below are
# illustrative choices that reproduce these qualitative patterns; they are not
# estimates from that study.
#
# Unlike `course_engagement`, whose attendance measure is built from the click
# measures, the five indicators here are conditionally independent within a
# state, so a diagonal covariance structure is the correct model. The
# trajectories differ both in the mix of states they visit (which a multilevel
# latent profile model recovers) and in their transitions (which a latent
# transition model recovers).

set.seed(2024)

n_students <- 106L
n_courses <- 32L
indicators <- c("browse", "lectures", "forum_read", "forum_post",
                "active_days")
states <- c("active", "average", "disengaged")
trajectories <- c("engaged", "fluctuating", "mostly_disengaged")

# ---- measurement ----------------------------------------------------------
# State means and standard deviations on a log-activity scale, drawn
# independently for every indicator, then standardized within course below.
state_means <- rbind(
  active     = c(3.90, 3.60, 3.90, 3.20, 3.40),
  average    = c(3.30, 3.05, 3.15, 2.60, 2.85),
  disengaged = c(2.70, 2.55, 2.40, 2.05, 2.30))
state_sds <- rbind(
  active     = c(0.50, 0.55, 0.52, 0.55, 0.42),
  average    = c(0.45, 0.50, 0.48, 0.50, 0.40),
  disengaged = c(0.42, 0.45, 0.45, 0.45, 0.38))

# ---- trajectories ---------------------------------------------------------
trajectory_share <- c(engaged = 0.35, fluctuating = 0.40,
                      mostly_disengaged = 0.25)
initial <- rbind(
  engaged           = c(0.60, 0.30, 0.10),
  fluctuating       = c(0.30, 0.40, 0.30),
  mostly_disengaged = c(0.10, 0.30, 0.60))
# Rows: the current state (active, average, disengaged); columns: the next.
transition <- list(
  # Momentum: activity persists, and the other states return to it.
  engaged = rbind(c(0.85, 0.12, 0.03),
                  c(0.45, 0.45, 0.10),
                  c(0.30, 0.40, 0.30)),
  # Frequent switching between all three states.
  fluctuating = rbind(c(0.35, 0.40, 0.25),
                      c(0.30, 0.40, 0.30),
                      c(0.25, 0.40, 0.35)),
  # Inertia: disengagement persists.
  mostly_disengaged = rbind(c(0.20, 0.30, 0.50),
                            c(0.10, 0.35, 0.55),
                            c(0.03, 0.12, 0.85)))

trajectory <- sample(trajectories, n_students, TRUE, trajectory_share)
sequence_length <- sample(10:15, n_students, replace = TRUE)

rows <- lapply(seq_len(n_students), function(student) {
  kind <- trajectory[student]
  n <- sequence_length[student]
  state <- integer(n)
  state[1L] <- sample.int(3L, 1L, prob = initial[kind, ])
  # A first-order chain over the student's courses: each state is drawn from
  # the row of the previous one, so the recursion is written as a loop.
  for (occasion in seq_len(n)[-1L]) {
    state[occasion] <- sample.int(3L, 1L,
                                  prob = transition[[kind]][state[occasion - 1L], ])
  }
  data.frame(student = student, sequence = seq_len(n),
             course = sample.int(n_courses, n, replace = FALSE),
             trajectory = kind, state = state, stringsAsFactors = FALSE)
})
engagement_momentum <- do.call(rbind, rows)

draws <- vapply(seq_len(nrow(engagement_momentum)), function(i) {
  k <- engagement_momentum$state[i]
  stats::rnorm(length(indicators), state_means[k, ], state_sds[k, ])
}, numeric(length(indicators)))
activity <- as.data.frame(t(draws))
names(activity) <- indicators

# ---- per-course standardization -------------------------------------------
# Each indicator is expressed in standard deviations from its own course's
# mean, as in `course_engagement`.
standardize_within <- function(values, by) {
  ave(values, by, FUN = function(z) {
    spread <- stats::sd(z)
    if (!is.finite(spread) || spread <= 0) return(z - mean(z))
    (z - mean(z)) / spread
  })
}
activity[] <- lapply(activity, standardize_within, by = engagement_momentum$course)
engagement_momentum <- cbind(engagement_momentum,
                             as.data.frame(lapply(activity, round, 2)))

# ---- achievement ----------------------------------------------------------
# The grade of each course depends on the engagement state in that course and,
# through momentum, on the grade of the previous course. `previous_grade` is
# the grade of the course before; at a student's first course an entry grade
# stands in its place.
state_effect <- c(active = 0.60, average = 0.00, disengaged = -0.70)
grade <- numeric(nrow(engagement_momentum))
previous <- numeric(nrow(engagement_momentum))
by_student <- split(seq_len(nrow(engagement_momentum)), engagement_momentum$student)
invisible(lapply(by_student, function(index) {
  last <- stats::rnorm(1L)
  # Momentum is a recursion over the student's courses, so a loop again.
  for (i in index) {
    previous[i] <<- last
    grade[i] <<- state_effect[[engagement_momentum$state[i]]] + 0.35 * last +
      stats::rnorm(1L, 0, 0.75)
    last <- grade[i]
  }
}))
grade_centre <- mean(grade)
grade_scale <- stats::sd(grade)
engagement_momentum$grade <- round((grade - grade_centre) / grade_scale, 2)
engagement_momentum$previous_grade <- round((previous - mean(previous)) /
                                              stats::sd(previous), 2)

engagement_momentum$state <- factor(states[engagement_momentum$state],
                                    levels = states)
engagement_momentum$trajectory <- factor(engagement_momentum$trajectory,
                                         levels = trajectories)
engagement_momentum <- engagement_momentum[c("student", "course", "sequence",
                                             indicators, "grade",
                                             "previous_grade", "state",
                                             "trajectory")]
row.names(engagement_momentum) <- NULL

stopifnot(
  "every student must have a clean 1..n sequence" =
    all(vapply(split(engagement_momentum$sequence, engagement_momentum$student),
               function(s) identical(s, seq_along(s)), logical(1))),
  "no missing values" = !anyNA(engagement_momentum),
  "each indicator must be standardized within every course" =
    all(vapply(indicators, function(v) {
      by_course <- split(engagement_momentum[[v]], engagement_momentum$course)
      all(vapply(by_course, function(z) {
        abs(mean(z)) < 0.02 && (length(z) < 2L || abs(stats::sd(z) - 1) < 0.02)
      }, logical(1)))
    }, logical(1))))

if (!identical(Sys.getenv("ENGAGEMENT_MOMENTUM_DRY_RUN"), "true")) {
  usethis::use_data(engagement_momentum, overwrite = TRUE)
}
