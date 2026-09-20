# Generates the two bundled example datasets. Run with:
#   Rscript data-raw/engagement.R
#
# Both are simulated, so the truth behind them is known and can be shipped
# beside the observations: `engaged` in school_engagement and `state` in
# engagement_panel are the generating labels, not measurements. A reader can
# therefore check whether the model recovers what produced the data, which is
# the one thing real data can never be used to check.

set.seed(2024)

n_schools <- 60L
per_school <- 12L

# Half the schools draw three quarters of their students from the engaged kind,
# half draw one quarter. This is the two-level structure the package estimates:
# the schools differ in their MIX, not in the kinds of student that exist.
engaged_share <- ifelse(rep(1:2, length.out = n_schools) == 1L, 0.75, 0.25)

school_engagement <- data.frame(
  school = rep(seq_len(n_schools), each = per_school),
  term   = rep(seq_len(per_school), times = n_schools)
)
school_engagement$engaged <-
  stats::runif(nrow(school_engagement)) < rep(engaged_share, each = per_school)

engaged <- school_engagement$engaged
school_engagement$homework_hours <-
  round(stats::rnorm(nrow(school_engagement), ifelse(engaged, 8.0, 4.2), 1.2), 1)
school_engagement$participation <-
  round(stats::rnorm(nrow(school_engagement), ifelse(engaged, 7.2, 4.0), 1.1), 1)
school_engagement$interest <-
  round(stats::rnorm(nrow(school_engagement), ifelse(engaged, 7.4, 4.3), 1.2), 1)

# The indicators are deliberately on wide, unbounded scales. An earlier draft
# bounded them to 1-7; the mass piling against the ceiling looked exactly like
# an extra class, and the criteria confidently chose three profiles out of a
# simulated two. A floor or ceiling effect is indistinguishable from a class.

set.seed(11)

n_students <- 120L
n_waves <- 4L
persistence <- 0.85

engagement_panel <- data.frame(
  student = rep(seq_len(n_students), each = n_waves),
  wave    = rep(seq_len(n_waves), times = n_students)
)

# Each wave depends on the one before, so the chain is accumulated rather than
# looped: Reduce(accumulate = TRUE) is the same idiom the package's own
# forward-backward recursion uses.
waves <- Reduce(function(previous, ignored) {
  ifelse(stats::runif(n_students) < persistence, previous, 3L - previous)
}, seq_len(n_waves - 1L),
   init = sample(1:2, n_students, replace = TRUE), accumulate = TRUE)

engagement_panel$state <- as.vector(t(do.call(cbind, waves)))
high <- engagement_panel$state == 1L
engagement_panel$homework_hours <-
  round(stats::rnorm(nrow(engagement_panel), ifelse(high, 8.0, 4.2), 1.2), 1)
engagement_panel$participation <-
  round(stats::rnorm(nrow(engagement_panel), ifelse(high, 7.2, 4.0), 1.1), 1)

stopifnot(
  "school_engagement has one row per student-term" =
    nrow(school_engagement) == n_schools * per_school,
  "engagement_panel has one row per student-wave" =
    nrow(engagement_panel) == n_students * n_waves,
  "no missing values" =
    !anyNA(school_engagement) && !anyNA(engagement_panel),
  "the schools differ in their mix, which is the point of the data" =
    diff(range(tapply(school_engagement$engaged, school_engagement$school, mean))) > 0.5)

save(school_engagement, file = "data/school_engagement.rda", compress = "xz")
save(engagement_panel, file = "data/engagement_panel.rda", compress = "xz")

cat("school_engagement:", nrow(school_engagement), "rows,",
    length(unique(school_engagement$school)), "schools\n")
cat("engagement_panel :", nrow(engagement_panel), "rows,",
    length(unique(engagement_panel$student)), "students\n")
cat("pooled engaged share:", round(mean(school_engagement$engaged), 3), "\n")
cat("per-school range    :",
    paste(round(range(tapply(school_engagement$engaged,
                             school_engagement$school, mean)), 2), collapse = " to "), "\n")
