# Generates the bundled simulated dataset `growth_schools`. Run from the
# package root with:
#   Rscript data-raw/growth-schools.R
#
# 600 students in 40 schools sit a reading test on five occasions (waves 0
# to 4). Two kinds of student follow two trajectories: improving (a low start
# and a steady rise) and plateauing (a higher start, nearly flat). Students
# vary around their trajectory (a random intercept and slope). Schools differ
# in how many students improve: in "supportive" schools most do, in
# "struggling" schools few do, and a school's support programme, measured
# once, makes it more likely to be supportive. About 8% of tests were
# missed, at random.
set.seed(20261003)
n_schools <- 40L
per_school <- 15L
waves <- 0:4
programme <- round(stats::rnorm(n_schools), 2)
supportive <- stats::runif(n_schools) < stats::plogis(0.8 * programme)
school_type <- factor(ifelse(supportive, "supportive", "struggling"),
                      levels = c("supportive", "struggling"))
improving_share <- ifelse(supportive, 0.8, 0.25)
n_students <- n_schools * per_school
school <- rep(seq_len(n_schools), each = per_school)
improving <- stats::runif(n_students) < improving_share[school]
trajectory <- factor(ifelse(improving, "improving", "plateauing"),
                     levels = c("improving", "plateauing"))
start <- ifelse(improving, 48, 56)
slope <- ifelse(improving, 3, 0.3)
covariance <- matrix(c(9, -0.2 * 3 * 0.6, -0.2 * 3 * 0.6, 0.36), 2L)
effects <- matrix(stats::rnorm(2L * n_students), n_students) %*% chol(covariance)
growth_schools <- data.frame(
  school = rep(school, each = length(waves)),
  student = rep(seq_len(n_students), each = length(waves)),
  wave = rep(waves, n_students),
  programme = rep(programme[school], each = length(waves)),
  school_type = rep(school_type[school], each = length(waves)),
  trajectory = rep(trajectory, each = length(waves)))
person <- growth_schools$student
growth_schools$score <- round(start[person] + effects[person, 1L] +
                                (slope[person] + effects[person, 2L]) *
                                growth_schools$wave +
                                stats::rnorm(nrow(growth_schools), 0, 2), 1)
growth_schools <- growth_schools[stats::runif(nrow(growth_schools)) > 0.08,
                                 c("school", "student", "wave", "score", "programme",
                                   "school_type", "trajectory")]
row.names(growth_schools) <- NULL
stopifnot(!anyNA(growth_schools),
          length(unique(growth_schools$student)) == n_students,
          length(unique(growth_schools$school)) == n_schools)
save(growth_schools, file = file.path("data", "growth_schools.rda"), compress = "xz")
