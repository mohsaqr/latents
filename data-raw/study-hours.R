# Generates the bundled simulated dataset `study_hours`. Run from the package
# root with:
#   Rscript data-raw/study-hours.R
#
# 150 students each report six study weeks. Every week a student studies with
# one of two strategies, and the strategy decides how study hours turn into a
# quiz score: "deep" weeks gain about 4.5 points per hour from a low base,
# "surface" weeks gain under one point per hour from a higher base. Students
# are of two kinds that differ in the *mix* of strategies -- a steady student
# studies deeply in about 85% of weeks, an erratic one in about 25% -- and
# motivation makes a student more likely to be steady. Sleep the night before
# the quiz makes a deep week more likely. The weekly pass/fail outcome and the
# count of questions asked follow the student's kind and the week's strategy
# respectively, so every nesting of mixture_regression() has a model that generated data.
set.seed(20260927)
n_students <- 150L
n_weeks <- 6L
student <- rep(seq_len(n_students), each = n_weeks)
week <- rep(seq_len(n_weeks), n_students)
motivation <- round(stats::rnorm(n_students), 2)
steady <- stats::runif(n_students) < stats::plogis(0.2 + 1.2 * motivation)
student_type <- factor(ifelse(steady, "steady", "erratic"),
                       levels = c("steady", "erratic"))
hours <- round(stats::runif(n_students * n_weeks, 0, 12), 1)
sleep <- round(stats::rnorm(n_students * n_weeks, 7, 1), 1)
base_logit <- ifelse(steady[student], stats::qlogis(0.85), stats::qlogis(0.25))
deep <- stats::runif(n_students * n_weeks) <
  stats::plogis(base_logit + 0.5 * (sleep - 7))
strategy <- factor(ifelse(deep, "deep", "surface"),
                   levels = c("deep", "surface"))
score <- round(ifelse(deep, 35 + 4.5 * hours, 55 + 0.8 * hours) +
                 stats::rnorm(n_students * n_weeks, 0, ifelse(deep, 5, 7)), 1)
passed <- as.integer(stats::runif(n_students * n_weeks) <
                       ifelse(steady[student], stats::plogis(-3 + 0.6 * hours),
                              stats::plogis(0.2 - 0.05 * hours)))
questions <- stats::rpois(n_students * n_weeks,
                          ifelse(deep, exp(0.3 + 0.12 * hours),
                                 exp(1.2 - 0.05 * hours)))
study_hours <- data.frame(
  student = student, week = week, hours = hours, sleep = sleep,
  motivation = motivation[student], score = score, passed = passed,
  questions = questions, strategy = strategy,
  student_type = student_type[student])
stopifnot(nrow(study_hours) == n_students * n_weeks,
          !anyNA(study_hours))
save(study_hours, file = file.path("data", "study_hours.rda"),
     compress = "xz", version = 2)
