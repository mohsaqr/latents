# Generates the bundled simulated dataset `growth_scores`. Run from the package
# root with:
#   Rscript data-raw/growth-scores.R
#
# 300 students sit an achievement test on up to seven occasions (waves 0 to 6).
# Three kinds of student follow three trajectories: improving (a low start and
# a steep rise), stable (a high start, nearly flat) and declining (a middling
# start and a fall). Within each kind, students vary around the trajectory:
# each has their own starting level and rate of change (a random intercept
# and slope, correlated -0.3, the same spread in every kind), so a model
# without random effects needs extra classes to absorb that spread. Motivation,
# measured once, makes a student more likely to be improving and less likely
# to be declining. About 12% of tests were missed, at random.
set.seed(20261002)
n_students <- 300L
waves <- 0:6
motivation <- round(stats::rnorm(n_students), 2)
logits <- cbind(improving = 0.4 + 0.9 * motivation, stable = 0.2,
                declining = 0 - 0.9 * motivation)
probabilities <- exp(logits) / rowSums(exp(logits))
kind <- apply(probabilities, 1L, function(p) sample(colnames(logits), 1L, prob = p))
trajectory <- factor(kind, levels = c("improving", "stable", "declining"))
start <- c(improving = 48, stable = 62, declining = 58)[kind]
slope <- c(improving = 3.5, stable = 0.5, declining = -2.5)[kind]
covariance <- matrix(c(16, -0.3 * 4 * 0.8, -0.3 * 4 * 0.8, 0.64), 2L)
effects <- matrix(stats::rnorm(2L * n_students), n_students) %*% chol(covariance)
growth_scores <- data.frame(
  student = rep(seq_len(n_students), each = length(waves)),
  wave = rep(waves, n_students),
  motivation = rep(motivation, each = length(waves)),
  trajectory = rep(trajectory, each = length(waves)))
person <- growth_scores$student
growth_scores$score <- round(start[person] + effects[person, 1L] +
                               (slope[person] + effects[person, 2L]) * growth_scores$wave +
                               stats::rnorm(nrow(growth_scores), 0, 3), 1)
growth_scores <- growth_scores[stats::runif(nrow(growth_scores)) > 0.12,
                               c("student", "wave", "score", "motivation", "trajectory")]
row.names(growth_scores) <- NULL
stopifnot(!anyNA(growth_scores), length(unique(growth_scores$student)) == n_students)
save(growth_scores, file = file.path("data", "growth_scores.rda"), compress = "xz")
