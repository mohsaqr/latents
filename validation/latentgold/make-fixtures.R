# Build the datasets for the Latent GOLD comparison, and record the values this
# package produces for them. Run from the project root:
#   Rscript validation/latentgold/make-fixtures.R
#
# The targets are the three things with no external anchor: bivariate residuals
# for continuous indicators, the BCH distal-outcome correction, and the
# covariate correction. Latent GOLD is the reference implementation for all
# three; its co-author wrote the three-step papers.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("validation", "latentgold")

## ---- dataset one: two-level Gaussian, with a planted local dependence -------
set.seed(4)
n_groups <- 40L; per <- 10L
group <- rep(seq_len(n_groups), each = per)
n <- length(group)
group_class <- rep(rep(c(1L, 2L), length.out = n_groups), each = per)
profile <- 1L + as.integer(stats::runif(n) < ifelse(group_class == 2L, 0.8, 0.2))
shared <- stats::rnorm(n)                      # y1 and y2 share this; y3 does not
bvr_data <- data.frame(
  id = seq_len(n), group = group,
  y1 = round(stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared, 6),
  y2 = round(stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared, 6),
  y3 = round(stats::rnorm(n, ifelse(profile == 2L, 1, -1), 0.8), 6))
bvr_groups <- n_groups
utils::write.table(bvr_data, file.path(directory, "bvr.dat"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

fit <- multilpa(bvr_data, c("y1", "y2", "y3"), "group", n_profiles = 2,
                n_group_classes = 2, n_starts = 20, seed = 1)
residuals <- get_data(fit, "residuals", data = bvr_data)

## ---- dataset two: a distal outcome and a covariate --------------------------
set.seed(21)
n_groups <- 50L; per <- 12L
group <- rep(seq_len(n_groups), each = per)
n <- length(group)
x <- round(stats::rnorm(n), 6)
truth <- 1L + as.integer(stats::runif(n) < stats::plogis(-0.3 + 1.2 * x))
step_data <- data.frame(
  id = seq_len(n), group = group, x = x,
  y1 = round(stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)), 6),
  y2 = round(stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)), 6),
  distal = round(stats::rnorm(n, ifelse(truth == 2L, 10, 0)), 6))
step_groups <- n_groups
utils::write.table(step_data, file.path(directory, "threestep.dat"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

step_fit <- multilpa(step_data, c("y1", "y2"), "group", n_profiles = 2,
                     n_group_classes = 1, n_starts = 20, seed = 1)

targets <- list(
  bvr = list(
    log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters,
    means = fit$means, profile_probabilities = fit$profile_probabilities,
    residuals = residuals),
  three_step = list(
    log_likelihood = step_fit$log_likelihood,
    errors = get_data(step_fit, "classification_errors", level = "individuals"),
    distal_bch = three_step(step_fit, step_data, "distal", method = "bch"),
    distal_modal = three_step(step_fit, step_data, "distal", method = "modal"),
    covariate = r3step(step_fit, step_data, "x")))
saveRDS(targets, file.path(directory, "multilpa-targets.rds"))

cat("=== dataset one: bvr.dat ===\n")
cat("rows", nrow(bvr_data), "| groups", bvr_groups,
    "| log likelihood", format(fit$log_likelihood, digits = 10),
    "| parameters", fit$n_parameters, "\n")
cat("worst bivariate residual:\n")
print(head(residuals[c("profile", "indicator_1", "indicator_2", "observed",
                       "expected", "residual")], 3), digits = 4, row.names = FALSE)

cat("\n=== dataset two: threestep.dat ===\n")
cat("rows", nrow(step_data), "| groups", step_groups,
    "| log likelihood", format(step_fit$log_likelihood, digits = 10), "\n")
cat("classification errors:\n")
print(targets$three_step$errors, digits = 5, row.names = FALSE)
cat("\ndistal outcome, BCH vs modal:\n")
print(data.frame(class = targets$three_step$distal_bch$class,
                 bch = targets$three_step$distal_bch$estimate,
                 modal = targets$three_step$distal_modal$estimate),
      digits = 5, row.names = FALSE)
cat("\ncovariate on membership (R3STEP):\n")
print(targets$three_step$covariate[c("outcome", "term", "estimate",
                                     "standard_error")],
      digits = 5, row.names = FALSE)
