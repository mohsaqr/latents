# Recovery against the generating parameters (Equivalence 2.0, E2.5).
#
# The Monte Carlo cases of equivalence/montecarlo/ were generated with known
# parameters (simulab stores them with each dataset). Agreement with Mplus and
# Latent GOLD shows that the package computes what they compute; this checks
# what all three might share: whether the estimates recover the parameters
# that generated the data, and how well the classes recover the true
# memberships. Writes validation/MONTECARLO_RECOVERY.md.
#
# Run from the project root:
#   Rscript validation/montecarlo-recovery.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
cases <- readRDS(file.path("equivalence", "montecarlo", "kit", "cases.rds"))
source(file.path("equivalence", "helper-montecarlo.R"))
cores <- max(1L, min(8L, parallel::detectCores() - 1L))

permutations <- function(k) if (k == 1L) list(1L) else
  do.call(c, lapply(seq_len(k), \(i) lapply(permutations(k - 1L),
                                           \(p) c(i, setdiff(seq_len(k), i)[p]))))

#' True profile means as a class x indicator matrix, from a tidy truth table
true_means <- function(case) {
  truth <- case$truth
  if (identical(case$config, "lpa_covariates")) return(truth$means)
  if (identical(case$config, "ml_lpa")) truth <- truth$parameters %||% truth[[1L]]
  stopifnot(all(c("profile", "variable", "mean") %in% names(truth)))
  wide <- tapply(truth$mean, list(truth$profile, truth$variable), identity)
  wide[, order(colnames(wide)), drop = FALSE]
}

#' Membership labels as 1..K in the order of the truth
true_classes <- function(case) as.integer(factor(case$membership))

recover_one <- function(name) {
  case <- cases[[name]]
  fit <- suppressWarnings(montecarlo_fit(case))
  k <- case$spec$k
  truth_class <- true_classes(case)
  observed <- rowSums(!is.na(case$data[case$spec$vars])) > 0L
  estimated <- max.col(fit$subject_posteriors, ties.method = "first")
  # Label alignment from the classification itself, then reused for every
  # parameter of the case.
  candidates <- permutations(k)
  agreement <- vapply(candidates, \(p) mean(p[estimated[observed]] == truth_class[observed]),
                      numeric(1))
  order <- candidates[[which.max(agreement)]]
  accuracy <- max(agreement)
  mean_error <- if (isTRUE(case$spec$categorical)) NA_real_ else {
    truth <- true_means(case)
    estimate <- fit$means[order(order), , drop = FALSE]
    as.vector(estimate - truth)
  }
  probability_error <- if (!isTRUE(case$spec$categorical)) NA_real_ else {
    truth <- case$truth
    first <- if ("probability_first" %in% names(truth)) {
      tapply(truth$probability_first, list(truth$class, truth$item), identity)
    } else {
      ones <- truth[truth$category == 1, ]
      wide <- tapply(ones$probability, list(ones$latent_class, ones$indicator), identity)
      colnames(wide) <- sub("item_", "u", colnames(wide))
      wide
    }
    estimate <- vapply(colnames(first), function(item)
      fit$response_probabilities[[item]][order(order), 1L], numeric(k))
    as.vector(estimate - first)
  }
  # Two classes: the package reports the logits of class 1 against class 2;
  # the truth is of class 2 against class 1, so it is negated when the
  # package's class 1 is the true class 1.
  logit_error <- if (!identical(case$config, "lpa_covariates")) rep(NA_real_, 3L) else {
    coefficients <- get_results(fit, "coefficients")
    estimate <- coefficients$estimate[match(c("(Intercept)", "x1", "x2"), coefficients$term)]
    true_first <- order(order)[1L] == 1L
    (if (true_first) -estimate else estimate) - case$truth$logit
  }
  data.frame(case = name, config = case$config, accuracy = accuracy,
             probability_bias = if (all(is.na(probability_error))) NA_real_ else
               mean(probability_error),
             probability_rmse = if (all(is.na(probability_error))) NA_real_ else
               sqrt(mean(probability_error^2)),
             intercept_error = logit_error[1L], x1_error = logit_error[2L],
             x2_error = logit_error[3L],
             mean_bias = if (all(is.na(mean_error))) NA_real_ else mean(mean_error),
             mean_rmse = if (all(is.na(mean_error))) NA_real_ else sqrt(mean(mean_error^2)),
             max_abs_mean_error = if (all(is.na(mean_error))) NA_real_ else max(abs(mean_error)))
}

rows <- do.call(rbind, parallel::mclapply(names(cases), recover_one, mc.cores = cores,
                                          mc.preschedule = FALSE))
by_config <- aggregate(cbind(accuracy, mean_bias, mean_rmse, max_abs_mean_error,
                             probability_bias, probability_rmse) ~ config,
                       data = rows, FUN = mean, na.action = na.pass)
logits <- subset(rows, config == "lpa_covariates")
logit_summary <- data.frame(
  term = c("(Intercept)", "x1", "x2"), truth = c(-0.3, 0.9, -0.6),
  bias = c(mean(logits$intercept_error), mean(logits$x1_error), mean(logits$x2_error)),
  rmse = sqrt(c(mean(logits$intercept_error^2), mean(logits$x1_error^2),
                mean(logits$x2_error^2))))
utils::write.csv(rows, file.path("validation", "montecarlo-recovery.csv"), row.names = FALSE)

format_number <- function(v) ifelse(is.na(v), "--", formatC(v, format = "f", digits = 3))
lines <- c(
  "# Monte Carlo recovery against the generating parameters",
  "",
  "Generated by `validation/montecarlo-recovery.R` from the 120 cases of",
  "`equivalence/montecarlo/` (8 configurations x 15 replicates; seed 1, 20 starts).",
  "Classes are aligned to the true labels by the permutation that maximizes",
  "agreement of the modal classes, and that alignment is reused for the means.",
  "Accuracy is over rows with an observed indicator; mean errors are over every",
  "profile-by-indicator mean of a case, then averaged over its replicates.",
  "",
  "| Configuration | Classification accuracy | Mean bias | Mean RMSE | Largest absolute mean error | Item probability bias | Item probability RMSE |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %s | %s | %s | %s |", by_config$config,
          format_number(by_config$accuracy), format_number(by_config$mean_bias),
          format_number(by_config$mean_rmse), format_number(by_config$max_abs_mean_error),
          format_number(by_config$probability_bias), format_number(by_config$probability_rmse)),
  "",
  "Class-membership logits of `lpa_covariates` (class 2 against class 1), over its",
  "15 replicates:",
  "",
  "| Term | Truth | Bias | RMSE |",
  "|---|---|---|---|",
  sprintf("| %s | %s | %s | %s |", logit_summary$term, format_number(logit_summary$truth),
          format_number(logit_summary$bias), format_number(logit_summary$rmse)),
  "",
  "Per-case results: `validation/montecarlo-recovery.csv`.")
writeLines(lines, file.path("validation", "MONTECARLO_RECOVERY.md"))
print(by_config)
