# Interval coverage of three_step() with the step-one estimates held fixed
# (step_one = "fixed") and with their uncertainty added (step_one =
# "estimated", the first-order correction of Bakk, Oberski and Vermunt 2014),
# for a categorical and a numeric distal outcome, by BCH and the ML
# step-three model. Run from the package root:
#   Rscript validation/three-step-step-one-coverage.R        # 200 per cell
#   Rscript validation/three-step-step-one-coverage.R 5      # a quick pilot
#
# 16 cells: K = 2, 3 classes; C = 2, 3 categories; relative entropy about 0.6
# and 0.8 (profile separation 1.7 and 2.35 for K = 2, 1.5 and 2.05 for K = 3,
# from a pilot); n = 300 and 1000. Two Gaussian indicators with class means
# (k - 1) * separation on both, equal class sizes, and the outcome's
# distribution in class k drawn from fixed tables below; the numeric outcome
# is normal with mean k - 1 and standard deviation 1 in class k. The fitted
# classes are matched to the true ones by their mean on the first indicator.
# The design is that of three-step-categorical-coverage.R (the numeric draws
# make the data sets themselves differ). Declared before running: with
# step_one = "estimated", 95% coverage of every class parameter between 0.92
# and 0.98 in every cell (two Monte Carlo SEs of 0.015 at 200 data sets; with
# up to nine parameters per cell an occasional miss is chance), and at least
# as close to 0.95 as step_one = "fixed" where that misses. Data sets where
# BCH leaves a share outside [0, 1] are counted (the ML model then reports no
# standard errors), as are data sets whose ML step-three information is
# singular (no standard errors at all).
pkgload::load_all(".", quiet = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(arguments) > 0L) as.integer(arguments[[1L]]) else 200L

outcome_tables <- list(
  `2` = list(`2` = rbind(c(0.7, 0.3), c(0.3, 0.7)),
             `3` = rbind(c(0.5, 0.3, 0.2), c(0.2, 0.3, 0.5))),
  `3` = list(`2` = rbind(c(0.75, 0.25), c(0.5, 0.5), c(0.25, 0.75)),
             `3` = rbind(c(0.6, 0.3, 0.1), c(0.2, 0.6, 0.2), c(0.1, 0.3, 0.6))))
separation <- list(`2` = c(`0.6` = 1.7, `0.8` = 2.35), `3` = c(`0.6` = 1.5, `0.8` = 2.05))
cells <- expand.grid(classes = 2:3, categories = 2:3, entropy = c("0.6", "0.8"),
                     n = c(300L, 1000L), stringsAsFactors = FALSE)

replicate_cell <- function(replication, cell) {
  k <- cell$classes
  truth_table <- outcome_tables[[as.character(k)]][[as.character(cell$categories)]]
  delta <- separation[[as.character(k)]][[cell$entropy]]
  truth <- sample.int(k, cell$n, replace = TRUE)
  data <- data.frame(y1 = stats::rnorm(cell$n, (truth - 1) * delta),
                     y2 = stats::rnorm(cell$n, (truth - 1) * delta))
  labels <- paste0("c", seq_len(cell$categories))
  data$outcome <- factor(vapply(truth, function(class) {
    sample(labels, 1L, prob = truth_table[class, ])
  }, character(1)), levels = labels)
  quiet <- function(expr) {
    withCallingHandlers(expr, warning = function(condition) invokeRestart("muffleWarning"))
  }
  fit <- quiet(lpa(data, c("y1", "y2"), n_profiles = k, model = "EEI", n_starts = 3, seed = 1))
  matched <- rank(fit$means[, "y1"])
  entropy <- get_results(fit, "entropy")$relative_entropy[[1L]]
  data$z <- stats::rnorm(cell$n, truth - 1)
  settings <- expand.grid(method = c("bch", "ml"), step_one = c("fixed", "estimated"),
                          outcome = c("outcome", "z"), stringsAsFactors = FALSE)
  rows <- Map(function(method, step_one, outcome) {
    flagged <- FALSE
    # A step-three ML model whose information is singular has no standard
    # errors (latents_singular_information); the data set is counted, not
    # dropped silently.
    table <- tryCatch(withCallingHandlers(
      three_step(fit, data, outcome, method = method, vcov_type = "independent",
                 step_one = step_one),
      latents_bch_out_of_range = function(condition) {
        flagged <<- TRUE
        invokeRestart("muffleWarning")
      },
      warning = function(condition) invokeRestart("muffleWarning")),
      latents_singular_information = function(condition) NULL)
    if (is.null(table)) {
      return(data.frame(cell_classes = k, cell_categories = cell$categories,
                        cell_entropy = cell$entropy, cell_n = cell$n,
                        replication = replication, method = method, step_one = step_one,
                        outcome = if (identical(outcome, "z")) "numeric" else "categorical",
                        relative_entropy = entropy, out_of_range = flagged, singular = TRUE,
                        true_class = NA_integer_, category = NA_character_,
                        truth = NA_real_, estimate = NA_real_, covered = NA))
    }
    table$true_class <- matched[table$class]
    numeric_outcome <- identical(outcome, "z")
    category <- if (numeric_outcome) rep("mean", nrow(table)) else table$category
    truth_value <- if (numeric_outcome) table$true_class - 1 else
      truth_table[cbind(table$true_class, match(table$category, labels))]
    data.frame(cell_classes = k, cell_categories = cell$categories, cell_entropy = cell$entropy,
               cell_n = cell$n, replication = replication, method = method,
               step_one = step_one, outcome = if (numeric_outcome) "numeric" else "categorical",
               relative_entropy = entropy, out_of_range = flagged, singular = FALSE,
               true_class = table$true_class, category = category,
               truth = truth_value, estimate = table$estimate,
               covered = table$conf_low <= truth_value & truth_value <= table$conf_high)
  }, settings$method, settings$step_one, settings$outcome)
  do.call(rbind, rows)
}

started <- proc.time()[["elapsed"]]
cluster <- parallel::makePSOCKcluster(2L)
parallel::clusterSetRNGStream(cluster, 20261102)
invisible(parallel::clusterEvalQ(cluster, {
  Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  pkgload::load_all(".", quiet = TRUE)
}))
parallel::clusterExport(cluster, c("outcome_tables", "separation", "replicate_cell"))
results <- do.call(rbind, lapply(seq_len(nrow(cells)), function(index) {
  do.call(rbind, parallel::parLapply(cluster, seq_len(datasets), replicate_cell,
                                     cell = cells[index, ]))
}))
parallel::stopCluster(cluster)

# One row per cell, outcome, method and step_one: the largest absolute bias
# over the class parameters, the range of their coverage, how often BCH went
# out of range and how often the ML information was singular.
summary_table <- do.call(rbind, lapply(
  split(results, with(results, list(cell_classes, cell_categories, cell_entropy, cell_n,
                                     outcome, method, step_one)), drop = TRUE),
  function(part) {
    # Coverage over the data sets that have an interval: an ML maximum on
    # the bound, or a BCH share outside [0, 1], has none.
    by_cell <- aggregate(cbind(bias = estimate - truth, coverage = covered) ~ true_class + category,
                         part, mean, na.rm = TRUE, na.action = stats::na.pass)
    with(part, data.frame(
      classes = cell_classes[1L], categories = cell_categories[1L], entropy = cell_entropy[1L],
      n = cell_n[1L], outcome = outcome[1L], method = method[1L], step_one = step_one[1L],
      mean_entropy = mean(relative_entropy),
      max_abs_bias = max(abs(by_cell$bias)), min_coverage = min(by_cell$coverage, na.rm = TRUE),
      max_coverage = max(by_cell$coverage, na.rm = TRUE),
      out_of_range = mean(tapply(out_of_range, replication, any)),
      singular = mean(tapply(singular, replication, any))))
  }))
summary_table <- summary_table[with(summary_table, order(
  outcome, classes, categories, entropy, n, method, step_one)), ]
summary_table$within_92_98 <- with(summary_table, min_coverage >= 0.92 & max_coverage <= 0.98)
print(summary_table, digits = 3, row.names = FALSE)
cat(sprintf("%.1f min\n", (proc.time()[["elapsed"]] - started) / 60))
if (datasets >= 100L) {
  utils::write.csv(summary_table, "validation/three-step-step-one-coverage.csv",
                   row.names = FALSE)
  utils::write.csv(results, "validation/three-step-step-one-coverage-replicates.csv",
                   row.names = FALSE)
}
