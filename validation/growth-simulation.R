# Predeclared simulation of the growth mixture model's estimates, standard
# errors and person-level R3STEP (phases G2 and G4 of todo.MD). The registry
# (growth-simulation-registry.csv) is written BEFORE any replicate runs and the
# script refuses to run if it changes. Gates as in ADDITIVE_SIMULATION.md:
# availability >= 98%, |bias| <= 0.1 empirical SD, observed SE ratio 0.9-1.1,
# observed 95% coverage 0.925-0.975.
#
# Class membership depends on a person covariate x, logit P(class 2) =
# 0.2 + 0.8 x. The growth mixture is fitted without it (its marginal class
# shares are then the truth, E[P(class | x)]), and r3step(fit, data, "x")
# estimates the slope, whose truth is -0.8 or 0.8 for `class_1` against the
# reference class depending on which true class the fit calls class 1.
#
# Run from the project root: Rscript validation/growth-simulation.R
# Outputs: growth-simulation-{replicates,summary}.csv, growth-simulation.txt

devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
out_dir <- Sys.getenv("GROWTH_SIM_OUT", "validation")
n_datasets <- as.integer(Sys.getenv("GROWTH_SIM_N", "1000"))
n_cores <- as.integer(Sys.getenv("GROWTH_SIM_CORES", "8"))

base_class <- list(intercept = c(2, 5), slope = c(0.6, -0.4),
                   variance = list(c(1, 0.09, 0.3), c(0.6, 0.05, -0.2)),
                   sigma = c(0.7, 0.8))
cells <- list(
  balanced = list(J = 300L, times = 0:5, jitter = 0, missing = 0,
                  classes = base_class),
  unbalanced = list(J = 300L, times = 0:5, jitter = 0.3, missing = 0.25,
                    classes = base_class),
  weak = list(J = 500L, times = 0:5, jitter = 0, missing = 0,
              classes = utils::modifyList(base_class,
                                          list(intercept = c(2, 3.2), slope = c(0.6, 0.1)))),
  small_variance = list(J = 300L, times = 0:5, jitter = 0, missing = 0,
                        classes = utils::modifyList(base_class, list(
                          variance = list(c(1, 0.01, 0.3), c(0.6, 0.01, -0.2))))))
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  data.frame(cell = names(cells)[k], persons = cell$J,
             times = paste(cell$times, collapse = ","), jitter = cell$jitter,
             missing = cell$missing,
             intercept = paste(cell$classes$intercept, collapse = ","),
             slope = paste(cell$classes$slope, collapse = ","),
             variance = paste(vapply(cell$classes$variance, paste, character(1),
                                     collapse = ","), collapse = " | "),
             sigma = paste(cell$classes$sigma, collapse = ","),
             membership = "0.2,0.8", n_starts = 5L, base_seed = 7e6L + 1e6L * k,
             stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "growth-simulation-registry.csv")
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- utils::read.csv(text = paste(utils::capture.output(
    utils::write.csv(registry, row.names = FALSE)), collapse = "\n"),
    stringsAsFactors = FALSE)
  stopifnot("registry changed after it was declared" = identical(stored, current))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}

covariance_of <- function(v) {
  # variance of intercept, variance of slope, their correlation
  off <- v[3L] * sqrt(v[1L] * v[2L])
  matrix(c(v[1L], off, off, v[2L]), 2L)
}

draw <- function(cell, seed) {
  set.seed(seed)
  x <- stats::rnorm(cell$J)
  class <- 1L + (stats::runif(cell$J) < stats::plogis(0.2 + 0.8 * x))
  effects <- t(vapply(seq_len(cell$J), function(i) {
    as.vector(t(chol(covariance_of(cell$classes$variance[[class[i]]]))) %*%
                stats::rnorm(2L))
  }, numeric(2)))
  n_times <- length(cell$times)
  frame <- data.frame(id = rep(seq_len(cell$J), each = n_times),
                      time = rep(cell$times, cell$J) +
                        stats::runif(cell$J * n_times, -cell$jitter, cell$jitter),
                      x = rep(x, each = n_times))
  k <- class[frame$id]
  frame$y <- cell$classes$intercept[k] + cell$classes$slope[k] * frame$time +
    effects[frame$id, 1L] + effects[frame$id, 2L] * frame$time +
    stats::rnorm(nrow(frame), 0, cell$classes$sigma[k])
  # Missing occasions at random; every person keeps their first two.
  keep <- stats::runif(nrow(frame)) >= cell$missing |
    stats::ave(frame$time, frame$id, FUN = seq_along) <= 2
  frame[keep, ]
}

truth_rows <- function(cell) {
  share_2 <- stats::integrate(function(x) stats::plogis(0.2 + 0.8 * x) * stats::dnorm(x),
                              -Inf, Inf)$value
  do.call(rbind, lapply(1:2, function(k) {
    v <- cell$classes$variance[[k]]
    label <- paste0("class_", k)
    data.frame(
      block = c("coefficient", "coefficient", "residual_sd", "random_variance",
                "random_variance", "random_sd", "random_sd", "random_correlation",
                "share"),
      class = label,
      term = c("(Intercept)", "time", "residual", "(Intercept)", "time",
               "(Intercept)", "time", "(Intercept):time", "share"),
      truth = c(cell$classes$intercept[k], cell$classes$slope[k], cell$classes$sigma[k],
                v[1L], v[2L], sqrt(v[1L]), sqrt(v[2L]), v[3L],
                if (k == 1L) 1 - share_2 else share_2),
      stringsAsFactors = FALSE)
  }))
}

one_replicate <- function(cell_name, i) {
  cell <- cells[[cell_name]]
  seed <- registry$base_seed[registry$cell == cell_name] + i
  data <- draw(cell, seed)
  truth <- truth_rows(cell)
  failure <- function(reason) {
    data.frame(cell = cell_name, replicate = i, truth[c("block", "class", "term")],
               truth = truth$truth, estimate = NA_real_, se_observed = NA_real_,
               se_robust = NA_real_, covered_observed = NA, covered_robust = NA,
               available = FALSE, reason = reason, stringsAsFactors = FALSE)
  }
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                       random = ~ 1 + time, n_starts = 5, seed = i, vcov_type = "none"),
    warning = function(w) {
      warnings <<- c(warnings, class(w)[1L])
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) return(list(rows = failure(conditionMessage(fit)), r3 = NULL))
  if (!isTRUE(fit$converged)) return(list(rows = failure("unconverged"), r3 = NULL))
  if (length(fit$degenerate_classes %||% character()) > 0L) {
    return(list(rows = failure("degenerate covariance"), r3 = NULL))
  }
  # Name the fitted classes by the true class they estimate.
  signature <- function(intercepts, slopes) cbind(intercepts, slopes)
  order <- .multilpa_match_order(
    signature(cell$classes$intercept, cell$classes$slope),
    signature(fit$params$beta[1L, ], fit$params$beta[2L, ]))
  true_label <- paste0("class_", match(seq_len(2L), order))
  tables <- tryCatch(list(
    observed = parameter_inference(fit, vcov_type = "observed"),
    robust = parameter_inference(fit, vcov_type = "robust")), error = function(e) e)
  if (inherits(tables, "error")) return(list(rows = failure(conditionMessage(tables)), r3 = NULL))
  relabel <- function(table) {
    table <- subset(table, block != "membership")
    table$class <- ifelse(table$class == "all", "all",
                          true_label[match(table$class, c("class_1", "class_2"))])
    merge(truth, table, by = c("block", "class", "term"), sort = FALSE)
  }
  observed <- relabel(tables$observed)
  robust <- relabel(tables$robust)
  robust <- robust[match(paste(observed$block, observed$class, observed$term),
                         paste(robust$block, robust$class, robust$term)), ]
  covered <- function(table) table$conf_low <= table$truth & table$truth <= table$conf_high
  rows <- data.frame(cell = cell_name, replicate = i,
                     observed[c("block", "class", "term", "truth", "estimate")],
                     se_observed = observed$standard_error,
                     se_robust = robust$standard_error,
                     covered_observed = covered(observed), covered_robust = covered(robust),
                     available = is.finite(observed$standard_error),
                     reason = NA_character_, stringsAsFactors = FALSE)
  # R3STEP: log odds of fitted class 1 against class 2 per unit of x.
  r3 <- tryCatch(r3step(fit, data, "x"), error = function(e) e)
  r3_row <- if (inherits(r3, "error")) {
    data.frame(cell = cell_name, replicate = i, truth = NA_real_, estimate = NA_real_,
               se = NA_real_, covered = NA, reason = conditionMessage(r3))
  } else {
    slope <- subset(r3, term == "x")
    truth_slope <- if (identical(true_label[1L], "class_1")) -0.8 else 0.8
    data.frame(cell = cell_name, replicate = i, truth = truth_slope,
               estimate = slope$estimate, se = slope$standard_error,
               covered = slope$conf_low <= truth_slope & truth_slope <= slope$conf_high,
               reason = NA_character_)
  }
  list(rows = rows, r3 = r3_row)
}

started <- Sys.time()
results <- unlist(lapply(names(cells), function(cell_name) {
  parallel::mclapply(seq_len(n_datasets), function(i) one_replicate(cell_name, i),
                     mc.cores = n_cores, mc.set.seed = TRUE)
}), recursive = FALSE)
errors <- vapply(results, inherits, logical(1), what = "try-error")
stopifnot("a worker failed" = !any(errors))
replicates <- do.call(rbind, lapply(results, `[[`, "rows"))
r3 <- do.call(rbind, lapply(results, `[[`, "r3"))
utils::write.csv(replicates, file.path(out_dir, "growth-simulation-replicates.csv"),
                 row.names = FALSE)
utils::write.csv(r3, file.path(out_dir, "growth-simulation-r3step.csv"), row.names = FALSE)

summarise <- function(rows) {
  ok <- rows[rows$available, ]
  sd_empirical <- stats::sd(ok$estimate)
  bias <- mean(ok$estimate) - rows$truth[1L]
  data.frame(availability = mean(rows$available),
             bias = bias, bias_sd = bias / sd_empirical, empirical_sd = sd_empirical,
             se_ratio_observed = mean(ok$se_observed) / sd_empirical,
             se_ratio_robust = mean(ok$se_robust) / sd_empirical,
             coverage_observed = mean(ok$covered_observed),
             coverage_robust = mean(ok$covered_robust))
}
keys <- unique(replicates[c("cell", "block", "class", "term")])
summary <- do.call(rbind, lapply(seq_len(nrow(keys)), function(r) {
  rows <- merge(keys[r, ], replicates, by = c("cell", "block", "class", "term"))
  cbind(keys[r, ], truth = rows$truth[1L], summarise(rows))
}))
summary$pass <- summary$availability >= 0.98 & abs(summary$bias_sd) <= 0.1 &
  summary$se_ratio_observed >= 0.9 & summary$se_ratio_observed <= 1.1 &
  summary$coverage_observed >= 0.925 & summary$coverage_observed <= 0.975
rownames(summary) <- NULL
utils::write.csv(summary, file.path(out_dir, "growth-simulation-summary.csv"),
                 row.names = FALSE)

r3_summary <- do.call(rbind, lapply(names(cells), function(cell_name) {
  rows <- r3[r3$cell == cell_name & is.finite(r3$estimate), ]
  aligned <- rows$estimate * sign(rows$truth)
  data.frame(cell = cell_name, available = nrow(rows) / n_datasets,
             bias_sd = (mean(aligned) - 0.8) / stats::sd(aligned),
             se_ratio = mean(rows$se) / stats::sd(aligned),
             coverage = mean(rows$covered))
}))
reasons <- table(replicates$reason[!duplicated(replicates[c("cell", "replicate")])],
                 useNA = "no")
sink(file.path(out_dir, "growth-simulation.txt"))
cat(sprintf("Growth mixture simulation: %d datasets per cell, %s elapsed\n\n",
            n_datasets, format(round(Sys.time() - started, 1))))
print(summary, digits = 3, row.names = FALSE)
cat("\nR3STEP slope on x (true |slope| 0.8, robust intervals)\n")
print(r3_summary, digits = 3, row.names = FALSE)
cat("\nUnavailable replicates by reason\n")
print(reasons)
cat("\n")
print(utils::sessionInfo())
sink()
cat(readLines(file.path(out_dir, "growth-simulation.txt")), sep = "\n")
