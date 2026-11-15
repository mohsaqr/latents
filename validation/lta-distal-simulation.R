# Predeclared recovery study of three_step(occasion =): a person-level
# distal outcome related to the profile an lta() chain occupies at one
# occasion, or to the pair of profiles at two occasions, by BCH and ML.
# The registry (lta-distal-simulation-registry.csv) is written BEFORE any
# replicate runs and the script refuses to run if it changes.
#
# Gates, fixed before the first replicate (those of the simulab recovery kit,
# equivalence/simulab-recovery/README.md):
# - bias: |mean(estimate - truth)| <= 3.5 Monte Carlo standard errors;
# - coverage: the share of 95% intervals containing the truth lies in the
#   central 99.9% binomial range for R replicates at 0.95;
# - every replicate fits (a failure is counted, not dropped).
# The uncorrected methods (modal, proportional) are recorded for contrast and
# are not gated: they are biased towards each other by design.
#
# Data are drawn in base R (simulab was not available on the machine that
# produced the retained results): a two-profile first-order chain, three
# normal indicators per occasion, a distal outcome normal given the profile
# at occasion 3 and another given the pair at occasions 2 and 3.
#
# Run from the project root:
#   Rscript validation/lta-distal-simulation.R            # 500 replicates per cell
#   LTA_DISTAL_SIM_N=100 Rscript validation/lta-distal-simulation.R
# Outputs: lta-distal-simulation-{replicates,summary}.csv, lta-distal-simulation.txt

devtools::load_all(".", quiet = TRUE)
out_dir <- "validation"
n_datasets <- as.integer(Sys.getenv("LTA_DISTAL_SIM_N", "500"))

cells <- list(
  separated = list(J = 500L, T = 4L, item_gap = 2.0, initial = c(0.55, 0.45),
                   stay = c(0.8, 0.85), distal = c(0, 1), distal_sd = c(1, 1.2),
                   pair = c(0, 0.5, 1.5, 2), pair_sd = 1),
  overlapping = list(J = 600L, T = 4L, item_gap = 1.4, initial = c(0.55, 0.45),
                     stay = c(0.8, 0.85), distal = c(0, 1), distal_sd = c(1, 1.2),
                     pair = c(0, 0.5, 1.5, 2), pair_sd = 1))
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  data.frame(cell = names(cells)[k], persons = cell$J, occasions = cell$T,
             item_gap = cell$item_gap, initial = paste(cell$initial, collapse = ","),
             stay = paste(cell$stay, collapse = ","),
             distal = paste(cell$distal, collapse = ","),
             distal_sd = paste(cell$distal_sd, collapse = ","),
             pair = paste(cell$pair, collapse = ","), pair_sd = cell$pair_sd,
             base_seed = 7e7L + 1e6L * k, stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "lta-distal-simulation-registry.csv")
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- utils::read.csv(text = paste(utils::capture.output(
    utils::write.csv(registry, row.names = FALSE)), collapse = "\n"),
    stringsAsFactors = FALSE)
  stopifnot("registry changed after it was declared" = identical(stored, current))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}

draw <- function(cell, seed) {
  set.seed(seed)
  state <- matrix(NA_integer_, cell$J, cell$T)
  state[, 1L] <- ifelse(stats::runif(cell$J) < cell$initial[1L], 1L, 2L)
  for (t in seq_len(cell$T)[-1L]) {  # simulation only
    stays <- stats::runif(cell$J) < cell$stay[state[, t - 1L]]
    state[, t] <- ifelse(stays, state[, t - 1L], 3L - state[, t - 1L])
  }
  frame <- data.frame(id = rep(seq_len(cell$J), each = cell$T),
                      time = rep(seq_len(cell$T), cell$J))
  s <- as.vector(t(state))
  half <- cell$item_gap / 2
  for (i in 1:3) frame[[paste0("y", i)]] <- stats::rnorm(length(s), c(-half, half)[s], 1)
  z <- stats::rnorm(cell$J, cell$distal[state[, 3L]], cell$distal_sd[state[, 3L]])
  pair <- (state[, 2L] - 1L) * 2L + state[, 3L]
  z23 <- stats::rnorm(cell$J, cell$pair[pair], cell$pair_sd)
  frame$z <- rep(z, each = cell$T)
  frame$z23 <- rep(z23, each = cell$T)
  frame
}

run_one <- function(cell_name, replicate) {
  cell <- cells[[cell_name]]
  seed <- registry$base_seed[registry$cell == cell_name] + replicate
  data <- draw(cell, seed)
  conditions <- character()
  record <- function(condition) {
    conditions <<- c(conditions, class(condition)[1L])
    invokeRestart(if (inherits(condition, "warning")) "muffleWarning" else "muffleMessage")
  }
  quiet <- function(expr) withCallingHandlers(expr, warning = record, message = record)
  fit <- quiet(lta(data, c("y1", "y2", "y3"), "id", n_profiles = 2, time = "time",
                   n_starts = 3, seed = seed))
  # Label switching: profile 1 has the lower first indicator.
  swap <- fit$means[1L, "y1"] > fit$means[2L, "y1"]
  one_label <- function(k) if (swap) 3L - k else k
  pair_label <- function(cell) {
    from <- one_label((cell - 1L) %/% 2L + 1L)
    to <- one_label((cell - 1L) %% 2L + 1L)
    (from - 1L) * 2L + to
  }
  rows <- list()
  for (method in c("bch", "ml", "modal", "proportional")) {
    single <- quiet(three_step(fit, data, "z", occasion = 3, method = method))
    paired <- quiet(three_step(fit, data, "z23", occasion = c(2, 3), method = method))
    rows[[length(rows) + 1L]] <- data.frame(
      target = "occasion 3", method = method,
      parameter = paste0("profile ", one_label(single$class)),
      truth = cell$distal[one_label(single$class)], estimate = single$estimate,
      se = single$standard_error, low = single$conf_low, high = single$conf_high)
    rows[[length(rows) + 1L]] <- data.frame(
      target = "occasions 2-3", method = method,
      parameter = paste0("pair ", c("1-1", "1-2", "2-1", "2-2")[pair_label(paired$class)]),
      truth = cell$pair[pair_label(paired$class)], estimate = paired$estimate,
      se = paired$standard_error, low = paired$conf_low, high = paired$conf_high)
  }
  data.frame(cell = cell_name, replicate = replicate, do.call(rbind, rows),
             converged = fit$converged, conditions = paste(unique(conditions), collapse = ";"),
             stringsAsFactors = FALSE)
}

started <- Sys.time()
jobs <- expand.grid(replicate = seq_len(n_datasets), cell = names(cells),
                    stringsAsFactors = FALSE)
replicates <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(k) {
  tryCatch(run_one(jobs$cell[k], jobs$replicate[k]), error = function(error) {
    data.frame(cell = jobs$cell[k], replicate = jobs$replicate[k], target = NA_character_,
               method = NA_character_, parameter = NA_character_, truth = NA_real_,
               estimate = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_,
               converged = NA, conditions = paste("error:", conditionMessage(error)),
               stringsAsFactors = FALSE)
  })
}))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
utils::write.csv(replicates, file.path(out_dir, "lta-distal-simulation-replicates.csv"),
                 row.names = FALSE)

usable <- subset(replicates, !is.na(parameter))
keys <- unique(usable[, c("cell", "target", "method", "parameter")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  rows <- merge(keys[k, ], usable)
  attempted <- sum(replicates$cell == keys$cell[k] &
                     !duplicated(replicates[c("cell", "replicate")]))
  r <- nrow(rows)
  errors <- rows$estimate - rows$truth
  mcse <- stats::sd(errors) / sqrt(r)
  covered <- mean(rows$low <= rows$truth & rows$truth <= rows$high, na.rm = TRUE)
  bounds <- stats::qbinom(c(0.0005, 0.9995), r, 0.95) / r
  gated <- keys$method[k] %in% c("bch", "ml")
  data.frame(keys[k, ], truth = rows$truth[1L], replicates = r, attempted = attempted,
             mean_estimate = mean(rows$estimate), bias = mean(errors), mcse = mcse,
             bias_ok = abs(mean(errors)) <= 3.5 * mcse,
             sd_estimate = stats::sd(rows$estimate), mean_se = mean(rows$se, na.rm = TRUE),
             coverage = covered, coverage_low = bounds[1L], coverage_high = bounds[2L],
             coverage_ok = covered >= bounds[1L] && covered <= bounds[2L],
             gated = gated,
             meets = !gated || (r == attempted && abs(mean(errors)) <= 3.5 * mcse &&
                                  covered >= bounds[1L] && covered <= bounds[2L]),
             row.names = NULL, stringsAsFactors = FALSE)
}))
utils::write.csv(summary_table, file.path(out_dir, "lta-distal-simulation-summary.csv"),
                 row.names = FALSE)
failed <- sum(is.na(replicates$parameter))
writeLines(c(
  sprintf("three_step(occasion =) recovery: %d replicates per cell, %d cells, %.1f minutes",
          n_datasets, length(cells), elapsed),
  sprintf("replicates that failed: %d", failed),
  sprintf("gated rows: %d, meeting every gate: %d", sum(summary_table$gated),
          sum(summary_table$gated & summary_table$meets)),
  "", utils::capture.output(print(summary_table[, c(
    "cell", "target", "method", "parameter", "truth", "mean_estimate", "bias", "mcse",
    "sd_estimate", "mean_se", "coverage", "meets")], row.names = FALSE, digits = 3))),
  file.path(out_dir, "lta-distal-simulation.txt"))
cat(readLines(file.path(out_dir, "lta-distal-simulation.txt")), sep = "\n")
