# Predeclared simulation of the general transition model's standard errors.
# The registry (lta-simulation-registry.csv) is written BEFORE any replicate
# runs and the script refuses to run if it changes. Gates as in
# ADDITIVE_SIMULATION.md: availability >= 98%, |bias| <= 0.1 SD, observed SE
# ratio 0.9-1.1, observed 95% coverage 0.925-0.975.
# Run from the project root: Rscript validation/lta-simulation.R
# Outputs: lta-simulation-{replicates,summary}.csv, lta-simulation.txt

devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
out_dir <- "validation"
n_datasets <- as.integer(Sys.getenv("LTA_SIM_N", "1000"))

cells <- list(
  covariate = list(J = 300L, T = 5L, items = "continuous", order = 1L,
                   transitions = "homogeneous", initial = c(0.3, 0.8),
                   # per origin: intercept, slope on time-varying z
                   move = rbind(c(-1.2, 0.8), c(-1.0, -0.6)), second = NULL),
  occasion = list(J = 300L, T = 4L, items = "continuous", order = 1L,
                  transitions = "occasion", initial = 0.2,
                  move = rbind(c(-1.5, -0.8, -0.3), c(-1.2, -1.2, -0.5)), second = NULL),
  second_order = list(J = 300L, T = 5L, items = "continuous", order = 2L,
                      transitions = "homogeneous", initial = 0,
                      move = rbind(-1.0, -1.0),
                      # origin pairs (1,1), (1,2), (2,1), (2,2): log odds of moving
                      second = c(-2.0, 0.3, 0.3, -2.0)),
  categorical = list(J = 400L, T = 3L, items = "binary", order = 1L,
                     transitions = "homogeneous", initial = c(0, 0.6),
                     move = rbind(c(-1.3, 0.7), c(-1.1, -0.5)), second = NULL))
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  data.frame(cell = names(cells)[k], groups = cell$J, occasions = cell$T,
             items = cell$items, order = cell$order, transitions = cell$transitions,
             initial = paste(cell$initial, collapse = ","),
             move = paste(apply(cell$move, 1L, paste, collapse = ","), collapse = " | "),
             second = paste(cell$second, collapse = ","), base_seed = 9e6L + 1e6L * k,
             stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "lta-simulation-registry.csv")
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
  w <- stats::rnorm(cell$J)
  z <- matrix(stats::rnorm(cell$J * cell$T), cell$J)
  initial_logit <- if (length(cell$initial) == 2L) cell$initial[1] + cell$initial[2] * w else
    rep(cell$initial[1], cell$J)
  state <- matrix(NA_integer_, cell$J, cell$T)
  state[, 1L] <- ifelse(stats::runif(cell$J) < stats::plogis(initial_logit), 1L, 2L)
  move_logit <- function(t, origin) {
    coef <- cell$move[origin, ]
    if (identical(cell$transitions, "occasion")) coef[t - 1L] else
      if (length(coef) == 2L) coef[1] + coef[2] * z[, t] else rep(coef[1], cell$J)
  }
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(cell$T)[-1L]) {  # simulation only
    logit <- if (cell$order == 2L && t >= 3L) {
      cell$second[(state[, t - 2L] - 1L) * 2L + state[, t - 1L]]
    } else ifelse(state[, t - 1L] == 1L, move_logit(t, 1L), move_logit(t, 2L))
    moves <- stats::runif(cell$J) < stats::plogis(logit)
    state[, t] <- ifelse(moves, 3L - state[, t - 1L], state[, t - 1L])
  }
  frame <- data.frame(id = rep(seq_len(cell$J), each = cell$T),
                      time = rep(seq_len(cell$T), cell$J),
                      w = rep(w, each = cell$T), z = as.vector(t(z)))
  s <- as.vector(t(state))
  if (identical(cell$items, "continuous")) {
    frame$y1 <- stats::rnorm(length(s), c(-1, 1)[s], 1)
    frame$y2 <- stats::rnorm(length(s), c(-0.8, 0.8)[s], 1)
  } else {
    probability <- rbind(c(0.15, 0.2, 0.25, 0.2, 0.15), c(0.8, 0.75, 0.85, 0.7, 0.8))
    invisible(lapply(1:5, function(i) {
      frame[[paste0("u", i)]] <<- factor(as.integer(stats::runif(length(s)) <
                                                      probability[s, i]))
    }))
  }
  frame
}

# True coefficients in the engine's parameterization (initial: log odds of
# profile 1 against 2; transitions: log odds of moving rather than staying).
truth_table <- function(cell) {
  initial_terms <- if (length(cell$initial) == 2L) c("(Intercept)", "w") else "(Intercept)"
  transition_terms <- if (identical(cell$transitions, "occasion"))
    paste0("occasion_", seq_len(ncol(cell$move)) + 1L) else
      if (ncol(cell$move) == 2L) c("(Intercept)", "z") else "(Intercept)"
  rbind(
    data.frame(block = "initial", key = paste("initial", initial_terms),
               truth = cell$initial),
    data.frame(block = "transition",
               key = paste("transition", rep(c("1->2", "2->1"), each = length(transition_terms)),
                           transition_terms),
               truth = as.vector(t(cell$move))),
    if (!is.null(cell$second)) data.frame(
      block = "second", key = paste("second", c("1-1->2", "1-2->1", "2-1->2", "2-2->1")),
      truth = cell$second))
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
  vars <- if (identical(cell$items, "continuous")) c("y1", "y2") else paste0("u", 1:5)
  fit <- withCallingHandlers(lta(
    data, vars, "id", n_profiles = 2, time = "time",
    categorical = if (identical(cell$items, "binary")) vars else character(),
    transitions = cell$transitions,
    transition_covariates = if (ncol(cell$move) == 2L && cell$transitions == "homogeneous") "z" else character(),
    initial_covariates = if (length(cell$initial) == 2L) "w" else character(),
    order = cell$order, n_starts = 5, seed = seed), warning = record, message = record)
  tables <- lapply(c("observed", "robust"), function(type) {
    withCallingHandlers(list(
      transition = get_results(fit, "transition_coefficients", vcov_type = type),
      initial = get_results(fit, "initial_coefficients", vcov_type = type),
      second = if (cell$order == 2L) get_results(fit, "second_order_coefficients",
                                                 vcov_type = type)),
      warning = record, message = record)
  })
  # Label switching: profile 1 has the lower first indicator.
  swap <- if (identical(cell$items, "continuous")) {
    fit$measurement[[1L]]$means[1L, "y1"] > fit$measurement[[1L]]$means[2L, "y1"]
  } else {
    fit$measurement[[1L]]$response_probabilities[[1L]][1L, 2L] >
      fit$measurement[[1L]]$response_probabilities[[1L]][2L, 2L]
  }
  relabel <- function(profile) {
    k <- as.integer(sub("profile_", "", profile))
    if (swap) 3L - k else k
  }
  keyed <- function(tables_one) {
    t <- tables_one$transition
    transition <- data.frame(key = paste("transition", paste0(relabel(t$from), "->",
                                                              relabel(t$to)), t$term),
                             estimate = t$estimate, se = t$standard_error,
                             low = t$conf_low, high = t$conf_high)
    i <- tables_one$initial
    # Initial logits are against the last profile: swapping flips their sign.
    initial <- data.frame(key = paste("initial", i$term),
                          estimate = if (swap) -i$estimate else i$estimate,
                          se = i$standard_error,
                          low = if (swap) -i$conf_high else i$conf_low,
                          high = if (swap) -i$conf_low else i$conf_high)
    second <- if (is.null(tables_one$second)) NULL else {
      s2 <- tables_one$second
      data.frame(key = paste("second", paste0(relabel(s2$previous), "-", relabel(s2$from),
                                              "->", relabel(s2$to))),
                 estimate = s2$estimate, se = s2$standard_error, low = s2$conf_low,
                 high = s2$conf_high)
    }
    rbind(transition, initial, second)
  }
  observed <- keyed(tables[[1L]])
  robust <- keyed(tables[[2L]])
  truth <- truth_table(cell)
  rows <- match(truth$key, observed$key)
  stopifnot("every true coefficient must be matched" = !anyNA(rows))
  data.frame(cell = cell_name, replicate = replicate, truth,
             estimate = observed$estimate[rows], se_observed = observed$se[rows],
             low_observed = observed$low[rows], high_observed = observed$high[rows],
             se_robust = robust$se[rows], low_robust = robust$low[rows],
             high_robust = robust$high[rows], converged = fit$converged,
             boundary = fit$boundary, conditions = paste(unique(conditions), collapse = ";"),
             stringsAsFactors = FALSE)
}

started <- Sys.time()
jobs <- expand.grid(cell = names(cells), replicate = seq_len(n_datasets),
                    stringsAsFactors = FALSE)
replicates <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  tryCatch(run_one(jobs$cell[k], jobs$replicate[k]), error = function(error) {
    data.frame(cell = jobs$cell[k], replicate = jobs$replicate[k], block = NA_character_,
               key = NA_character_, truth = NA_real_, estimate = NA_real_,
               se_observed = NA_real_, low_observed = NA_real_, high_observed = NA_real_,
               se_robust = NA_real_, low_robust = NA_real_, high_robust = NA_real_,
               converged = NA, boundary = NA,
               conditions = paste("error:", conditionMessage(error)),
               stringsAsFactors = FALSE)
  })
}, mc.cores = max(1L, parallel::detectCores() - 1L), mc.set.seed = TRUE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
utils::write.csv(replicates, file.path(out_dir, "lta-simulation-replicates.csv"),
                 row.names = FALSE)
usable <- subset(replicates, !is.na(key))
keys <- unique(usable[, c("cell", "key")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  rows <- subset(usable, cell == keys$cell[k] & key == keys$key[k])
  ok <- is.finite(rows$se_observed)
  okr <- is.finite(rows$se_robust)
  sd_estimate <- stats::sd(rows$estimate[ok])
  data.frame(cell = keys$cell[k], key = keys$key[k], truth = rows$truth[1L],
             availability = sum(ok) / n_datasets,
             bias_over_sd = (mean(rows$estimate[ok]) - rows$truth[1L]) / sd_estimate,
             se_ratio_observed = mean(rows$se_observed[ok]) / sd_estimate,
             coverage_observed = mean(rows$truth[ok] >= rows$low_observed[ok] &
                                        rows$truth[ok] <= rows$high_observed[ok]),
             se_ratio_robust = mean(rows$se_robust[okr]) / stats::sd(rows$estimate[okr]),
             coverage_robust = mean(rows$truth[okr] >= rows$low_robust[okr] &
                                      rows$truth[okr] <= rows$high_robust[okr]),
             stringsAsFactors = FALSE)
}))
summary_table$meets_release <- with(summary_table,
  availability >= 0.98 & abs(bias_over_sd) <= 0.1 & se_ratio_observed >= 0.9 &
    se_ratio_observed <= 1.1 & coverage_observed >= 0.925 & coverage_observed <= 0.975)
utils::write.csv(summary_table, file.path(out_dir, "lta-simulation-summary.csv"),
                 row.names = FALSE)
status <- do.call(rbind, lapply(names(cells), function(name) {
  per <- unique(subset(replicates, cell == name)[, c("replicate", "converged", "boundary",
                                                     "conditions")])
  data.frame(cell = name, attempted = n_datasets,
             errors = sum(startsWith(per$conditions, "error:")),
             converged = sum(per$converged %in% TRUE), boundary = sum(per$boundary %in% TRUE))
}))
sink(file.path(out_dir, "lta-simulation.txt"))
cat(sprintf("LTA extension simulation: %d datasets per cell, %.1f minutes; %s\n\n",
            n_datasets, elapsed, R.version.string))
print(status, row.names = FALSE)
cat("\n")
print(summary_table, digits = 3, row.names = FALSE)
sink()
cat(sprintf("done in %.1f minutes\n", elapsed))
