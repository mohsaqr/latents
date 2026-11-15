# Predeclared recovery study of the three-step latent Markov model with group
# classes: r3step(initial_covariates =, transition_covariates =) on
# lta(n_group_classes = 2). Data are drawn from a mixture of two Markov
# chains with known coefficients (a stable and a mobile group class), the
# step-one lta() is fitted without covariates, and step three's estimates are
# compared with the generating values.
#
# The registry (lta-three-step-groups-simulation-registry.csv) is written
# BEFORE any replicate runs and the script refuses to run if it changes.
# Criteria, fixed in advance (those of the verification article, Layer 4):
#   - availability: every replicate fits (a failure is counted);
#   - bias: |mean estimate - truth| <= 3.5 Monte Carlo SE (sd / sqrt(R));
#   - coverage: the share of 95% intervals containing the truth lies in the
#     central 99.9% binomial range for R draws at 0.95.
# Robust (the default) and observed intervals are both reported; the gate is
# applied to the robust ones. Step three treats step one as fixed, so a
# shortfall in coverage would be a property of the method (Bakk, Oberski and
# Vermunt 2014), reported as such, not hidden.
#
# Run from the project root:
#   Rscript validation/lta-three-step-groups-simulation.R
# LTA3G_SIM_N sets the replicate count (default 200); serial, one core.
# LTA3G_SIM_START (default 1) is the first replicate, so a run can be split
# over processes and the replicate files combined; seeds are per replicate.
# Outputs: lta-three-step-groups-simulation-{replicates,summary}.csv and .txt

devtools::load_all(".", quiet = TRUE)
out_dir <- "validation"
n_datasets <- as.integer(Sys.getenv("LTA3G_SIM_N", "200"))
# LTA3G_OBSERVED=0 skips the second, observed-information step three (it
# halves the run time; the gate is on the robust intervals either way).
with_observed <- !identical(Sys.getenv("LTA3G_OBSERVED", "1"), "0")

design <- list(J = 800L, T = 6L, shares = c(0.5, 0.5),
               means = rbind(c(-1.5, -1.3), c(1.5, 1.3)),
               # per group class: initial log odds of profile 1 (intercept, w)
               initial = rbind(c(1, 0.8), c(-1, 0.8)),
               # per group class and origin: log odds of moving (intercept, z)
               moves = list(rbind(c(-2.0, 0.8), c(-1.8, -0.6)),
                            rbind(c(-0.4, 0.8), c(-0.2, -0.6))),
               base_seed = 31000000L)
registry <- data.frame(
  cell = "groups_k2_h2", persons = design$J, occasions = design$T,
  shares = paste(design$shares, collapse = ","),
  means = paste(apply(design$means, 1L, paste, collapse = ","), collapse = " | "),
  initial = paste(apply(design$initial, 1L, paste, collapse = ","), collapse = " | "),
  moves = paste(vapply(design$moves, function(m) paste(apply(m, 1L, paste, collapse = ","),
                                                        collapse = " / "), ""),
                collapse = " | "),
  base_seed = design$base_seed, stringsAsFactors = FALSE)
registry_file <- file.path(out_dir, "lta-three-step-groups-simulation-registry.csv")
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- utils::read.csv(text = paste(utils::capture.output(
    utils::write.csv(registry, row.names = FALSE)), collapse = "\n"),
    stringsAsFactors = FALSE)
  stopifnot("registry changed after it was declared" = identical(stored, current))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}

draw <- function(seed) {
  set.seed(seed)
  J <- design$J
  g <- 1L + (stats::runif(J) > design$shares[1L])
  w <- stats::rnorm(J)
  z <- matrix(stats::rnorm(J * design$T), J)
  state <- matrix(NA_integer_, J, design$T)
  first <- design$initial[g, 1L] + design$initial[g, 2L] * w
  state[, 1L] <- ifelse(stats::runif(J) < stats::plogis(first), 1L, 2L)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(design$T)[-1L]) {  # simulation only
    from <- state[, t - 1L]
    coefficients <- t(vapply(seq_len(J), function(i) design$moves[[g[i]]][from[i], ],
                             numeric(2)))
    moved <- stats::runif(J) < stats::plogis(coefficients[, 1L] + coefficients[, 2L] * z[, t])
    state[, t] <- ifelse(moved, 3L - from, from)
  }
  frame <- data.frame(id = rep(seq_len(J), each = design$T),
                      time = rep(seq_len(design$T), J),
                      w = rep(w, each = design$T), z = as.vector(t(z)))
  s <- as.vector(t(state))
  frame$y1 <- stats::rnorm(length(s), design$means[s, 1L], 1)
  frame$y2 <- stats::rnorm(length(s), design$means[s, 2L], 1)
  frame
}

# The true coefficients in r3step()'s labels: group logit of class 1 against
# class 2; per class, initial log odds of profile 1 against 2 and, per
# origin, of moving rather than staying.
truth <- rbind(
  data.frame(key = "group 1 (Intercept)", truth = log(design$shares[1L] / design$shares[2L])),
  do.call(rbind, lapply(1:2, function(h) rbind(
    data.frame(key = paste("initial", h, c("(Intercept)", "w")), truth = design$initial[h, ]),
    data.frame(key = paste("transition", h, rep(c("1->2", "2->1"), each = 2L),
                           c("(Intercept)", "z")),
               truth = as.vector(t(design$moves[[h]])))))))

run_one <- function(replicate) {
  seed <- design$base_seed + replicate
  data <- draw(seed)
  conditions <- character()
  record <- function(condition) {
    conditions <<- c(conditions, class(condition)[1L])
    invokeRestart(if (inherits(condition, "warning")) "muffleWarning" else "muffleMessage")
  }
  fit <- withCallingHandlers(
    lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", n_group_classes = 2,
        n_starts = 5, seed = seed), warning = record, message = record)
  results <- lapply(c("robust", if (with_observed) "observed"), function(type) {
    withCallingHandlers(
      r3step(fit, data, initial_covariates = "w", transition_covariates = "z",
             vcov_type = type), warning = record, message = record)
  })
  # Label switching. Profile 1 has the lower first indicator; group class 1
  # is the stable one (lower mean move intercept).
  swap_profiles <- fit$means[1L, "y1"] > fit$means[2L, "y1"]
  r <- results[[1L]]
  move_intercepts <- vapply(c("group_class_1", "group_class_2"), function(h) {
    mean(r$estimate[r$block == "transition" & r$group_class == h & r$term == "(Intercept)"])
  }, numeric(1))
  swap_groups <- move_intercepts[[1L]] > move_intercepts[[2L]]
  keyed <- function(r) {
    profile <- function(p) {
      k <- as.integer(sub("profile_", "", p))
      if (swap_profiles) 3L - k else k
    }
    group <- function(h) {
      k <- as.integer(sub("group_class_", "", h))
      if (swap_groups) 3L - k else k
    }
    flip <- (r$block == "group" & swap_groups) | (r$block == "initial" & swap_profiles)
    key <- ifelse(r$block == "group", "group 1 (Intercept)",
                  ifelse(r$block == "initial", paste("initial", group(r$group_class), r$term),
                         paste("transition", group(r$group_class),
                               paste0(profile(r$from), "->", profile(r$to)), r$term)))
    data.frame(key = key, estimate = ifelse(flip, -r$estimate, r$estimate),
               se = r$standard_error,
               low = ifelse(flip, -r$conf_high, r$conf_low),
               high = ifelse(flip, -r$conf_low, r$conf_high), stringsAsFactors = FALSE)
  }
  robust <- keyed(results[[1L]])
  observed <- if (with_observed) keyed(results[[2L]]) else
    transform(robust, se = NA_real_, low = NA_real_, high = NA_real_)
  rows <- match(truth$key, robust$key)
  stopifnot("every true coefficient must be matched" = !anyNA(rows))
  data.frame(replicate = replicate, truth, estimate = robust$estimate[rows],
             se_robust = robust$se[rows], low_robust = robust$low[rows],
             high_robust = robust$high[rows],
             se_observed = observed$se[match(truth$key, observed$key)],
             low_observed = observed$low[match(truth$key, observed$key)],
             high_observed = observed$high[match(truth$key, observed$key)],
             step_one_converged = fit$converged,
             conditions = paste(unique(conditions), collapse = ";"),
             stringsAsFactors = FALSE)
}

first <- as.integer(Sys.getenv("LTA3G_SIM_START", "1"))
parts_dir <- file.path(out_dir, "lta-three-step-groups-simulation-parts")
started <- Sys.time()
if (n_datasets > 0L) {
  part <- do.call(rbind, lapply(first - 1L + seq_len(n_datasets), function(replicate) {
    tryCatch(run_one(replicate), error = function(error) {
      data.frame(replicate = replicate, key = NA_character_, truth = NA_real_,
                 estimate = NA_real_, se_robust = NA_real_, low_robust = NA_real_,
                 high_robust = NA_real_, se_observed = NA_real_, low_observed = NA_real_,
                 high_observed = NA_real_, step_one_converged = NA,
                 conditions = paste("error:", conditionMessage(error)),
                 stringsAsFactors = FALSE)
    })
  }))
  dir.create(parts_dir, showWarnings = FALSE)
  utils::write.csv(part, file.path(parts_dir, sprintf("replicates-%04d-%04d.csv", first,
                                                      first + n_datasets - 1L)),
                   row.names = FALSE)
}
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
# The summary reads every part written so far (LTA3G_SIM_N=0 only summarizes).
replicates <- do.call(rbind, lapply(list.files(parts_dir, "^replicates-.*[.]csv$",
                                               full.names = TRUE),
                                    utils::read.csv, stringsAsFactors = FALSE))
replicates <- replicates[order(replicates$replicate), ]
# A part with no condition reads the column back as logical NA.
replicates$conditions <- ifelse(is.na(replicates$conditions), "",
                                as.character(replicates$conditions))
n_datasets <- length(unique(replicates$replicate))
utils::write.csv(replicates, file.path(out_dir, "lta-three-step-groups-simulation-replicates.csv"),
                 row.names = FALSE)
usable <- replicates[!is.na(replicates$key), ]
fitted <- length(unique(usable$replicate))
coverage_range <- stats::qbinom(c(0.0005, 0.9995), fitted, 0.95) / fitted
summary_table <- do.call(rbind, lapply(truth$key, function(key) {
  rows <- usable[usable$key == key, ]
  mcse <- stats::sd(rows$estimate) / sqrt(nrow(rows))
  bias <- mean(rows$estimate) - rows$truth[1L]
  covered <- function(low, high) mean(rows$truth >= low & rows$truth <= high)
  data.frame(key = key, truth = rows$truth[1L], replicates = nrow(rows),
             mean_estimate = mean(rows$estimate), bias = bias, mcse = mcse,
             bias_over_mcse = bias / mcse,
             se_ratio_robust = mean(rows$se_robust) / stats::sd(rows$estimate),
             coverage_robust = covered(rows$low_robust, rows$high_robust),
             coverage_observed = covered(rows$low_observed, rows$high_observed),
             stringsAsFactors = FALSE)
}))
summary_table$bias_ok <- abs(summary_table$bias_over_mcse) <= 3.5
summary_table$coverage_ok <- summary_table$coverage_robust >= coverage_range[1L] &
  summary_table$coverage_robust <= coverage_range[2L]
utils::write.csv(summary_table, file.path(out_dir, "lta-three-step-groups-simulation-summary.csv"),
                 row.names = FALSE)
sink(file.path(out_dir, "lta-three-step-groups-simulation.txt"))
cat(sprintf(paste("Three-step latent Markov model with group classes: %d datasets",
                  "attempted, %d fitted; %s\n"), n_datasets, fitted, R.version.string))
cat(sprintf("Coverage range (central 99.9%% binomial, %d draws at 0.95): %.3f-%.3f\n\n",
            fitted, coverage_range[1L], coverage_range[2L]))
failures <- unique(replicates$conditions[startsWith(replicates$conditions, "error:")])
if (length(failures)) cat("Failed replicates:", failures, sep = "\n  ")
print(summary_table, digits = 3, row.names = FALSE)
cat(sprintf("\nBias gate met: %d of %d; coverage gate met: %d of %d\n",
            sum(summary_table$bias_ok), nrow(summary_table),
            sum(summary_table$coverage_ok), nrow(summary_table)))
sink()
cat(sprintf("done in %.1f minutes\n", elapsed))
