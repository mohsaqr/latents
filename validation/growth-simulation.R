# Predeclared simulation of the growth mixture model's estimates, standard
# errors and person-level R3STEP (phases G2 and G4 of todo.MD), extended on
# 2026-11-15 with correlated residuals (phase G5, `residual_correlation`) and
# the coverage of the person-bootstrap intervals (`parameter_inference(method
# = "bootstrap")`). The registry (growth-simulation-registry.csv) is written
# BEFORE any replicate of a cell runs, and the script refuses to run if a
# declared row changes; rows are only ever added.
#
# Three studies share the generator:
#   wald        the four cells of 2026-10-17 (balanced, unbalanced, weak,
#               small_variance): observed and robust Wald intervals, R3STEP.
#               Gates as in ADDITIVE_SIMULATION.md: availability >= 98%,
#               |bias| <= 0.1 empirical SD, observed SE ratio (mean SE /
#               empirical SD) 0.9-1.1, observed 95% coverage 0.925-0.975.
#   correlated  ar and brownian: the balanced design with a residual process
#               shared by the classes, fitted with the matching
#               `residual_correlation`; observed and robust Wald intervals.
#   bootstrap   bootstrap_balanced and bootstrap_weak: the balanced and weak
#               designs; percentile intervals of the person bootstrap beside
#               the Wald intervals of the same fits.
# The correlated and bootstrap studies are judged by the recovery criteria
# of the verification article (Layer 4), fixed before their first replicate:
# availability >= 98%; |mean estimate - truth| <= 3.5 Monte Carlo SE
# (3.5 sd / sqrt(R)); the share of 95% intervals (observed Wald for
# correlated, bootstrap percentile for bootstrap) covering the truth inside
# the central 99.9% binomial range for R draws at 0.95, R the replicates
# available. The wald study's gates are also reported for them.
#
# Class membership depends on a person covariate x, logit P(class 2) =
# 0.2 + 0.8 x. The growth mixture is fitted without it (its marginal class
# shares are then the truth, E[P(class | x)]), and r3step(fit, data, "x")
# estimates the slope (wald study only), whose truth is -0.8 or 0.8 for
# `class_1` against the reference class depending on which true class the
# fit calls class 1.
#
# Every replicate sets its own seeds (the data from base_seed + i, the fit
# and the bootstrap from i), so its result does not depend on the worker
# that runs it: the fork (mclapply) and socket (PSOCK) backends give the
# same replicates.
#
# Run from the project root:
#   Rscript validation/growth-simulation.R                       # every study
#   GROWTH_SIM_STUDY=correlated Rscript validation/growth-simulation.R
#   GROWTH_SIM_STUDY=bootstrap GROWTH_SIM_N=100 GROWTH_SIM_BACKEND=psock \
#     GROWTH_SIM_CORES=6 Rscript validation/growth-simulation.R  # reduced run
# GROWTH_SIM_N, when set, overrides the registry's n_datasets (a reduced run;
# the output records the count used); GROWTH_SIM_CELLS restricts a study to
# the cells named; GROWTH_SIM_DECLARE_ONLY=1 writes the registry and stops. Outputs: for the wald study
# growth-simulation-{replicates,summary,r3step}.csv and growth-simulation.txt;
# for the others growth-simulation-<study>-{replicates,summary}.csv and
# growth-simulation-<study>.txt.

source(file.path("validation", "fork-runtime.R"))
devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
out_dir <- Sys.getenv("GROWTH_SIM_OUT", "validation")
n_override <- Sys.getenv("GROWTH_SIM_N", "")
n_cores <- as.integer(Sys.getenv("GROWTH_SIM_CORES", "8"))
backend <- Sys.getenv("GROWTH_SIM_BACKEND", "fork")
studies <- strsplit(Sys.getenv("GROWTH_SIM_STUDY", "wald,correlated,bootstrap"), ",")[[1L]]
only_cells <- strsplit(Sys.getenv("GROWTH_SIM_CELLS", ""), ",")[[1L]]
stopifnot(backend %in% c("fork", "psock"),
          all(studies %in% c("wald", "correlated", "bootstrap")))

base_class <- list(intercept = c(2, 5), slope = c(0.6, -0.4),
                   variance = list(c(1, 0.09, 0.3), c(0.6, 0.05, -0.2)),
                   sigma = c(0.7, 0.8))
weak_class <- utils::modifyList(base_class, list(intercept = c(2, 3.2), slope = c(0.6, 0.1)))
# Replace whole elements of a class list (modifyList() would merge, and leave
# an unnamed list such as `variance` unchanged).
with_class <- function(base, ...) {
  value <- list(...)
  base[names(value)] <- value
  base
}
wald_cell <- function(...) c(list(study = "wald", inference = "wald", process = NULL,
                                  n_datasets = 1000L, random = "1 + time"), list(...))
cells <- list(
  balanced = wald_cell(J = 300L, times = 0:5, jitter = 0, missing = 0,
                       classes = base_class),
  unbalanced = wald_cell(J = 300L, times = 0:5, jitter = 0.3, missing = 0.25,
                         classes = base_class),
  weak = wald_cell(J = 500L, times = 0:5, jitter = 0, missing = 0, classes = weak_class),
  # Declared as "random-slope SD 0.1", but utils::modifyList() leaves an
  # unnamed list unchanged, so the cell as declared and run on 2026-10-17
  # drew the base variances: a second balanced cell with its own seeds (found
  # 2026-11-15; the registry records what ran, and so does this line).
  small_variance = wald_cell(J = 300L, times = 0:5, jitter = 0, missing = 0,
                             classes = base_class),
  # Correlated residuals: the balanced design over 8 occasions with a random
  # intercept, its residual variance split into a class-specific white noise
  # and a process shared by the classes. (A pilot of 6 replicates with a
  # random intercept and slope over 6 occasions left the autoregressive
  # process on its bound in 1 of 6 and its rate barely identified, standard
  # errors as large as the estimate: the process and the random slope compete
  # for the same covariance.)
  ar = list(study = "correlated", inference = "wald", n_datasets = 500L, random = "1",
            J = 300L, times = 0:7, jitter = 0, missing = 0,
            classes = with_class(base_class, sigma = c(0.5, 0.6), variance = list(1, 0.6)),
            process = list(type = "ar", rate = 0.7, sd = 0.6)),
  brownian = list(study = "correlated", inference = "wald", n_datasets = 500L, random = "1",
                  J = 300L, times = 0:7, jitter = 0, missing = 0,
                  classes = with_class(base_class, sigma = c(0.5, 0.6), variance = list(1, 0.6)),
                  process = list(type = "brownian", sd = 0.3)),
  # Bootstrap intervals: percentile intervals of 49 person resamples, each
  # refitted from 1 start (the cost: one replicate is 50 fits).
  bootstrap_balanced = list(study = "bootstrap", inference = "bootstrap", n_datasets = 200L,
                            iter = 49L, boot_starts = 1L,
                            J = 300L, times = 0:5, jitter = 0, missing = 0,
                            classes = base_class, process = NULL, random = "1 + time"),
  bootstrap_weak = list(study = "bootstrap", inference = "bootstrap", n_datasets = 200L,
                        iter = 49L, boot_starts = 1L,
                        J = 500L, times = 0:5, jitter = 0, missing = 0,
                        classes = weak_class, process = NULL, random = "1 + time"))
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
             study = cell$study, inference = cell$inference, random = cell$random,
             residual_correlation = cell$process$type %||% "none",
             process = if (is.null(cell$process)) "" else
               paste(names(cell$process)[-1L], unlist(cell$process[-1L]), sep = "=",
                     collapse = ","),
             bootstrap_iter = cell$iter %||% NA_integer_,
             bootstrap_starts = cell$boot_starts %||% NA_integer_,
             n_datasets = cell$n_datasets,
             stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "growth-simulation-registry.csv")
round_trip <- function(x) utils::read.csv(text = paste(utils::capture.output(
  utils::write.csv(x, row.names = FALSE)), collapse = "\n"), stringsAsFactors = FALSE)
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- round_trip(registry)
  # Declared rows never change; a row is added (and declared) before its
  # cell's first replicate.
  stopifnot("registry changed after it was declared" =
              identical(stored, current[seq_len(nrow(stored)), names(stored)]))
  if (nrow(current) > nrow(stored) || !identical(names(stored), names(current))) {
    utils::write.csv(registry, registry_file, row.names = FALSE)
  }
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}
if (nzchar(Sys.getenv("GROWTH_SIM_DECLARE_ONLY"))) quit(save = "no")

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
    v <- cell$classes$variance[[class[i]]]
    # A random intercept alone: its variance, and no slope.
    if (length(v) == 1L) return(c(sqrt(v) * stats::rnorm(1L), 0))
    as.vector(t(chol(covariance_of(v))) %*% stats::rnorm(2L))
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
  if (!is.null(cell$process)) {
    # The shared residual process, drawn after the cells of 2026-10-17 drew
    # everything else, from its covariance at each person's times (the
    # parameterization of lcmm::hlme(cor = AR(time)) and BM(time)).
    p <- cell$process
    frame$y <- frame$y + unlist(lapply(split(frame$time, frame$id), function(t) {
      w <- if (identical(p$type, "ar")) p$sd^2 * exp(-p$rate * abs(outer(t, t, "-"))) else
        p$sd^2 * outer(t, t, pmin)
      if (identical(p$type, "brownian") && t[1L] == 0) {
        # Brownian motion is zero at time 0.
        c(0, as.vector(t(chol(w[-1L, -1L, drop = FALSE])) %*% stats::rnorm(length(t) - 1L)))
      } else as.vector(t(chol(w)) %*% stats::rnorm(length(t)))
    }), use.names = FALSE)
  }
  # Missing occasions at random; every person keeps their first two.
  keep <- stats::runif(nrow(frame)) >= cell$missing |
    stats::ave(frame$time, frame$id, FUN = seq_along) <= 2
  frame[keep, ]
}

truth_rows <- function(cell) {
  share_2 <- stats::integrate(function(x) stats::plogis(0.2 + 0.8 * x) * stats::dnorm(x),
                              -Inf, Inf)$value
  rows <- do.call(rbind, lapply(1:2, function(k) {
    v <- cell$classes$variance[[k]]
    label <- paste0("class_", k)
    rows <- data.frame(
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
    # A random intercept alone has no slope variance and no correlation.
    random_slope <- startsWith(rows$block, "random_") & rows$term != "(Intercept)"
    if (length(v) == 1L) rows <- rows[!random_slope, ]
    rows
  }))
  if (!is.null(cell$process)) {
    values <- unlist(cell$process[-1L])
    rows <- rbind(rows, data.frame(block = "residual_process", class = "all",
                                   term = names(values), truth = unname(values),
                                   stringsAsFactors = FALSE))
  }
  rows
}

one_replicate <- function(cell_name, i) {
  cell <- cells[[cell_name]]
  seed <- registry$base_seed[registry$cell == cell_name] + i
  data <- draw(cell, seed)
  truth <- truth_rows(cell)
  extended <- !identical(cell$study, "wald")
  failure <- function(reason) {
    rows <- data.frame(cell = cell_name, replicate = i, truth[c("block", "class", "term")],
                       truth = truth$truth, estimate = NA_real_, se_observed = NA_real_,
                       se_robust = NA_real_, covered_observed = NA, covered_robust = NA,
                       available = FALSE, reason = reason, stringsAsFactors = FALSE)
    if (extended) {
      rows$se_bootstrap <- NA_real_
      rows$covered_bootstrap <- NA
    }
    rows
  }
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                       random = stats::as.formula(paste("~", cell$random)), n_starts = 5,
                       seed = i, vcov_type = "none",
                       residual_correlation = cell$process$type %||% "none"),
    warning = function(w) {
      warnings <<- c(warnings, class(w)[1L])
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) return(list(rows = failure(conditionMessage(fit)), r3 = NULL))
  if (!isTRUE(fit$converged)) return(list(rows = failure("unconverged"), r3 = NULL))
  if (length(fit$degenerate_classes %||% character()) > 0L) {
    return(list(rows = failure("degenerate covariance"), r3 = NULL))
  }
  # `residual_boundary` names a process coordinate or a class's white-noise
  # SD held at its floor; the reason below covers both (the 2026-11-15 run's
  # autoregressive misses were all a residual SD at its floor).
  if (length(fit$residual_boundary %||% character()) > 0L) {
    return(list(rows = failure("residual process on its bound"), r3 = NULL))
  }
  # Name the fitted classes by the true class they estimate.
  signature <- function(intercepts, slopes) cbind(intercepts, slopes)
  order <- .multilpa_match_order(
    signature(cell$classes$intercept, cell$classes$slope),
    signature(fit$params$beta[1L, ], fit$params$beta[2L, ]))
  true_label <- paste0("class_", match(seq_len(2L), order))
  tables <- tryCatch(c(list(
    observed = parameter_inference(fit, vcov_type = "observed"),
    robust = parameter_inference(fit, vcov_type = "robust")),
    if (identical(cell$inference, "bootstrap")) list(
      bootstrap = suppressWarnings(parameter_inference(
        fit, method = "bootstrap", iter = cell$iter, n_starts = cell$boot_starts,
        seed = i), classes = "latents_bootstrap_dropped"))),
    error = function(e) e)
  if (inherits(tables, "error")) return(list(rows = failure(conditionMessage(tables)), r3 = NULL))
  relabel <- function(table) {
    table <- subset(table, block != "membership")
    table$class <- ifelse(table$class == "all", "all",
                          true_label[match(table$class, c("class_1", "class_2"))])
    merge(truth, table, by = c("block", "class", "term"), sort = FALSE)
  }
  observed <- relabel(tables$observed)
  key <- paste(observed$block, observed$class, observed$term)
  aligned <- function(table) {
    table <- relabel(table)
    table[match(key, paste(table$block, table$class, table$term)), ]
  }
  robust <- aligned(tables$robust)
  covered <- function(table) table$conf_low <= table$truth & table$truth <= table$conf_high
  rows <- data.frame(cell = cell_name, replicate = i,
                     observed[c("block", "class", "term", "truth", "estimate")],
                     se_observed = observed$standard_error,
                     se_robust = robust$standard_error,
                     covered_observed = covered(observed), covered_robust = covered(robust),
                     available = is.finite(observed$standard_error),
                     reason = NA_character_, stringsAsFactors = FALSE)
  if (extended) {
    rows$se_bootstrap <- NA_real_
    rows$covered_bootstrap <- NA
    if (!is.null(tables$bootstrap)) {
      bootstrap <- aligned(tables$bootstrap)
      rows$se_bootstrap <- bootstrap$standard_error
      rows$covered_bootstrap <- covered(bootstrap)
      rows$available <- rows$available & is.finite(bootstrap$conf_low)
    }
    return(list(rows = rows, r3 = NULL))
  }
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

run_replicates <- function(jobs) {
  if (identical(backend, "psock") && n_cores > 1L) {
    cluster <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(root) {
      setwd(root)
      suppressMessages(devtools::load_all(".", quiet = TRUE))
      RNGkind("L'Ecuyer-CMRG")
      NULL
    }, getwd())
    parallel::clusterExport(cluster, c("cells", "registry", "draw", "truth_rows",
                                       "covariance_of", "one_replicate"),
                            envir = globalenv())
    parallel::parLapplyLB(cluster, jobs, function(job) {
      try(one_replicate(job$cell, job$i), silent = TRUE)
    })
  } else {
    parallel::mclapply(jobs, function(job) one_replicate(job$cell, job$i),
                       mc.cores = .validation_fork_cores(n_cores), mc.set.seed = TRUE)
  }
}

# Recovery criteria of the verification article (Layer 4).
binomial_range <- function(n, p = 0.95, level = 0.999) {
  stats::qbinom(c((1 - level) / 2, (1 + level) / 2), n, p) / n
}

summarise <- function(rows, methods) {
  ok <- rows[rows$available, ]
  sd_empirical <- stats::sd(ok$estimate)
  bias <- mean(ok$estimate) - rows$truth[1L]
  out <- data.frame(availability = mean(rows$available), n_available = nrow(ok),
                    bias = bias, bias_sd = bias / sd_empirical,
                    bias_mcse = bias / (sd_empirical / sqrt(nrow(ok))),
                    empirical_sd = sd_empirical)
  for (m in methods) out[[paste0("se_ratio_", m)]] <- mean(ok[[paste0("se_", m)]]) / sd_empirical
  for (m in methods) out[[paste0("coverage_", m)]] <- mean(ok[[paste0("covered_", m)]])
  out
}

run_study <- function(study) {
  names_in <- names(cells)[vapply(cells, `[[`, character(1), "study") == study]
  if (length(only_cells) > 0L) names_in <- intersect(names_in, only_cells)
  n_of <- function(cell_name) {
    if (nzchar(n_override)) as.integer(n_override) else cells[[cell_name]]$n_datasets
  }
  jobs <- unlist(lapply(names_in, function(cell_name) {
    lapply(seq_len(n_of(cell_name)), function(i) list(cell = cell_name, i = i))
  }), recursive = FALSE)
  started <- Sys.time()
  results <- run_replicates(jobs)
  errors <- vapply(results, inherits, logical(1), what = "try-error")
  stopifnot("a worker failed" = !any(errors))
  replicates <- do.call(rbind, lapply(results, `[[`, "rows"))
  prefix <- if (identical(study, "wald")) "growth-simulation" else
    paste0("growth-simulation-", study)
  utils::write.csv(replicates, file.path(out_dir, paste0(prefix, "-replicates.csv")),
                   row.names = FALSE)
  methods <- c("observed", "robust", if (identical(study, "bootstrap")) "bootstrap")
  keys <- unique(replicates[c("cell", "block", "class", "term")])
  summary <- do.call(rbind, lapply(seq_len(nrow(keys)), function(r) {
    rows <- merge(keys[r, ], replicates, by = c("cell", "block", "class", "term"))
    cbind(keys[r, ], truth = rows$truth[1L], summarise(rows, methods))
  }))
  gates <- summary$availability >= 0.98 & abs(summary$bias_sd) <= 0.1 &
    summary$se_ratio_observed >= 0.9 & summary$se_ratio_observed <= 1.1 &
    summary$coverage_observed >= 0.925 & summary$coverage_observed <= 0.975
  if (identical(study, "wald")) {
    # The columns of the 2026-10-17 run.
    summary <- summary[setdiff(names(summary), c("n_available", "bias_mcse"))]
    summary$pass <- gates
  } else {
    primary <- if (identical(study, "bootstrap")) "coverage_bootstrap" else "coverage_observed"
    range <- t(vapply(summary$n_available, binomial_range, numeric(2)))
    summary$coverage_low <- range[, 1L]
    summary$coverage_high <- range[, 2L]
    summary$criterion <- primary
    summary$pass_gates <- gates
    summary$pass <- summary$availability >= 0.98 & abs(summary$bias_mcse) <= 3.5 &
      summary[[primary]] >= summary$coverage_low & summary[[primary]] <= summary$coverage_high
  }
  rownames(summary) <- NULL
  utils::write.csv(summary, file.path(out_dir, paste0(prefix, "-summary.csv")),
                   row.names = FALSE)
  reasons <- table(replicates$reason[!duplicated(replicates[c("cell", "replicate")])],
                   useNA = "no")
  r3_summary <- NULL
  if (identical(study, "wald")) {
    r3 <- do.call(rbind, lapply(results, `[[`, "r3"))
    utils::write.csv(r3, file.path(out_dir, "growth-simulation-r3step.csv"), row.names = FALSE)
    r3_summary <- do.call(rbind, lapply(names_in, function(cell_name) {
      rows <- r3[r3$cell == cell_name & is.finite(r3$estimate), ]
      aligned <- rows$estimate * sign(rows$truth)
      data.frame(cell = cell_name, available = nrow(rows) / n_of(cell_name),
                 bias_sd = (mean(aligned) - 0.8) / stats::sd(aligned),
                 se_ratio = mean(rows$se) / stats::sd(aligned),
                 coverage = mean(rows$covered))
    }))
  }
  report <- file.path(out_dir, paste0(prefix, ".txt"))
  sink(report)
  cat(sprintf("Growth mixture simulation, study %s: %s datasets per cell (%s), %s elapsed, %s backend, %d cores\n\n",
              study, paste(vapply(names_in, n_of, integer(1)), collapse = "/"),
              paste(names_in, collapse = "/"), format(round(Sys.time() - started, 1)),
              backend, n_cores))
  print(summary, digits = 3, row.names = FALSE)
  if (!is.null(r3_summary)) {
    cat("\nR3STEP slope on x (true |slope| 0.8, robust intervals)\n")
    print(r3_summary, digits = 3, row.names = FALSE)
  }
  cat("\nUnavailable replicates by reason\n")
  print(reasons)
  cat("\n")
  print(utils::sessionInfo())
  sink()
  cat(readLines(report), sep = "\n")
}

for (study in studies) run_study(study)
