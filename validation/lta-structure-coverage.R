# Predeclared coverage study of the Wald intervals for the constrained
# covariance structures in the homogeneous lta() (`multilpa_transitions`).
# The registry (lta-structure-coverage-registry.csv) is written BEFORE any
# replicate runs and the script refuses to run if it changes.
#
# Gates, fixed with the registry (as in LTA_SIMULATION.md, with the coverage
# gate of the verification article): availability >= 98%; |bias| <= 0.1
# empirical SD; observed SE ratio (mean SE / empirical SD) 0.9-1.1; observed
# 95% coverage inside the central 99.9% binomial range for the number of
# usable replicates.
#
# Each replicate draws its data from its own seed (base_seed + replicate), so
# the result does not depend on the number of workers or on the platform's
# parallel back end. Run from the project root:
#
#   Rscript validation/lta-structure-coverage.R                 # 1000 per cell
#   LTA_STRUCT_N=200 LTA_STRUCT_CORES=12 Rscript validation/lta-structure-coverage.R
#
# Outputs: lta-structure-coverage-{replicates,summary}.csv and
# lta-structure-coverage.txt in validation/.

out_dir <- "validation"
n_datasets <- as.integer(Sys.getenv("LTA_STRUCT_N", "1000"))
cores <- as.integer(Sys.getenv("LTA_STRUCT_CORES", "1"))

rotation <- function(degrees) {
  r <- degrees * base::pi / 180
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2L)
}
cells <- list(
  EII = list(sigma = list(diag(2), diag(2))),
  VEI = list(sigma = list(0.8 * diag(c(1.5, 1 / 1.5)), 1.6 * diag(c(1.5, 1 / 1.5)))),
  EVE = list(sigma = list(rotation(30) %*% diag(c(2, 0.5)) %*% t(rotation(30)),
                          rotation(30) %*% diag(c(0.6, 1 / 0.6)) %*% t(rotation(30)))),
  VEE = list(sigma = list(0.7 * matrix(c(1.25, 0.5, 0.5, 1), 2L),
                          1.4 * matrix(c(1.25, 0.5, 0.5, 1), 2L))))
design <- list(groups = 250L, occasions = 5L, means = rbind(c(0, 0), c(2, 1.6)),
               initial = c(0.6, 0.4),
               transition = matrix(c(0.85, 0.15, 0.25, 0.75), 2L, byrow = TRUE))
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  sigma <- cells[[k]]$sigma
  data.frame(cell = names(cells)[k], groups = design$groups,
             occasions = design$occasions,
             means = "0,0 | 2,1.6", initial = "0.6,0.4", transition = "0.85,0.15 | 0.25,0.75",
             sigma_1 = paste(signif(as.vector(sigma[[1L]]), 10), collapse = ","),
             sigma_2 = paste(signif(as.vector(sigma[[2L]]), 10), collapse = ","),
             base_seed = 31e6L + 1e6L * k, n_starts = 3L,
             gates = "availability>=0.98; |bias|<=0.1 SD; SE ratio 0.9-1.1; coverage in central 99.9% binomial range",
             stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "lta-structure-coverage-registry.csv")
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- utils::read.csv(text = paste(utils::capture.output(
    utils::write.csv(registry, row.names = FALSE)), collapse = "\n"),
    stringsAsFactors = FALSE)
  stopifnot("registry changed after it was declared" = isTRUE(all.equal(stored, current)))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}

# The true natural parameters tracked, keyed as parameter_inference() labels
# them once profile 1 is the profile with the lower mean of `a`.
truth_of <- function(structure) {
  sigma <- cells[[structure]]$sigma
  spread <- if (endsWith(structure, "I")) {
    rbind(data.frame(parameter = "variance", outcome = "profile_1", term = "a",
                     truth = sigma[[1L]][1L, 1L]),
          data.frame(parameter = "variance", outcome = "profile_2", term = "b",
                     truth = sigma[[2L]][2L, 2L]))
  } else {
    rbind(data.frame(parameter = "covariance", outcome = "profile_1", term = "a:a",
                     truth = sigma[[1L]][1L, 1L]),
          data.frame(parameter = "covariance", outcome = "profile_1", term = "b:a",
                     truth = sigma[[1L]][2L, 1L]),
          data.frame(parameter = "covariance", outcome = "profile_2", term = "b:b",
                     truth = sigma[[2L]][2L, 2L]))
  }
  if (identical(structure, "EII")) {
    spread <- data.frame(parameter = "variance", outcome = "shared", term = "a", truth = 1)
  }
  rbind(
    data.frame(parameter = "initial_probability", outcome = "profile_1",
               term = "group_class_1", truth = design$initial[1L]),
    data.frame(parameter = "transition_probability", outcome = c("profile_1", "profile_2"),
               term = c("group_class_1:profile_1", "group_class_1:profile_2"),
               truth = diag(design$transition)),
    data.frame(parameter = "mean", outcome = c("profile_1", "profile_2"), term = "a",
               truth = design$means[, 1L]),
    spread)
}

draw <- function(structure, seed) {
  set.seed(seed)
  n <- design$groups
  waves <- design$occasions
  state <- matrix(NA_integer_, n, waves)
  state[, 1L] <- sample.int(2L, n, replace = TRUE, prob = design$initial)
  for (t in seq_len(waves)[-1L]) {  # simulation only: each state depends on the last
    state[, t] <- vapply(state[, t - 1L], function(previous) {
      sample.int(2L, 1L, prob = design$transition[previous, ])
    }, integer(1))
  }
  s <- as.vector(t(state))
  factors <- lapply(cells[[structure]]$sigma, chol)
  noise <- matrix(stats::rnorm(2L * length(s)), length(s), 2L)
  values <- t(vapply(seq_along(s), function(row) {
    design$means[s[row], ] + as.vector(noise[row, ] %*% factors[[s[row]]])
  }, numeric(2)))
  data.frame(id = rep(seq_len(n), each = waves), time = rep(seq_len(waves), n),
             a = values[, 1L], b = values[, 2L])
}

run_one <- function(structure, replicate) {
  seed <- registry$base_seed[registry$cell == structure] + replicate
  data <- draw(structure, seed)
  conditions <- character()
  record <- function(condition) {
    conditions <<- c(conditions, class(condition)[1L])
    invokeRestart(if (inherits(condition, "warning")) "muffleWarning" else "muffleMessage")
  }
  truth <- truth_of(structure)
  fit <- withCallingHandlers(
    lta(data, c("a", "b"), "id", 2L, time = "time", model = structure,
        n_starts = 3L, seed = seed, tol = 1e-10, max_iter = 20000L),
    warning = record, message = record)
  table <- parameter_inference(fit)
  # Label switching: profile 1 has the lower mean of `a`.
  swap <- fit$means[1L, "a"] > fit$means[2L, "a"]
  relabel <- function(text) {
    if (!swap) return(text)
    text <- gsub("profile_1", "profile_X", text, fixed = TRUE)
    text <- gsub("profile_2", "profile_1", text, fixed = TRUE)
    gsub("profile_X", "profile_2", text, fixed = TRUE)
  }
  key <- paste(table$parameter, relabel(table$outcome), relabel(table$term))
  rows <- match(paste(truth$parameter, truth$outcome, truth$term), key)
  stopifnot("every true parameter must be matched" = !anyNA(rows))
  data.frame(cell = structure, replicate = replicate, truth,
             estimate = table$estimate[rows], se = table$standard_error[rows],
             low = table$conf_low[rows], high = table$conf_high[rows],
             converged = fit$converged, boundary = isTRUE(fit$boundary),
             conditions = paste(unique(conditions), collapse = ";"),
             stringsAsFactors = FALSE)
}

run_safely <- function(job) {
  tryCatch(run_one(job$cell, job$replicate), error = function(error) {
    data.frame(cell = job$cell, replicate = job$replicate, parameter = NA_character_,
               outcome = NA_character_, term = NA_character_, truth = NA_real_,
               estimate = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_,
               converged = NA, boundary = NA,
               conditions = paste("error:", conditionMessage(error)),
               stringsAsFactors = FALSE)
  })
}

suppressMessages(pkgload::load_all(".", quiet = TRUE))
jobs <- expand.grid(replicate = seq_len(n_datasets), cell = names(cells),
                    stringsAsFactors = FALSE)
jobs <- lapply(seq_len(nrow(jobs)), function(k) jobs[k, ])
started <- Sys.time()
replicates_file <- file.path(out_dir, "lta-structure-coverage-replicates.csv")
# LTA_STRUCT_FROM_REPLICATES=1 summarises the retained replicates of a run
# with the same N instead of refitting them.
if (identical(Sys.getenv("LTA_STRUCT_FROM_REPLICATES"), "1")) {
  replicates <- utils::read.csv(replicates_file, stringsAsFactors = FALSE)
  stopifnot("the retained replicates must be of this N" =
              all(table(replicates$cell[!duplicated(replicates[c("cell", "replicate")])]) ==
                    n_datasets))
} else if (cores > 1L) {
  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(root) {
    suppressMessages(pkgload::load_all(root, quiet = TRUE))
    NULL
  }, getwd())
  parallel::clusterExport(cluster, c("cells", "design", "registry", "truth_of", "draw",
                                     "run_one"))
  replicates <- do.call(rbind, parallel::parLapplyLB(cluster, jobs, run_safely))
} else {
  replicates <- do.call(rbind, lapply(jobs, run_safely))
}
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
if (!identical(Sys.getenv("LTA_STRUCT_FROM_REPLICATES"), "1")) {
  utils::write.csv(replicates, replicates_file, row.names = FALSE)
}

attempted <- table(factor(vapply(jobs, `[[`, character(1), "cell"), names(cells)))
usable <- subset(replicates, !is.na(parameter) & converged & !boundary)
keys <- unique(usable[c("cell", "parameter", "outcome", "term")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  rows <- usable[usable$cell == keys$cell[k] & usable$parameter == keys$parameter[k] &
                   usable$outcome == keys$outcome[k] & usable$term == keys$term[k], ]
  n <- nrow(rows)
  covered <- mean(rows$low <= rows$truth & rows$truth <= rows$high)
  range <- stats::qbinom(c(0.0005, 0.9995), n, 0.95) / n
  empirical_sd <- stats::sd(rows$estimate)
  bias_sd <- (mean(rows$estimate) - rows$truth[1L]) / empirical_sd
  ratio <- mean(rows$se) / empirical_sd
  availability <- n / attempted[[keys$cell[k]]]
  data.frame(keys[k, ], truth = rows$truth[1L], usable = n, availability = availability,
             bias_sd = bias_sd, se_ratio = ratio, coverage = covered,
             coverage_low = range[1L], coverage_high = range[2L],
             pass = availability >= 0.98 && abs(bias_sd) <= 0.1 && ratio >= 0.9 &&
               ratio <= 1.1 && covered >= range[1L] && covered <= range[2L],
             row.names = NULL, stringsAsFactors = FALSE)
}))
utils::write.csv(summary_table, file.path(out_dir, "lta-structure-coverage-summary.csv"),
                 row.names = FALSE)
report <- c(
  if (identical(Sys.getenv("LTA_STRUCT_FROM_REPLICATES"), "1")) {
    sprintf("lta() constrained-structure coverage: %d datasets per cell, summarised from the retained replicates",
            n_datasets)
  } else sprintf("lta() constrained-structure coverage: %d datasets per cell, %d workers, %.1f minutes",
                 n_datasets, cores, elapsed),
  sprintf("latents %s, %s", as.character(utils::packageVersion("latents")), R.version.string),
  sprintf("errors: %d; unconverged or boundary: %d",
          sum(is.na(replicates$parameter)),
          length(unique(paste(replicates$cell, replicates$replicate)[
            !is.na(replicates$parameter) & !(replicates$converged & !replicates$boundary)]))),
  sprintf("%d of %d parameters pass every gate", sum(summary_table$pass),
          nrow(summary_table)),
  "",
  utils::capture.output(print(summary_table, digits = 3, row.names = FALSE)))
writeLines(report, file.path(out_dir, "lta-structure-coverage.txt"))
cat(report, sep = "\n")
