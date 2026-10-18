# Starts study: how many random starts reach the best maximum, at the true
# number of classes and above it, and does the number of starts change which
# number of classes BIC selects?
#
# Design (fully crossed): family (LPA VVI on 4 indicators, LCA on 6 binary
# items) x true K (2, 3, 4) x separation (low, moderate, high) x n (300, 1000)
# x 6 replicates = 216 datasets. Each dataset is fitted at k = max(1, K - 1)
# to K + 2 with 500 starts.
#
# Why one 500-start fit per (dataset, k) is enough: start 1 is the
# deterministic Ward start and starts 2..500 are independent k-means starts
# drawn in sequence (R/multilpa.R, .multilpa_initialize; R/kernel-em.R,
# .latents_run_starts). A fit with s starts is the best of the Ward start and
# s - 1 independent k-means starts, so disjoint blocks of the 499 k-means
# starts, each joined with the Ward start, are independent draws of an
# s-start fit. check-emulation.R confirms this against real 10- and 20-start
# fits.
#
# Run from the project root, one shard per process:
#   Rscript equivalence/starts-study/simulate.R <shard> <n_shards>

starts_design <- function(reps = 6L) {
  grid <- expand.grid(family = c("lpa", "lca"), true_k = 2:4,
                      separation = c("low", "moderate", "high"), n = c(300L, 1000L),
                      rep = seq_len(reps), stringsAsFactors = FALSE)
  grid$dataset <- sprintf("%s_k%d_%s_n%d_r%d", grid$family, grid$true_k,
                          grid$separation, grid$n, grid$rep)
  grid$data_seed <- 20261018L + seq_len(nrow(grid))
  grid
}

# Distinct 0/1 profile patterns: the first K rows, so each class differs from
# every other on at least one indicator.
starts_patterns <- function(k, p) {
  patterns <- list(
    lpa = rbind(c(0, 0, 0, 0), c(1, 1, 0, 0), c(0, 1, 1, 1), c(1, 0, 1, 0)),
    lca = rbind(c(0, 0, 0, 0, 0, 0), c(1, 1, 1, 0, 0, 0), c(0, 0, 1, 1, 1, 1),
                c(1, 1, 1, 1, 1, 1))
  )
  patterns[[if (p == 4L) "lpa" else "lca"]][seq_len(k), , drop = FALSE]
}

starts_proportions <- function(k) {
  list(c(0.6, 0.4), c(0.5, 0.3, 0.2), c(0.4, 0.3, 0.2, 0.1))[[k - 1L]]
}

starts_simulate <- function(row) {
  set.seed(row$data_seed)
  k <- row$true_k
  n <- row$n
  class <- sample.int(k, n, replace = TRUE, prob = starts_proportions(k))
  if (identical(row$family, "lpa")) {
    shift <- c(low = 1.5, moderate = 2.5, high = 3.5)[[row$separation]]
    means <- starts_patterns(k, 4L) * shift
    sds <- matrix(sqrt(stats::runif(k * 4L, 0.7, 1.3)), k, 4L)
    values <- means[class, ] + matrix(stats::rnorm(n * 4L), n, 4L) * sds[class, ]
    data <- as.data.frame(values)
    names(data) <- paste0("y", 1:4)
  } else {
    high <- c(low = 0.7, moderate = 0.8, high = 0.9)[[row$separation]]
    probabilities <- ifelse(starts_patterns(k, 6L) == 1, high, 1 - high)
    values <- matrix(stats::rbinom(n * 6L, 1L, probabilities[class, ]), n, 6L)
    data <- as.data.frame(values)
    names(data) <- paste0("u", 1:6)
  }
  data$id <- seq_len(n)
  list(data = data, class = class)
}

starts_fit <- function(data, family, k, n_starts, seed) {
  vars <- setdiff(names(data), "id")
  warnings <- character()
  fit <- withCallingHandlers(
    if (identical(family, "lpa")) {
      lpa(data, vars, n_profiles = k, model = "VVI", n_starts = n_starts, seed = seed)
    } else {
      lca(data, vars, n_classes = k, n_starts = n_starts, seed = seed)
    },
    warning = function(w) {
      warnings <<- c(warnings, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = paste(unique(warnings), collapse = ";"))
}

starts_one <- function(row, n_starts = 500L) {
  simulated <- starts_simulate(row)
  ks <- seq(max(1L, row$true_k - 1L), row$true_k + 2L)
  do.call(rbind, lapply(ks, function(k) {
    timing <- system.time(result <- starts_fit(simulated$data, row$family, k,
                                               n_starts, seed = row$data_seed))
    fit <- result$fit
    starts <- fit$starts
    data.frame(dataset = row$dataset, family = row$family, true_k = row$true_k,
               separation = row$separation, n = row$n, rep = row$rep, k = k,
               start = starts$start, log_likelihood = starts$log_likelihood,
               converged = starts$converged, boundary = starts$boundary,
               final_log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters,
               smallest_class = min(tabulate(fit$subject_profiles, k)),
               entropy = fit$entropy %||% NA_real_, seconds = timing[["elapsed"]],
               warnings = result$warnings, stringsAsFactors = FALSE)
  }))
}

if (sys.nframe() == 0L) {
  arguments <- as.integer(commandArgs(trailingOnly = TRUE))
  shard <- arguments[1L]
  n_shards <- arguments[2L]
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  mine <- subset(starts_design(), (data_seed - 20261018L) %% n_shards == shard %% n_shards)
  out_dir <- file.path("equivalence", "starts-study", "fits")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  invisible(lapply(seq_len(nrow(mine)), function(i) {
    path <- file.path(out_dir, paste0(mine$dataset[i], ".rds"))
    if (file.exists(path)) return(NULL)
    saveRDS(starts_one(mine[i, ]), path)
    message(sprintf("shard %d: %s done", shard, mine$dataset[i]))
  }))
}
