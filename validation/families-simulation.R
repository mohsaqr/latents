# Predeclared simulation of the dispersion and additive-dispersion families.
# Same procedure and release gates as additive-simulation.R; its own registry
# (families-simulation-registry.csv) is written BEFORE any replicate runs and
# the script refuses to run if it changes. Run from the project root:
#   Rscript validation/families-simulation.R
# Outputs: families-simulation-{replicates,summary}.csv, families-simulation.txt

devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
out_dir <- "validation"
n_datasets <- 1000L

two <- function(a, b) rbind(a, b)
cells <- list(
  dispersion_reference = list(family = "dispersion", between_variance = "equal",
    H = 2L, J = 200L, sizes = 10L, weights = c(0.5, 0.5),
    means = two(c(0, 0), c(0, 0)), between = two(c(0.4, 0.4), c(0.4, 0.4)),
    within = two(c(0.5, 0.5), c(2, 2))),
  dispersion_unequal_sizes = list(family = "dispersion", between_variance = "equal",
    H = 2L, J = 200L, sizes = c(2L, 5L, 10L, 20L), weights = c(0.5, 0.5),
    means = two(c(0, 0), c(0, 0)), between = two(c(0.4, 0.4), c(0.4, 0.4)),
    within = two(c(0.5, 0.5), c(2, 2))),
  additive_dispersion_varying = list(family = "additive_dispersion",
    between_variance = "varying", H = 2L, J = 200L, sizes = 10L,
    weights = c(0.5, 0.5), means = two(c(-1, -1), c(1, 1)),
    between = two(c(0.25, 0.5), c(0.5, 0.25)), within = two(c(0.6, 0.8), c(1.5, 1.2))),
  additive_dispersion_equal = list(family = "additive_dispersion",
    between_variance = "equal", H = 2L, J = 200L, sizes = 10L,
    weights = c(0.5, 0.5), means = two(c(-1, -1), c(1, 1)),
    between = two(c(0.4, 0.4), c(0.4, 0.4)), within = two(c(0.6, 0.8), c(1.5, 1.2))))
flat <- function(m) paste(apply(m, 1L, paste, collapse = ","), collapse = " | ")
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  data.frame(cell = names(cells)[k], family = cell$family,
             between_variance = cell$between_variance, n_groups = cell$J,
             sizes = paste(cell$sizes, collapse = "/"),
             weights = paste(cell$weights, collapse = "/"), means = flat(cell$means),
             between = flat(cell$between), within = flat(cell$within),
             base_seed = 5e6L + 1e6L * k, stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "families-simulation-registry.csv")
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
  sizes <- rep_len(cell$sizes, cell$J)
  classes <- sample.int(cell$H, cell$J, replace = TRUE, prob = cell$weights)
  do.call(rbind, lapply(seq_len(cell$J), function(j) {
    h <- classes[j]
    intercept <- stats::rnorm(2L, cell$means[h, ], sqrt(cell$between[h, ]))
    ratings <- sweep(matrix(stats::rnorm(sizes[j] * 2L), sizes[j], 2L) %*%
                       diag(sqrt(cell$within[h, ])), 2L, intercept, "+")
    data.frame(group = j, y1 = ratings[, 1L], y2 = ratings[, 2L])
  }))
}

# True value of every reported parameter, keyed as the parameter table is
# after alignment: <parameter>.<level>.<class or shared>.<indicator>.
truth_table <- function(cell) {
  structure <- .additive_structure(cell$family, cell$between_variance)
  block <- function(values, mode, name, level) {
    rows <- if (mode == "equal") 1L else seq_len(cell$H)
    labels <- if (mode == "equal") "shared" else as.character(rows)
    data.frame(parameter_id = sprintf("%s.%s.%s.%s", name, level,
                                      rep(labels, each = 2L), c("y1", "y2")),
               truth = as.vector(t(values[rows, , drop = FALSE])))
  }
  rbind(block(cell$means, structure$means, "mean", "between"),
        block(cell$within, structure$within, "variance", "within"),
        block(cell$between, structure$between, "variance", "between"),
        data.frame(parameter_id = sprintf("weight.group_class.%d.NA", seq_len(cell$H)),
                   truth = cell$weights))
}

run_one <- function(cell_name, replicate) {
  cell <- cells[[cell_name]]
  seed <- registry$base_seed[registry$cell == cell_name] + replicate
  data <- draw(cell, seed)
  conditions <- character()
  record <- function(condition) {
    conditions <<- c(conditions, class(condition)[1L])
    invokeRestart(if (inherits(condition, "warning")) "muffleWarning" else
      "muffleMessage")
  }
  fit <- withCallingHandlers(
    multilpa(data, c("y1", "y2"), "group", n_group_classes = cell$H,
             family = cell$family, between_variance = cell$between_variance,
             seed = seed), warning = record, message = record)
  tables <- lapply(c("observed", "robust"), function(type) {
    withCallingHandlers(get_results(fit, "parameters", vcov_type = type),
                        warning = record, message = record)
  })
  # Align classes by means and log within variances jointly: dispersion
  # classes differ only in spread.
  signature <- cbind(fit$means, log(fit$within_variances))
  target <- cbind(cell$means, log(cell$within))
  swap <- sum((signature[2:1, ] - target)^2) < sum((signature - target)^2)
  class_index <- if (swap) c(2L, 1L) else c(1L, 2L)
  key <- function(table) {
    label <- ifelse(table$group_class == "shared", "shared",
                    as.character(class_index[match(table$group_class,
                                                   c("group_class_1", "group_class_2"))]))
    sprintf("%s.%s.%s.%s", table$parameter, table$level, label,
            ifelse(is.na(table$indicator), "NA", table$indicator))
  }
  truth <- truth_table(cell)
  observed <- tables[[1L]]
  robust <- tables[[2L]]
  rows <- match(truth$parameter_id, key(observed))
  stopifnot("every true parameter must be matched once" = !anyNA(rows))
  data.frame(cell = cell_name, replicate = replicate, seed = seed, truth,
             estimate = observed$estimate[rows],
             se_observed = observed$standard_error[rows],
             se_robust = robust$standard_error[rows],
             low_observed = observed$conf_low[rows],
             high_observed = observed$conf_high[rows],
             low_robust = robust$conf_low[rows], high_robust = robust$conf_high[rows],
             converged = fit$converged,
             boundary = any(fit$between_zero) || any(fit$within_floor),
             conditions = paste(unique(conditions), collapse = ";"),
             stringsAsFactors = FALSE)
}

started <- Sys.time()
jobs <- expand.grid(cell = names(cells), replicate = seq_len(n_datasets),
                    stringsAsFactors = FALSE)
replicates <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  tryCatch(run_one(jobs$cell[k], jobs$replicate[k]), error = function(error) {
    data.frame(cell = jobs$cell[k], replicate = jobs$replicate[k], seed = NA_real_,
               parameter_id = NA_character_, truth = NA_real_, estimate = NA_real_,
               se_observed = NA_real_, se_robust = NA_real_, low_observed = NA_real_,
               high_observed = NA_real_, low_robust = NA_real_, high_robust = NA_real_,
               converged = NA, boundary = NA, conditions = paste("error:", conditionMessage(error)),
               stringsAsFactors = FALSE)
  })
}, mc.cores = max(1L, parallel::detectCores() - 1L), mc.set.seed = TRUE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
utils::write.csv(replicates, file.path(out_dir, "families-simulation-replicates.csv"),
                 row.names = FALSE)

usable <- subset(replicates, !is.na(parameter_id))
keys <- unique(usable[, c("cell", "parameter_id")])
summary_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  rows <- subset(usable, cell == keys$cell[k] & parameter_id == keys$parameter_id[k])
  ok_observed <- is.finite(rows$se_observed)
  ok_robust <- is.finite(rows$se_robust)
  cover <- function(low, high, ok) mean(rows$truth[ok] >= low[ok] & rows$truth[ok] <= high[ok])
  sd_estimate <- stats::sd(rows$estimate)
  data.frame(cell = keys$cell[k], parameter_id = keys$parameter_id[k],
             truth = rows$truth[1L], attempted = n_datasets,
             availability = sum(ok_observed) / n_datasets,
             bias_over_sd = (mean(rows$estimate) - rows$truth[1L]) / sd_estimate,
             se_ratio_observed = mean(rows$se_observed[ok_observed]) /
               stats::sd(rows$estimate[ok_observed]),
             coverage_observed = cover(rows$low_observed, rows$high_observed, ok_observed),
             se_ratio_robust = mean(rows$se_robust[ok_robust]) /
               stats::sd(rows$estimate[ok_robust]),
             coverage_robust = cover(rows$low_robust, rows$high_robust, ok_robust),
             stringsAsFactors = FALSE)
}))
summary_table$meets_release <- with(summary_table,
  availability >= 0.98 & abs(bias_over_sd) <= 0.1 &
    se_ratio_observed >= 0.9 & se_ratio_observed <= 1.1 &
    coverage_observed >= 0.925 & coverage_observed <= 0.975)
utils::write.csv(summary_table, file.path(out_dir, "families-simulation-summary.csv"),
                 row.names = FALSE)
status <- do.call(rbind, lapply(names(cells), function(name) {
  per <- unique(subset(replicates, cell == name)[, c("replicate", "converged",
                                                     "boundary", "conditions")])
  data.frame(cell = name, attempted = n_datasets,
             errors = sum(startsWith(per$conditions, "error:")),
             converged = sum(per$converged %in% TRUE),
             boundary = sum(per$boundary %in% TRUE),
             weak_class_warning = sum(grepl("latents_weak_class", per$conditions)))
}))
sink(file.path(out_dir, "families-simulation.txt"))
cat(sprintf("Families simulation: %d datasets per cell, %.1f minutes; %s\n\n",
            n_datasets, elapsed, R.version.string))
print(status, row.names = FALSE)
cat("\n")
print(summary_table, digits = 3, row.names = FALSE)
sink()
cat(sprintf("done in %.1f minutes\n", elapsed))
