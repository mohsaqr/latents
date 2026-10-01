# Shared setup for the additive-family simulation scripts: the cell registry
# (verified unchanged against additive-simulation-registry.csv), data
# generation, class alignment and the per-replicate fit. Sourced by
# additive-simulation.R and additive-diagnostic.R from the project root.

devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
out_dir <- "validation"

# ---- Registry -----------------------------------------------------------------
two <- function(a, b) rbind(a, b)
regular <- list(H = 2L, J = 200L, sizes = 10L, weights = c(0.5, 0.5),
                means = two(c(-1, -1), c(1, 1)), within = c(1, 1),
                between = two(c(0.25, 0.5), c(0.5, 0.25)), fit = "varying",
                residual = "normal")
cells <- list(
  regular_varying = regular,
  regular_shared = modifyList(regular, list(between = two(c(0.4, 0.4), c(0.4, 0.4)),
                                            fit = "equal")),
  small_groups = modifyList(regular, list(J = 50L)),
  unequal_sizes = modifyList(regular, list(sizes = c(2L, 5L, 10L, 20L))),
  singleton_admixture = modifyList(regular, list(sizes = c(1L, 10L, 10L, 10L))),
  three_classes = list(H = 3L, J = 300L, sizes = 10L, weights = rep(1 / 3, 3),
                       means = rbind(c(-1.5, -1.5), c(0, 0), c(1.5, 1.5)),
                       within = c(1, 1), between = matrix(0.3, 3, 2),
                       fit = "equal", residual = "normal"),
  rare_weak = modifyList(regular, list(weights = c(0.9, 0.1),
                                       means = two(c(-0.4, -0.4), c(0.4, 0.4)),
                                       between = two(c(0.4, 0.4), c(0.4, 0.4)),
                                       fit = "equal")),
  boundary_zero = modifyList(regular, list(between = matrix(0, 2, 2), fit = "equal")),
  boundary_small = modifyList(regular, list(between = matrix(0.01, 2, 2),
                                            fit = "equal")),
  misspecified_correlated = modifyList(regular, list(residual = "correlated")),
  misspecified_t5 = modifyList(regular, list(residual = "t5"))
)
registry <- do.call(rbind, lapply(seq_along(cells), function(k) {
  cell <- cells[[k]]
  data.frame(cell = names(cells)[k], n_classes = cell$H, n_groups = cell$J,
             sizes = paste(cell$sizes, collapse = "/"),
             weights = paste(signif(cell$weights, 4), collapse = "/"),
             means = paste(apply(cell$means, 1L, paste, collapse = ","), collapse = " | "),
             between = paste(apply(cell$between, 1L, paste, collapse = ","),
                             collapse = " | "),
             within = paste(cell$within, collapse = ","),
             residual = cell$residual, fitted_between_variance = cell$fit,
             base_seed = 1e6L * k, stringsAsFactors = FALSE)
}))
registry_file <- file.path(out_dir, "additive-simulation-registry.csv")
if (file.exists(registry_file)) {
  stopifnot("registry changed after it was declared" =
              identical(utils::read.csv(registry_file, stringsAsFactors = FALSE),
                        utils::read.csv(text = paste(utils::capture.output(
                          utils::write.csv(registry, row.names = FALSE)),
                          collapse = "\n"), stringsAsFactors = FALSE)))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}

# ---- Data generation ------------------------------------------------------------
draw <- function(cell, seed) {
  set.seed(seed)
  sizes <- rep_len(cell$sizes, cell$J)
  classes <- sample.int(cell$H, cell$J, replace = TRUE, prob = cell$weights)
  d <- ncol(cell$means)
  blocks <- lapply(seq_len(cell$J), function(j) {
    intercept <- stats::rnorm(d, cell$means[classes[j], ],
                              sqrt(cell$between[classes[j], ]))
    noise <- switch(cell$residual,
      normal = matrix(stats::rnorm(sizes[j] * d), sizes[j], d),
      # Correlation 0.4 between indicators, unit variances.
      correlated = matrix(stats::rnorm(sizes[j] * d), sizes[j], d) %*%
        chol(matrix(c(1, 0.4, 0.4, 1), 2)),
      # t with 5 df scaled to unit variance.
      t5 = matrix(stats::rt(sizes[j] * d, 5) * sqrt(3 / 5), sizes[j], d))
    ratings <- sweep(noise %*% diag(sqrt(cell$within), d), 2L, intercept, "+")
    data.frame(group = j, y1 = ratings[, 1L], y2 = ratings[, 2L])
  })
  do.call(rbind, blocks)
}

truth_table <- function(cell) {
  shared <- identical(cell$fit, "equal")
  between <- if (shared) cell$between[1L, ] else as.vector(t(cell$between))
  data.frame(
    parameter_id = c(sprintf("mean.%d.%s", rep(seq_len(cell$H), each = 2), c("y1", "y2")),
                     sprintf("within.%s", c("y1", "y2")),
                     if (shared) sprintf("between.shared.%s", c("y1", "y2")) else
                       sprintf("between.%d.%s", rep(seq_len(cell$H), each = 2),
                               c("y1", "y2")),
                     sprintf("weight.%d", seq_len(cell$H))),
    truth = c(as.vector(t(cell$means)), cell$within, between, cell$weights),
    stringsAsFactors = FALSE)
}

# All permutations of 1..H (H <= 3 here).
permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  smaller <- permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))
  }))
}

# Align fitted classes to the generating ones by the permutation that
# minimizes the squared distance between class means.
align <- function(table, cell) {
  means <- subset(table, parameter == "mean")
  fitted <- matrix(means$estimate, cell$H, 2L, byrow = TRUE)
  options <- permutations(cell$H)
  loss <- vapply(seq_len(nrow(options)), function(k) {
    sum((fitted[options[k, ], , drop = FALSE] - cell$means)^2)
  }, numeric(1))
  options[which.min(loss), ]  # fitted class options[h] plays true class h
}

failure_row <- function(base, status, reason) {
  cbind(base, parameter_id = NA_character_, truth = NA_real_,
        estimate = NA_real_, se_observed = NA_real_, se_robust = NA_real_,
        low_observed = NA_real_, high_observed = NA_real_,
        low_robust = NA_real_, high_robust = NA_real_,
        fitted = FALSE, converged = NA, boundary = NA,
        inference = status, conditions = reason)
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
  fit <- tryCatch(withCallingHandlers(
    multilpa(data, c("y1", "y2"), "group", n_group_classes = cell$H,
             family = "additive", between_variance = cell$fit, seed = seed),
    warning = record, message = record),
    error = function(error) error)
  base <- data.frame(cell = cell_name, replicate = replicate, seed = seed,
                     stringsAsFactors = FALSE)
  if (inherits(fit, "error")) return(failure_row(base, "fit_failed", class(fit)[1L]))
  tables <- lapply(c("observed", "robust"), function(type) {
    tryCatch(withCallingHandlers(get_results(fit, "parameters", vcov_type = type),
                                 message = record, warning = record),
             error = function(error) {
               conditions <<- c(conditions, class(error)[1L])
               NULL
             })
  })
  observed <- tables[[1L]]
  robust <- tables[[2L]]
  mapping <- align(observed, cell)
  class_index <- match(observed$group_class, sprintf("group_class_%d", mapping))
  observed$parameter_id <- ifelse(
    observed$parameter == "mean", sprintf("mean.%d.%s", class_index, observed$indicator),
    ifelse(observed$parameter == "weight", sprintf("weight.%d", class_index),
      ifelse(observed$level == "within", sprintf("within.%s", observed$indicator),
        ifelse(observed$group_class == "shared",
               sprintf("between.shared.%s", observed$indicator),
               sprintf("between.%d.%s", class_index, observed$indicator)))))
  truth <- truth_table(cell)
  rows <- match(truth$parameter_id, observed$parameter_id)
  stopifnot("every true parameter must be matched once" = !anyNA(rows))
  boundary <- any(fit$between_zero) || any(fit$within_floor)
  has_observed <- any(is.finite(observed$standard_error))
  cbind(base, truth,
        estimate = observed$estimate[rows],
        se_observed = observed$standard_error[rows],
        se_robust = robust$standard_error[rows],
        low_observed = observed$conf_low[rows], high_observed = observed$conf_high[rows],
        low_robust = robust$conf_low[rows], high_robust = robust$conf_high[rows],
        fitted = TRUE, converged = fit$converged, boundary = boundary,
        inference = if (has_observed) "available" else "withheld",
        conditions = paste(unique(conditions), collapse = ";"))
}

