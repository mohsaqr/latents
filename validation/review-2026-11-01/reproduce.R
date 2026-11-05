# Review probes against commit 4c8f5aa7. Run from the repository root.
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
quiet <- function(expr) withCallingHandlers(expr,
  warning = function(w) invokeRestart("muffleWarning"),
  message = function(m) invokeRestart("muffleMessage"))
probe <- function(name, expr) {
  cat("\nPROBE:", name, "\n")
  tryCatch(force(expr), error = function(e) {
    cat("ERROR CLASS:", class(e)[1L], "\nMESSAGE:", conditionMessage(e), "\n")
  })
}
probe("LRT shared versus group_class slopes", {
  set.seed(917)
  g <- rep(1:60, each = 8)
  z <- rnorm(length(g))
  group <- rep(rep(1:2, each = 30), each = 8)
  state <- 1L + rbinom(length(g), 1, plogis(c(-1.8, 1.8)[group] + 0.5 * z))
  d <- data.frame(g, z, a = rnorm(length(g), c(-2, 2)[state]),
                  b = rnorm(length(g), c(-1.5, 1.5)[state]))
  fit <- function(slopes) quiet(multilpa(d, c("a", "b"), "g", 2, 2,
    profile_covariates = "z", profile_slopes = slopes, model = "EEI",
    n_starts = 3, seed = 1, tol = 1e-9, max_iter = 5000))
  shared <- fit("shared")
  separate <- fit("group_class")
  cat("Converged:", shared$converged, separate$converged, "\n")
  print(as.data.frame(compare_models(shared, separate, test = "lrt")))
})
probe("BCH out of range pairwise inference", {
  set.seed(10)
  truth <- sample.int(2L, 300L, replace = TRUE, prob = c(0.8, 0.2))
  d <- data.frame(y = rnorm(300L, c(0, 1.2)[truth]), z = rnorm(300L, c(0, 1.2)[truth]))
  d$outcome <- factor(ifelse(runif(300L) < c(0.02, 0.3)[truth], "rare", "common"))
  fit <- quiet(lpa(d, c("y", "z"), 2, n_starts = 2, seed = 1))
  print(as.data.frame(quiet(three_step(fit, d, "outcome"))))
  print(as.data.frame(quiet(three_step(fit, d, "outcome", contrast = "pairs"))))
})
probe("align ordinal category labels", {
  set.seed(12)
  n <- 400
  state <- sample(1:2, n, TRUE)
  d <- data.frame(a = rnorm(n, c(-2, 2)[state]),
                  u = ordered(sample(c("low", "medium", "high"), n, TRUE),
                              levels = c("low", "medium", "high")))
  e <- d
  levels(e$u) <- c("red", "green", "blue")
  first <- quiet(lpa(d, c("a", "u"), 2, ordinal = "u", n_starts = 1, seed = 1))
  second <- quiet(lpa(e, c("a", "u"), 2, ordinal = "u", n_starts = 1, seed = 1))
  cat("Reference levels:", first$extra_data$ordinal_levels$u, "\n")
  cat("Fit levels:", second$extra_data$ordinal_levels$u, "\n")
  cat("Accepted class:", class(align_labels(second, first)), "\n")
})
probe("covariate alignment posterior column labels", {
  set.seed(3)
  n <- 450
  z <- rnorm(n)
  k <- vapply(z, function(v) sample.int(3, 1, prob = exp(c(0.8, 0, -0.8) * v)), integer(1))
  means <- rbind(c(-3, 0), c(0, 3), c(3, -3))
  d <- data.frame(a = means[k, 1] + rnorm(n), b = means[k, 2] + rnorm(n), z)
  fit <- quiet(lpa(d, c("a", "b"), 3, model = "VVV", profile_covariates = "z",
                   n_starts = 2, seed = 1))
  swapped <- .multilpa_cov_permute(fit, c(3L, 1L, 2L), 1L)
  aligned <- align_labels(fit, swapped)
  cat("Measurement row names:", rownames(aligned$means), "\n")
  cat("Posterior column names:", colnames(aligned$subject_posteriors), "\n")
  cat("Covariance dimension names:", dimnames(aligned$covariances)[[3]], "\n")
  cat("Posterior names absent before alignment too:", is.null(colnames(fit$subject_posteriors)), "\n")
})
probe("confint for fixed covariate measurement", {
  set.seed(3)
  n <- 600
  z <- rnorm(n)
  state <- vapply(z, function(v) sample.int(3, 1, prob = exp(c(0.8, 0, -0.8) * v)), integer(1))
  centres <- rbind(c(-2, 0), c(0, 2), c(2, -2))
  d <- data.frame(a = centres[state, 1] + rnorm(n), b = centres[state, 2] + rnorm(n), z)
  plain <- quiet(lpa(d, c("a", "b"), 3, model = "VVI", n_starts = 3, seed = 1))
  held <- quiet(lpa(d, c("a", "b"), 3, model = "VVI", profile_covariates = "z",
                    start = starting_values(plain), fixed = "measurement", n_starts = 3, seed = 1))
  cat("Free parameters:", held$n_parameters, "\n")
  cat("Inference rows:", nrow(parameter_inference(held)), "\n")
  cat("vcov dimension:", dim(vcov(held)), "\n")
  print(confint(held))
})
probe("LRT shared predictor with changed values", {
  set.seed(300)
  n <- 600
  z <- rnorm(n)
  w <- rnorm(n)
  state <- 1L + rbinom(n, 1, plogis(2 * z))
  d <- data.frame(a = rnorm(n, c(-1.5, 1.5)[state]),
                  b = rnorm(n, c(-1, 1)[state]), z, w)
  e <- d
  e$z <- sample(d$z)
  small <- quiet(lpa(e, c("a", "b"), 2, model = "EEI", profile_covariates = "z",
                     n_starts = 3, seed = 1))
  large <- quiet(lpa(d, c("a", "b"), 2, model = "EEI", profile_covariates = c("z", "w"),
                     n_starts = 3, seed = 1))
  cat("Same predictor name, equal values:", identical(d$z, e$z), "\n")
  cat("Converged:", small$converged, large$converged, "\n")
  print(as.data.frame(compare_models(small, large, test = "lrt"))[
    c("model", "lr_statistic", "lr_df", "lr_p_value")])
})
probe("bootstrap plot in small measurement units", {
  d <- transform(datasets::faithful, eruptions = eruptions * 1e-12, waiting = waiting * 1e-12)
  fit <- quiet(lpa(d, c("eruptions", "waiting"), 2, model = "VVI", n_starts = 2, seed = 1))
  boot <- quiet(parameter_inference(fit, method = "bootstrap", iter = 6, n_starts = 1, seed = 1))
  means <- which(boot$parameter == "mean")
  cat("Mean replicate ranges:", apply(attr(boot, "replicates")[, means], 2, function(x) diff(range(x))), "\n")
  drawn <- plot(boot, what = "mean")
  cat("Drawn mean parameters:", nrow(drawn$data), "\n")
})
probe("bootstrap row reorder followed by row name normalization", {
  fit <- quiet(lpa(datasets::faithful, c("eruptions", "waiting"), 2,
                   model = "VVI", n_starts = 2, seed = 1))
  boot <- quiet(parameter_inference(fit, method = "bootstrap", iter = 6, n_starts = 1, seed = 1))
  changed <- boot[rev(seq_len(nrow(boot))), , drop = FALSE]
  rownames(changed) <- NULL
  tidy <- as.data.frame(changed, what = "replicates")
  key <- function(x) with(x, paste(level, outcome, term, parameter))
  observed_sd <- vapply(key(changed), function(k) sd(tidy$estimate[key(tidy) == k]), numeric(1))
  cat("Maximum difference from the matching reported standard error:",
      max(abs(observed_sd - changed$standard_error)), "\n")
  if (requireNamespace("dplyr", quietly = TRUE)) {
    changed <- dplyr::arrange(boot, dplyr::desc(estimate))
    cat("dplyr::arrange preserves row names:", identical(rownames(changed), rownames(boot)), "\n")
    tidy <- as.data.frame(changed, what = "replicates")
    observed_sd <- vapply(key(changed), function(k) sd(tidy$estimate[key(tidy) == k]), numeric(1))
    cat("dplyr::arrange maximum standard error mismatch:",
        max(abs(observed_sd - changed$standard_error)), "\n")
  }
})
probe("categorical outcome category present only in excluded rows", {
  d <- datasets::faithful
  d$outcome <- factor(rep(c("yes", "no"), length.out = nrow(d)))
  fit <- quiet(lpa(d, c("eruptions", "waiting"), 2, n_starts = 2, seed = 1))
  base <- quiet(three_step(fit, d, "outcome", method = "ml"))
  e <- rbind(transform(d, outcome = as.character(outcome)),
              data.frame(eruptions = NA, waiting = NA, outcome = "excluded"))
  e$outcome <- factor(e$outcome)
  empty <- quiet(lpa(e, c("eruptions", "waiting"), 2, n_starts = 2, seed = 1,
                     missing = "fiml"))
  added <- quiet(three_step(empty, e, "outcome", method = "ml"))
  cat("Finite errors before:", sum(is.finite(base$standard_error)), "\n")
  cat("Finite errors after:", sum(is.finite(added$standard_error)), "\n")
  print(as.data.frame(added))
})
probe("LRT ordinal category order mismatch", {
  set.seed(931)
  n <- 600
  z <- rnorm(n)
  state <- 1L + rbinom(n, 1, plogis(3 * z))
  labels <- c("low", "middle", "high")
  probabilities <- t(vapply(c(-1, 1), function(location) {
    value <- exp(c(0, 0.1, -0.5) + 0:2 * location)
    value / sum(value)
  }, numeric(3)))
  u <- vapply(state, function(k) sample(labels, 1, prob = probabilities[k, ]), character(1))
  d <- data.frame(a = rnorm(n, c(-1.5, 1.5)[state]),
                  b = rnorm(n, c(-1, 1)[state]), u = ordered(u, levels = labels), z)
  e <- d
  e$u <- ordered(e$u, levels = c("low", "high", "middle"))
  first <- quiet(lpa(d, c("a", "b", "u"), 2, model = "EEI", ordinal = "u",
                     n_starts = 3, seed = 1))
  second <- quiet(lpa(e, c("a", "b", "u"), 2, model = "EEI", ordinal = "u",
                      profile_covariates = "z", n_starts = 3, seed = 1))
  cat("Null category order:", first$extra_data$ordinal_levels$u, "\n")
  cat("Alternative category order:", second$extra_data$ordinal_levels$u, "\n")
  cat("Nesting checker returns NULL:", is.null(.multilpa_lrt_nesting(first, second)), "\n")
  print(as.data.frame(compare_models(first, second, test = "lrt"))[
    c("model", "lr_statistic", "lr_df", "lr_p_value")])
})
