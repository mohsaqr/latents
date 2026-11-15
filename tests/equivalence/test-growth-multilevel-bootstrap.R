# The cluster bootstrap of the multilevel growth mixture
# (parameter_inference(method = "bootstrap") on mixture_regression(random =,
# cluster =, n_group_classes =)). Two references that share no code with the
# matching it checks:
#
# - label matching against brute force: on real cluster resamples, the
#   matching by posterior memberships (persons for the classes, clusters for
#   the group classes) must choose the same relabelling as an exhaustive
#   search over every pair of class and group-class permutations for the one
#   nearest the original's natural parameters;
# - known parameters: the retained, predeclared simulation
#   (validation/growth-simulation-multilevel-bootstrap.R) is re-summarised
#   from its replicates, and its criteria, fixed in the script before the
#   first replicate, must hold: bias within 3.5 Monte Carlo standard errors,
#   coverage of the bootstrap and Wald 95% intervals within the central 99.9%
#   binomial range, bootstrap / Hessian standard-error ratio in [0.8, 1.25],
#   and the matching keeping the class and group-class order in >= 99% of
#   resamples.

skip_on_cran()

ml_boot_data <- function(seed, n_schools = 30L, per_school = 8L) {
  withr::with_seed(seed, {
    group_class <- 1L + (stats::runif(n_schools) >= 0.5)
    do.call(rbind, lapply(seq_len(n_schools), function(s) {
      do.call(rbind, lapply(seq_len(per_school), function(i) {
        k <- if (stats::runif(1) < c(0.85, 0.2)[group_class[s]]) 1L else 2L
        effects <- c(stats::rnorm(1, 0, 2), stats::rnorm(1, 0, 0.4))
        data.frame(school = s, student = (s - 1L) * per_school + i, wave = 0:3,
                   score = c(50, 55)[k] + c(3, 0)[k] * 0:3 + effects[1L] +
                     effects[2L] * 0:3 + stats::rnorm(4, 0, 1.5))
      }))
    }))
  })
}

ml_boot_fit <- function(data, n_starts = 3L) {
  mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                     class_level = "group", random = ~ 1 + wave,
                     random_covariance = "equal", n_group_classes = 2,
                     cluster = "school", n_starts = n_starts, seed = 1,
                     vcov_type = "none")
}

test_that("the matching chooses the brute-force nearest relabelling on cluster resamples", {
  skip_if_not_installed("withr")
  data <- ml_boot_data(11L)
  fit <- ml_boot_fit(data)
  natural <- latents:::.growth_natural_parameters(fit)
  original <- natural$values(fit$params, fit$spec)
  # Rows whose value moves under a relabelling (logits change sign, shares
  # and coefficients swap); the shared random effects do not.
  rows_by_cluster <- split(seq_len(nrow(fit$spec$model_data)),
                           fit$spec$cluster_index[fit$spec$group_index])
  permutations <- list(1:2, 2:1)
  flips <- 0L
  for (r in 1:5) {
    drawn <- withr::with_seed(100L + r, sample.int(fit$spec$n_clusters, replace = TRUE))
    resampled <- latents:::.multilpa_resample_groups(fit$spec$model_data, rows_by_cluster,
                                                     "school", drawn)
    resampled <- latents:::.growth_renumber_persons(resampled, "school", "student")
    refit <- suppressWarnings(ml_boot_fit(resampled, n_starts = 2L))
    expect_true(refit$converged)
    # Relabel the refit at random, so that both levels are exercised.
    scramble <- withr::with_seed(200L + r, list(sample(2L), sample(2L)))
    refit$params <- latents:::.growth_permute_params(refit$spec, refit$params,
                                                     scramble[[1L]], scramble[[2L]])
    refit$expectation <- latents:::.growth_expectation(refit$spec, refit$stats, refit$params)
    flips <- flips + !identical(scramble, list(1:2, 1:2))
    # Brute force: every class order with every group-class order, the
    # nearest to the original's natural parameters (on their own scale,
    # standardised by the original's Wald-free magnitude).
    candidates <- expand.grid(class = 1:2, group = 1:2)
    distance <- vapply(seq_len(nrow(candidates)), function(j) {
      params <- latents:::.growth_permute_params(
        refit$spec, refit$params, permutations[[candidates$class[j]]],
        permutations[[candidates$group[j]]])
      values <- natural$values(params, refit$spec)
      sum(((values - original) / pmax(abs(original), 1))^2)
    }, numeric(1))
    best <- candidates[which.min(distance), ]
    matched <- latents:::.growth_match_classes(fit, refit)
    expect_identical(matched$order, permutations[[best$class]])
    expect_identical(matched$group_order, permutations[[best$group]])
    # The nearest relabelling is clearly nearest: no tie decides it.
    expect_lt(min(distance), 0.5 * sort(distance)[2L])
  }
  expect_gte(flips, 1L)
})

test_that("the predeclared multilevel growth bootstrap simulation meets its criteria", {
  directory <- test_path("..", "..", "validation")
  files <- file.path(directory, paste0("growth-simulation-multilevel-bootstrap-",
                                       c("registry", "replicates", "summary"), ".csv"))
  skip_if_not(all(file.exists(files)), "the retained simulation results are not available")
  registry <- utils::read.csv(files[1L], stringsAsFactors = FALSE)
  replicates <- utils::read.csv(files[2L], stringsAsFactors = FALSE)
  summary <- utils::read.csv(files[3L], stringsAsFactors = FALSE)
  # The script declares the same registry it retained.
  script <- readLines(file.path(directory, "growth-simulation-multilevel-bootstrap.R"))
  expect_true(any(grepl(sprintf("base_seed = %dL", registry$base_seed), script, fixed = TRUE)))
  n_datasets <- length(unique(replicates$replicate))
  expect_identical(sort(unique(replicates$replicate)), seq_len(n_datasets))
  ok <- replicates[replicates$available, ]
  # Availability, counted over every dataset drawn.
  expect_gte(length(unique(ok$replicate)) / n_datasets, 0.98)
  # Re-summarise from the replicates: the stored summary is theirs.
  keys <- paste(summary$block, summary$class, summary$term)
  for (r in seq_len(nrow(summary))) {
    rows <- ok[paste(ok$block, ok$class, ok$term) == keys[r], ]
    n <- nrow(rows)
    expect_identical(n, summary$n[r])
    expect_equal(mean(rows$estimate) - rows$truth[1L], summary$bias[r], tolerance = 1e-8)
    range <- stats::qbinom(c(0.0005, 0.9995), n, 0.95) / n
    boot <- mean(rows$low_boot <= rows$truth & rows$truth <= rows$high_boot)
    wald <- mean(rows$low_wald <= rows$truth & rows$truth <= rows$high_wald)
    expect_equal(boot, summary$coverage_boot[r], tolerance = 1e-12)
    # The criteria fixed before the first replicate.
    label <- keys[r]
    expect_lte(abs(mean(rows$estimate) - rows$truth[1L]),
               3.5 * stats::sd(rows$estimate) / sqrt(n), label = paste(label, "bias"))
    expect_true(boot >= range[1L] && boot <= range[2L],
                label = sprintf("%s bootstrap coverage %.3f", label, boot))
    expect_true(wald >= range[1L] && wald <= range[2L],
                label = sprintf("%s Wald coverage %.3f", label, wald))
    ratio <- stats::median(rows$se_boot / rows$se_wald)
    expect_true(ratio >= 0.8 && ratio <= 1.25,
                label = sprintf("%s bootstrap/Wald SE ratio %.3f", label, ratio))
  }
  # Every natural parameter of the design is summarised.
  expect_identical(nrow(summary), 18L)
  # Label matching kept the class and group-class order across resamples.
  first <- ok[!duplicated(ok$replicate), ]
  expect_gte(mean(first$class_kept), 0.99)
  expect_gte(mean(first$group_kept), 0.99)
})
