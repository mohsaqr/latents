#' A small separated Gaussian panel for public diagnostic-table checks
#' @return Group IDs, occasions, a covariate and two continuous indicators.
skip_on_cran()

gaussian_table_data <- function() {
  set.seed(3)
  n <- 90L
  location <- rep(c(-2, 2), each = n / 2L)
  data.frame(id = rep(seq_len(n / 3L), each = 3L), time = rep(seq_len(3L), n / 3L),
             a = rnorm(n) + location, b = rnorm(n) + location, z = rnorm(n))
}

test_that("ordinary and covariate fits expose the same Gaussian table in summaries", {
  data <- gaussian_table_data()
  fits <- list(
    lpa(data, c("a", "b"), 2L, n_starts = 1L, seed = 1),
    lpa(data, c("a", "b"), 2L, profile_covariates = "z", n_starts = 1L, seed = 1))
  invisible(lapply(fits, function(fit) {
    table <- get_results(fit, "gaussian_diagnostics")
    expect_identical(table, fit$gaussian_diagnostics)
    expect_identical(get_results(summary(fit), "gaussian_diagnostics"), table)
    expect_equal(table$context, rep("profiles", 2L))
    expect_equal(table$profile, paste0("profile_", seq_len(2L)))
    expect_equal(table$gaussian_mass, unname(colSums(fit$subject_posteriors)))
    reference_variance <- vapply(data[c("a", "b")], function(x) mean((x - mean(x))^2),
                                 numeric(1))
    expected <- apply(sweep(fit$variances, 2L, reference_variance, "/"), 1L, min)
    expect_equal(table$min_standardized_eigenvalue, unname(expected))
    expect_output(print(fit), "Smallest standardized Gaussian eigenvalue")
    legacy <- fit
    legacy$gaussian_diagnostics <- NULL
    expect_false(any(grepl("Smallest standardized Gaussian eigenvalue",
                           capture.output(print(legacy)), fixed = TRUE)))
  }))
})

test_that("homogeneous and general LTA Gaussian diagnostics follow measurement occasions", {
  data <- gaussian_table_data()
  shared <- lta(data, c("a", "b"), "id", 2L, time = "time", n_starts = 1L, seed = 1)
  occasion <- quietly(lta(
    data, c("a", "b"), "id", 2L, time = "time", transitions = "occasion",
    measurement = "occasion", n_starts = 1L, seed = 1, max_iter = 0L),
    "latents_unconverged")
  shared_table <- get_results(shared, "gaussian_diagnostics")
  expect_identical(shared_table, shared$gaussian_diagnostics)
  expect_identical(get_results(suppressMessages(summary(shared)), "gaussian_diagnostics"),
                   shared_table)
  expect_equal(shared_table$context, rep("profiles", 2L))
  expect_equal(shared_table$gaussian_mass, unname(colSums(shared$subject_posteriors)))
  occasion_table <- get_results(occasion, "gaussian_diagnostics")
  expect_identical(occasion_table, occasion$gaussian_diagnostics)
  # A general LTA summary serves its tables through get_results(), as every
  # other summary does (release review API-6).
  expect_identical(get_results(suppressMessages(summary(occasion)), "gaussian_diagnostics"),
                   occasion_table)
  expect_equal(occasion_table$context, rep(paste("occasion", 1:3), each = 2L))
  expect_equal(occasion_table$profile, rep(paste0("profile_", 1:2), 3L))
  invisible(lapply(seq_len(3L), function(time) {
    rows <- data$time == time
    index <- occasion_table$context == paste("occasion", time)
    expect_equal(occasion_table$gaussian_mass[index],
                 unname(colSums(occasion$subject_posteriors[rows, ])))
    variances <- vapply(data[rows, c("a", "b")], function(x) mean((x - mean(x))^2),
                        numeric(1))
    expected <- apply(sweep(occasion$measurement[[time]]$variances, 2L, variances, "/"),
                      1L, min)
    expect_equal(occasion_table$min_standardized_eigenvalue[index], unname(expected))
  }))
  expect_output(print(shared), "Smallest standardized Gaussian eigenvalue")
  expect_output(print(occasion), "Smallest standardized Gaussian eigenvalue")
})

test_that("weighted Gaussian tables use observed support and weighted reference variances", {
  data <- gaussian_table_data()
  data$weight <- rep(c(.25, 1, 2), 30L)
  data$a[1L] <- NA_real_
  data$b[2L] <- NA_real_
  data[3L, c("a", "b")] <- NA_real_
  weights <- data$weight * nrow(data) / sum(data$weight)
  observed <- rowSums(!is.na(data[c("a", "b")])) > 0L
  reference_variance <- vapply(data[c("a", "b")], function(values) {
    available <- !is.na(values)
    center <- sum(values[available] * weights[available]) / sum(weights[available])
    sum(weights[available] * (values[available] - center)^2) / sum(weights[available])
  }, numeric(1))
  invisible(lapply(c(FALSE, TRUE), function(covariate) {
    fit <- lpa(data, c("a", "b"), 2L, weights = "weight", missing = "fiml",
               profile_covariates = if (covariate) "z" else character(),
               n_starts = 1L, seed = 1L)
    table <- get_results(fit, "gaussian_diagnostics")
    contributions <- fit$subject_posteriors * weights * observed
    expect_equal(table$gaussian_mass, unname(colSums(contributions)))
    expect_equal(table$weight_effective_n,
                 unname(colSums(contributions)^2 / colSums(contributions^2)))
    expected <- apply(sweep(fit$variances, 2L, reference_variance, "/"), 1L, min)
    expect_equal(table$min_standardized_eigenvalue, unname(expected))
    expect_equal(table$positive_rows, unname(colSums(contributions > 0)))
  }))
})

test_that("cross-level tables distinguish individual and group Gaussian densities", {
  data <- gaussian_table_data()
  invisible(lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    fit <- quietly(multilpa(
      data, c("a", "b"), "id", 2L, 2L, family = family,
      between_variance = "varying", n_starts = 1L, seed = 1L, max_iter = 0L),
      "latents_unconverged")
    table <- get_results(fit, "gaussian_diagnostics")
    expect_identical(table, fit$gaussian_diagnostics)
    expect_identical(summary(fit)$gaussian_diagnostics, table)
    expect_equal(table$context, rep(c("individual profiles", "group classes"), each = 2L))
    individual <- table$context == "individual profiles"
    group <- table$context == "group classes"
    expect_equal(table$gaussian_mass[individual], unname(colSums(fit$subject_posteriors)))
    expect_equal(table$gaussian_mass[group], unname(colSums(fit$group_posteriors)))
    expect_equal(table$profile[individual], rownames(fit$profile_variances))
    expect_equal(table$profile[group], rownames(fit$group_mean_variances))
    expect_output(print(fit), "Smallest standardized Gaussian eigenvalue")
  }))
})

test_that("moderate flat components remain visible without a new selection threshold", {
  expect_no_warning(fit <- lpa(srl, names(srl), 4L, model = "VVV", seed = 1))
  table <- get_results(fit, "gaussian_diagnostics")
  expect_true(fit$converged)
  expect_false(fit$boundary)
  expect_gt(min(table$gaussian_mass), 20)
  expect_lt(min(table$min_standardized_eigenvalue), .002)
  expect_gt(min(table$min_standardized_eigenvalue), .001)
  expect_output(print(fit), "Smallest standardized Gaussian eigenvalue: 0.001")
  expect_equal(fit$log_likelihood, max(fit$starts$log_likelihood, na.rm = TRUE),
               tolerance = 1e-6)
})

test_that("categorical models and old objects have no invented Gaussian diagnostics", {
  data <- gaussian_table_data()
  data$c <- factor(ifelse(data$a > 0, "high", "low"))
  data$d <- factor(ifelse(data$b > 0, "high", "low"))
  # Two classes of two binary items are not identified (CORE-5).
  expect_warning(profile_fit <- quietly(lpa(data, c("c", "d"), 2L, categorical = c("c", "d"),
                                            n_starts = 1L, seed = 1L, max_iter = 0L),
                                        "latents_unconverged"),
                 class = "latents_unidentified")
  fits <- list(
    profile_fit,
    quietly(lta(data, c("c", "d"), "id", 2L, categorical = c("c", "d"),
                 time = "time", transitions = "occasion", n_starts = 1L,
                 seed = 1L, max_iter = 0L), "latents_unconverged"))
  invisible(lapply(fits, function(fit) {
    expect_equal(nrow(fit$gaussian_diagnostics), 0L)
    expect_false("gaussian_diagnostics" %in% names(suppressMessages(get_results(fit, "all"))))
    expect_false(any(grepl("Smallest standardized Gaussian eigenvalue",
                           capture.output(print(fit)), fixed = TRUE)))
    expect_true(is.na(.latents_gaussian_minimum(fit, "min_standardized_eigenvalue")))
  }))
  legacy <- fits[[1L]]
  legacy$gaussian_diagnostics <- NULL
  expect_true(is.na(.latents_gaussian_minimum(legacy, "gaussian_mass")))
})

test_that("enumeration exposes continuous minima from the returned Gaussian fits", {
  data <- gaussian_table_data()
  grid <- enumerate_lpa(data, c("a", "b"), n_profiles = 1:2, model = "VVI",
                        n_starts = 1L, seed = 1L)
  expect_true(all(c("min_standardized_eigenvalue", "min_gaussian_mass") %in%
                    names(grid$table)))
  expected_eigenvalue <- vapply(grid$fits, function(fit) {
    min(fit$gaussian_diagnostics$min_standardized_eigenvalue)
  }, numeric(1))
  expected_mass <- vapply(grid$fits, function(fit) min(fit$gaussian_diagnostics$gaussian_mass),
                          numeric(1))
  expect_equal(grid$table$min_standardized_eigenvalue, unname(expected_eigenvalue))
  expect_equal(grid$table$min_gaussian_mass, unname(expected_mass))
})
