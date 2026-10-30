# Review 2026-10-26, inference partition: Wald information, robust and OPG
# covariances, Rubin's rules, label alignment and seed scoping, checked against
# independent computations rather than against the package's own arithmetic.
# The comparisons with independent references are in tests/equivalence/test-review-1026-inference.R.

skip_on_cran()

test_that("score-difference information recovers an analytic Hessian under both schemes", {
  set.seed(3)
  x <- stats::rnorm(40, 1, 2)
  # Normal log likelihood in (mu, eta = log sigma); score and -Hessian by hand.
  score <- function(theta) {
    sigma2 <- exp(2 * theta[2])
    c(sum(x - theta[1]) / sigma2, -length(x) + sum((x - theta[1])^2) / sigma2)
  }
  theta <- c(mu = 0.7, eta = 0.4)
  sigma2 <- exp(2 * theta[[2]])
  analytic <- matrix(c(length(x) / sigma2, 2 * sum(x - theta[[1]]) / sigma2,
                       2 * sum(x - theta[[1]]) / sigma2,
                       2 * sum((x - theta[[1]])^2) / sigma2), 2L)
  invisible(lapply(c("central", "five_point"), function(scheme) {
    invisible(lapply(c("plus1", "max1", "absolute"), function(rule) {
      information <- .inference_score_information(score, theta, 1e-4, scheme, rule)
      expect_equal(unname(information), analytic, tolerance = 1e-7)
      expect_identical(dimnames(information), list(names(theta), names(theta)))
    }))
  }))
  expect_error(.inference_invert(diag(c(1, -1))), class = "latents_singular_information")
  expect_error(.inference_invert(diag(c(1, 1e-12))), class = "latents_singular_information")
  expect_equal(.inference_invert(matrix(c(2, 1, 1, 2), 2L)), solve(matrix(c(2, 1, 1, 2), 2L)))
})

test_that("alignment recovers a known permutation and relabelling leaves the SE table invariant", {
  skip_on_cran()
  set.seed(8)
  n <- 240L
  state <- sample.int(3L, n, TRUE, prob = c(0.3, 0.3, 0.4))
  data <- data.frame(g = rep(1:40, each = 6),
                     a = c(-2, 0, 2)[state] + stats::rnorm(n, sd = c(1, 0.5, 1.5)[state]),
                     b = c(1, -1, 0)[state] + stats::rnorm(n))
  permutation <- c(3L, 1L, 2L)
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 3, n_group_classes = 1,
                            variance_model = "varying", covariance_model = covariance_model,
                            n_starts = 3, seed = 1))
    shuffled <- .multilpa_permute_profiles(fit, permutation)
    aligned <- .multilpa_pool_align(shuffled, fit)
    expect_identical(aligned$profile_order, order(permutation))
    expect_identical(coef(aligned$fit), coef(fit))
    expect_identical(aligned$fit$subject_profiles, fit$subject_profiles)
    original <- parameter_inference(fit)
    relabelled <- parameter_inference(shuffled)
    old_label <- function(outcome) {
      profile <- grepl("^profile_", outcome)
      outcome[profile] <- paste0("profile_",
                                 permutation[as.integer(sub("^profile_", "", outcome[profile]))])
      outcome
    }
    position <- match(paste(original$parameter, original$outcome, original$term),
                      paste(relabelled$parameter, old_label(relabelled$outcome),
                            relabelled$term))
    expect_false(anyNA(position))
    expect_equal(relabelled$standard_error[position], original$standard_error,
                 tolerance = 1e-6)
  }))
  # Equal means: only the variance features can tell the profiles apart.
  set.seed(9)
  state <- sample.int(2L, 300L, TRUE)
  spread <- data.frame(g = rep(1:50, each = 6),
                       a = stats::rnorm(300, 0, c(0.5, 3)[state]),
                       b = stats::rnorm(300, 0, c(0.5, 3)[state]))
  fit <- quietly(multilpa(spread, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          variance_model = "varying", n_starts = 3, seed = 1))
  signature <- .multilpa_profile_signature(fit)
  expect_true(all(is.finite(signature)))
  expect_identical(.multilpa_pool_align(.multilpa_permute_profiles(fit, 2:1), fit)$profile_order,
                   2:1)
})

test_that("standard errors are translation invariant and scale equivariant", {
  skip_on_cran()
  set.seed(17)
  n <- 240L
  state <- sample.int(2L, n, TRUE)
  data <- data.frame(g = rep(1:40, each = 6))
  data$a <- c(-1.5, 1.5)[state] + stats::rnorm(n)
  data$b <- c(1, -1)[state] + stats::rnorm(n)
  moved <- data
  moved$a <- moved$a * 1000 + 1e6
  fit_of <- function(frame) quietly(multilpa(frame, c("a", "b"), "g", n_profiles = 2,
                                             n_group_classes = 1, n_starts = 3, seed = 1,
                                             tol = 1e-11, variance_model = "varying",
                                             covariance_model = "full"))
  base <- parameter_inference(fit_of(data))
  scaled <- parameter_inference(fit_of(moved))
  power <- ifelse(base$parameter == "mean" & base$term %in% "a", 1,
                  ifelse(base$term %in% "a:a", 2, ifelse(base$term %in% "b:a", 1, 0)))
  expect_equal(scaled$standard_error / 1000^power, base$standard_error, tolerance = 1e-4)
})

test_that("a seeded bootstrap restores the caller's stream and is reproducible", {
  skip_on_cran()
  set.seed(1)
  data <- data.frame(g = rep(1:12, each = 6),
                     a = c(stats::rnorm(36, -1.5), stats::rnorm(36, 1.5)),
                     b = stats::rnorm(72))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          n_starts = 2, seed = 1))
  set.seed(99)
  expected <- stats::runif(3)
  set.seed(99)
  first <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3,
                                       n_starts = 1, seed = 7))
  expect_identical(stats::runif(3), expected)
  second <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3,
                                        n_starts = 1, seed = 7))
  expect_identical(attr(second, "replicates"), attr(first, "replicates"))
  set.seed(5)
  after_error <- stats::runif(1)
  set.seed(5)
  failing <- function() {
    .latents_local_seed(3)
    stop("boom")
  }
  expect_error(failing(), "boom")
  expect_identical(stats::runif(1), after_error)
})

test_that("starting variances keep profiles in rows with one indicator", {
  start <- list(covariances = array(c(2, 5), c(1, 1, 2)))
  expect_identical(.multilpa_start_variances(start), matrix(c(2, 5), 2L, 1L))
  start <- list(covariances = array(c(1, 0, 0, 3, 4, 0, 0, 6), c(2, 2, 2)))
  expect_identical(.multilpa_start_variances(start), matrix(c(1, 4, 3, 6), 2L, 2L))
})

# Regressions from the 2026-10-26 inference review (INF-1, INF-2).

test_that("covariate fits keep response standard errors for reserved category names", {
  set.seed(5)
  n <- 200L
  data <- data.frame(g = rep(1:40, each = 5), w = stats::rnorm(n))
  state <- stats::rbinom(n, 1, 0.5)
  data$a <- stats::rnorm(n, 2 * state)
  data$c <- factor(ifelse(stats::runif(n) < ifelse(state == 1, 0.7, 0.3), "x, y",
                          ifelse(stats::runif(n) < 0.5, "p:q", "z")))
  fit <- quietly(multilpa(data, c("a", "c"), "g", n_profiles = 2, n_group_classes = 1,
                          categorical = "c", profile_covariates = "w", n_starts = 3,
                          seed = 1))
  expect_false(anyNA(get_results(fit, "responses")$probability_standard_error))
})

test_that("a CSV round trip of a person-centred fit's data still aligns", {
  set.seed(1)
  data <- data.frame(g = rep(1:20, each = 5), a = stats::rnorm(100), b = stats::rnorm(100))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          centering = "person", n_starts = 1, seed = 1))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data, path, row.names = FALSE)
  expect_no_error(.multilpa_check_alignment(fit, utils::read.csv(path)))
})

test_that("a CSV round trip still refuses a changed person-centred value", {
  set.seed(1)
  data <- data.frame(g = rep(1:20, each = 5), a = stats::rnorm(100), b = stats::rnorm(100))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          centering = "person", n_starts = 1, seed = 1))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data, path, row.names = FALSE)
  changed <- utils::read.csv(path)
  changed$a[7L] <- changed$a[7L] + 1e-9
  expect_error(.multilpa_check_alignment(fit, changed),
               class = "latents_bad_inference_data")
  swapped <- utils::read.csv(path)
  swapped[c(1L, 2L), ] <- swapped[c(2L, 1L), ]
  expect_error(.multilpa_check_alignment(fit, swapped),
               class = "latents_bad_inference_data")
})
