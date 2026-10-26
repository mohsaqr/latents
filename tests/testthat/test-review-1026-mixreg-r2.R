# Round-2 review of mixture_regression() (2026-10-26): N1 to N8.
# validation/review-2026-10-26/round-2/mixreg.md

r2_membership_data <- function(n = 400, seed = 14) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), t = stats::rnorm(n))
  cls <- stats::rbinom(n, 1, stats::plogis(0.3 + 1.2 * d$t)) + 1
  d$y <- ifelse(cls == 1, 1 + 2 * d$x, -1 - d$x) + stats::rnorm(n, sd = 0.6)
  d
}

test_that("a membership covariate far from its origin reaches the same maximum (N1)", {
  d <- r2_membership_data()
  fit_at <- function(shift) {
    shifted <- d
    shifted$m <- shift + d$t
    mixture_regression(y ~ x, shifted, 2, membership = "m", n_starts = 3, seed = 1)
  }
  base <- fit_at(0)
  base_membership <- get_results(base, "membership")
  for (shift in c(2010, 1e4)) {
    fit <- fit_at(shift)
    expect_true(fit$converged)
    expect_equal(fit$log_likelihood, base$log_likelihood, tolerance = 1e-9)
    membership <- get_results(fit, "membership")
    # The slope is origin-free; compare it in its standard errors.
    expect_lt(abs(membership$estimate[2] - base_membership$estimate[2]) /
                base_membership$std_error[2], 1e-6)
    expect_equal(membership$std_error[2], base_membership$std_error[2], tolerance = 1e-6)
  }
})

test_that("one membership M-step reaches the weighted logistic maximum at any origin (N1)", {
  set.seed(15)
  n <- 400
  t <- stats::rnorm(n)
  p2 <- stats::plogis(0.3 - 1.2 * t)
  counts <- cbind(1 - p2, p2) * stats::runif(n, 0.6, 1)
  counts <- counts / rowSums(counts)
  for (shift in c(0, 2010, 1e4)) {
    design <- cbind(1, shift + t)
    step <- .mixture_multinomial(design, counts, matrix(0, 2, 2))
    expect_equal(step$coefficients[2, 2], -1.2, tolerance = 1e-8)
  }
})

test_that("membership separation is read on fitted probabilities, not raw coefficients (N2)", {
  d <- r2_membership_data()
  for (setting in list(c(0, 0.01), c(2010, 6), c(30, 1))) {
    shifted <- d
    shifted$m <- setting[1] + setting[2] * d$t
    expect_no_warning(mixture_regression(y ~ x, shifted, 2, membership = "m",
                                         n_starts = 3, seed = 1),
                      class = "latents_separation")
  }
  # Genuine separation is still flagged: membership decides the class exactly.
  set.seed(3)
  groups <- 80
  g <- data.frame(id = rep(seq_len(groups), each = 6), x = stats::rnorm(groups * 6))
  side <- rep(c(-1, 1), length.out = groups)
  g$gm <- (side * stats::runif(groups, 0.5, 2))[g$id]
  g$y <- ifelse(g$gm > 0, 3 + 2 * g$x, -3 - g$x) + stats::rnorm(nrow(g), sd = 0.3)
  expect_warning(mixture_regression(y ~ x, g, 2, id = "id", class_level = "group",
                                    membership = "gm", n_starts = 3, seed = 1),
                 class = "latents_separation")
})

test_that("the bootstrap test gives a p-value with a calendar-year membership covariate (N2)", {
  skip_on_cran()
  d <- r2_membership_data(n = 200)
  p_value_at <- function(shift) {
    shifted <- d
    shifted$m <- shift + 2 * d$t
    grid <- suppressWarnings(enumerate_regressions(
      y ~ x, shifted, n_classes = 1:2, membership = "m", n_starts = 2, seed = 2,
      bootstrap = 3, bootstrap_starts = 2))
    as.data.frame(grid)$blrt_p_value[2]
  }
  at_2010 <- p_value_at(2010)
  expect_false(is.na(at_2010))
  expect_equal(at_2010, p_value_at(0))
})

test_that("the observation nesting counts rows, with or without weights and id (N3)", {
  set.seed(4)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), id = rep(1:60, each = 5))
  cls <- stats::rbinom(n, 1, 0.5) + 1
  d$y <- ifelse(cls == 1, 1 + 2 * d$x, -1 - d$x) + stats::rnorm(n, sd = 0.7)
  d$one <- 1
  plain <- mixture_regression(y ~ x, d, 2, id = "id", n_starts = 3, seed = 1,
                              vcov_type = "none")
  weighted <- mixture_regression(y ~ x, d, 2, id = "id", n_starts = 3, seed = 1,
                                 weights = "one", vcov_type = "none")
  a <- get_results(plain, "fit")
  b <- get_results(weighted, "fit")
  expect_equal(nobs(weighted), 300L)
  expect_equal(attr(logLik(weighted), "nobs"), attr(logLik(plain), "nobs"))
  expect_equal(b$bic, a$bic, tolerance = 1e-10)
  expect_equal(b$sabic, a$sabic, tolerance = 1e-10)
  expect_equal(b$icl, a$icl, tolerance = 1e-10)
  # A row of weight zero leaves the sample.
  d$w0 <- 1
  d$w0[d$id == 60] <- 0
  dropped <- mixture_regression(y ~ x, d, 2, id = "id", n_starts = 3, seed = 1,
                                weights = "w0", vcov_type = "none")
  expect_equal(nobs(dropped), 295L)
})

test_that("a negative-binomial dispersion at the Poisson limit lands on its floor and is held (N4)", {
  set.seed(102)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  cls <- stats::rbinom(n, 1, 0.5) + 1
  d$y <- stats::rpois(n, exp(ifelse(cls == 1, 0.3 + 0.8 * d$x, 1.5 - 0.4 * d$x)))
  fit <- NULL
  expect_warning(fit <- mixture_regression(y ~ x, d, 2, family = "negative_binomial",
                                           n_starts = 5, seed = 1, tol = 1e-12),
                 class = "latents_boundary")
  expect_true(all(fit$params$sigma2 <= .latents_min_dispersion * (1 + 1e-6)))
  classes <- get_results(fit, "classes")
  expect_true(all(is.na(classes$dispersion_std_error)))
  expect_true(all(is.finite(as.data.frame(fit)$std_error)))
  # The density near the limit is the exact Poisson-plus-expansion value.
  alpha <- matrix(1e-8 * (1 + (0:10) * 1e-5), 50, 11, byrow = TRUE)
  counts <- stats::rpois(50, 3)
  density <- .mixture_negative_binomial_log_density(counts, matrix(3, 50, 11), alpha)
  totals <- colSums(density)
  linear <- sum(stats::dpois(counts, 3, log = TRUE)) +
    alpha[1, ] * sum(((counts - 3)^2 - counts) / 2)
  expect_lt(max(abs(totals - linear)), 1e-12)
})

test_that("an ordinal class with an empty middle category is a boundary fit (N5)", {
  set.seed(24)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), m = stats::rnorm(n))
  cls <- stats::rbinom(n, 1, stats::plogis(d$m)) + 1
  d$count <- stats::rpois(n, exp(ifelse(cls == 1, 0.3 + 0.8 * d$x, 1.5 - 0.4 * d$x)))
  d$o <- cut(ifelse(cls == 1, 2 * d$x, -d$x) + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
             labels = FALSE)
  d$w <- sample(1:3, n, TRUE)
  long <- d[rep(seq_len(n), d$w), ]
  fit <- NULL
  expect_warning(fit <- mixture_regression(o ~ x, long, 2, family = "ordinal",
                                           membership = "m", n_starts = 4, seed = 1,
                                           tol = 1e-13, max_iter = 5000),
                 class = "latents_boundary")
  theta <- coef(fit)
  gap <- "log_threshold_gap.class_2.2|3"
  expect_lt(theta[[gap]], .ordinal_min_log_gap)
  expect_true(is.na(sqrt(vcov(fit)[gap, gap])))
  table <- as.data.frame(fit)
  class_2 <- table[table$class == "class_2", ]
  expect_false(is.na(class_2$std_error[class_2$term == "threshold:1|2"]))
  expect_true(all(is.na(class_2$std_error[class_2$term %in% c("threshold:2|3",
                                                              "threshold:3|4")])))
  expect_true(all(is.finite(class_2$std_error[class_2$term == "x"])))
})

test_that("standard errors do not depend on a predictor's origin (N6)", {
  set.seed(3)
  n <- 300
  d <- data.frame(t = stats::rnorm(n))
  d$y <- 1 + 2 * d$t + stats::rnorm(n)
  d$count <- stats::rpois(n, exp(0.5 + 0.3 * d$t))
  for (shift in c(1e4, 5e4, 1e6)) {
    shifted <- d
    shifted$x <- d$t + shift
    gaussian <- mixture_regression(y ~ x, shifted, 1, n_starts = 0)
    reference <- stats::glm(y ~ x, data = shifted)
    expect_equal(as.data.frame(gaussian)$std_error[2],
                 sqrt(stats::vcov(reference)[2, 2] * (n - 2) / n), tolerance = 1e-7)
    poisson <- mixture_regression(count ~ x, shifted, 1, family = "poisson", n_starts = 0)
    reference <- stats::glm(count ~ x, data = shifted, family = stats::poisson)
    expect_equal(as.data.frame(poisson)$std_error[2],
                 sqrt(stats::vcov(reference)[2, 2]), tolerance = 1e-6)
  }
  # Two classes: the slope's errors are the same at origin 0 and 1e5.
  set.seed(2)
  cls <- stats::rbinom(n, 1, 0.4) + 1
  d$z <- ifelse(cls == 1, 1 + 2 * d$t, -1 - d$t) + stats::rnorm(n, sd = 0.5)
  errors_at <- function(shift) {
    shifted <- d
    shifted$x <- d$t + shift
    fit <- mixture_regression(z ~ x, shifted, 2, n_starts = 4, seed = 1, tol = 1e-12)
    table <- as.data.frame(fit)
    table$std_error[table$term == "x"]
  }
  expect_equal(errors_at(1e5), errors_at(0), tolerance = 1e-8)
})

test_that("units of weight zero leave the classes and classification tables (N7)", {
  set.seed(19)
  groups <- 40
  d <- data.frame(id = rep(seq_len(groups), each = 5), x = stats::rnorm(groups * 5))
  gc <- (stats::rbinom(groups, 1, 0.5) + 1)[d$id]
  d$y <- ifelse(gc == 1, 1 + 2 * d$x, -1 - d$x) + stats::rnorm(nrow(d), sd = 0.7)
  d$w <- stats::runif(groups, 0.5, 2)[d$id]
  extra <- d[d$id <= 6, ]
  extra$id <- extra$id + 1000
  extra$w <- 0
  extra$y <- extra$y + 5
  with_zero <- rbind(d, extra)
  for (level in c("group", "observation")) {
    a <- mixture_regression(y ~ x, d, 2, id = "id", class_level = level,
                            n_starts = 3, seed = 1, weights = "w")
    b <- mixture_regression(y ~ x, with_zero, 2, id = "id", class_level = level,
                            n_starts = 3, seed = 1, weights = "w")
    expect_identical(get_results(b, "classes")$n_assigned,
                     get_results(a, "classes")$n_assigned)
    expect_equal(get_results(b, "classes")$mean_posterior,
                 get_results(a, "classes")$mean_posterior, tolerance = 1e-4)
    expect_equal(get_results(b, "classification")$mean_posterior,
                 get_results(a, "classification")$mean_posterior, tolerance = 1e-4)
  }
})

test_that("row classes of a two-level model are counted against rows (N8)", {
  set.seed(22)
  groups <- 4
  n <- 400
  d <- data.frame(id = rep(seq_len(groups), each = n / groups), x = stats::rnorm(n))
  cls <- sample(1:4, n, TRUE)
  d$y <- c(-3, -1, 1, 3)[cls] + d$x + stats::rnorm(n, sd = 0.3)
  fit <- mixture_regression(y ~ x, d, 4, id = "id", n_group_classes = 2,
                            n_starts = 2, seed = 1, vcov_type = "none")
  expect_s3_class(fit, "latents_mixture_regression")
  # Group classes still need more groups than group classes.
  expect_error(mixture_regression(y ~ x, d, 2, id = "id", n_group_classes = 4,
                                  n_starts = 1, seed = 1, vcov_type = "none"),
               class = "latents_bad_data")
})

test_that("the bootstrap p-value uses the valid replicates unless fewer than half are valid", {
  data <- data.frame(y = 1:5, x = 1:5)
  null <- list(log_likelihood = -10, converged = TRUE,
               spec = list(response_name = "y", model_data = data, formula = y ~ x))
  alternative <- null
  alternative$log_likelihood <- -8
  grid <- data.frame(n_classes = 1:2, n_group_classes = 1L)
  # Replicates whose number is in `failing` fail their two-class refit; the
  # others give the statistic 2 * (k - 1) * 0.5 = 1, below the observed 4.
  run_with <- function(failing) {
    calls <- 0L
    testthat::local_mocked_bindings(
      .mixture_draw = function(object) data$y,
      .mixture_try_fit = function(formula, data, n, ...) {
        if (n == 2L) {
          calls <<- calls + 1L
          if (calls %in% failing) return(errorCondition("refit failed"))
        }
        list(log_likelihood = -10 + 0.5 * (n - 1), converged = TRUE)
      }, .package = "latents")
    .mixture_blrt(list(null, alternative), grid, 2L, 3L, 0L, list())
  }
  one <- NULL
  expect_warning(one <- run_with(2L), class = "latents_failed_replicates")
  expect_equal(one$blrt_p_value, 1 / 3)
  expect_identical(one$blrt_replicates, 2L)
  expect_identical(one$blrt_n_invalid, 1L)
  expect_identical(one$blrt_flagged, 1L)
  two <- NULL
  expect_warning(two <- run_with(c(1L, 3L)), class = "latents_failed_replicates")
  expect_true(is.na(two$blrt_p_value))
  expect_identical(two$blrt_replicates, 1L)
  expect_identical(two$blrt_n_invalid, 2L)
  none <- NULL
  expect_no_warning(none <- run_with(integer()))
  expect_equal(none$blrt_p_value, 1 / 4)
  expect_identical(none$blrt_n_invalid, 0L)
})

test_that("the shared bootstrap p-value uses valid replicates unless fewer than half are valid", {
  replicates <- data.frame(valid = c(TRUE, FALSE, TRUE), statistic = c(1, NA, 5))
  one <- NULL
  expect_warning(one <- .latents_blrt_p_value(replicates, 4, 3L),
                 class = "latents_failed_replicates")
  expect_equal(one$p_value, (1 + 1) / (2 + 1))
  expect_equal(one$monte_carlo_se, sqrt(one$p_value * (1 - one$p_value) / 3))
  expect_identical(one$n_valid, 2L)
  expect_identical(one$n_invalid, 1L)
  replicates$valid[3] <- FALSE
  two <- NULL
  expect_warning(two <- .latents_blrt_p_value(replicates, 4, 3L),
                 class = "latents_failed_replicates")
  expect_identical(two$p_value, NA_real_)
  expect_identical(two$n_invalid, 2L)
})

# One invalid replicate of three, as the replicate loop records it, for the
# family verbs' wiring of the shared rule.
r2_replicates_one_invalid <- function(iter, replicate, muffle_warnings) {
  data.frame(replicate = seq_len(3L), statistic = c(0.5, NA, 100), valid = c(TRUE, FALSE, TRUE),
             boundary = FALSE, logits_settled = NA, null_replications = NA_integer_,
             alternative_replications = NA_integer_, warnings = NA_character_,
             error = c(NA, "Nonconvergence or reversed likelihood", NA))
}

r2_expect_bootstrap_rule <- function(result) {
  test <- get_results(result, "test")
  expect_identical(test$n_valid, 2L)
  expect_identical(test$n_invalid, 1L)
  expected <- (1 + (100 >= test$statistic) + (0.5 >= test$statistic)) / 3
  expect_equal(test$p_value, expected)
  expect_output(print(result), "1 invalid")
}

test_that("bootstrap_lrt() for profile fits takes the p-value over valid replicates", {
  set.seed(1)
  data <- data.frame(school = rep(seq_len(10), each = 10),
                     score = stats::rnorm(100, rep(c(-2, 2), each = 50)))
  smaller <- multilpa(data, "score", "school", n_profiles = 1, n_group_classes = 1,
                      n_starts = 2, seed = 1)
  larger <- multilpa(data, "score", "school", n_profiles = 2, n_group_classes = 1,
                     n_starts = 2, seed = 1)
  testthat::local_mocked_bindings(.latents_blrt_replicates = r2_replicates_one_invalid,
                                  .package = "latents")
  result <- NULL
  expect_warning(result <- bootstrap_lrt(smaller, larger, iter = 3L, seed = 1),
                 class = "latents_failed_replicates")
  r2_expect_bootstrap_rule(result)
})

test_that("bootstrap_lrt() for additive fits takes the p-value over valid replicates", {
  set.seed(31)
  groups <- 40L
  sizes <- sample(2:7, groups, replace = TRUE)
  classes <- sample(1:2, groups, replace = TRUE)
  data <- do.call(rbind, lapply(seq_len(groups), function(j) {
    centre <- stats::rnorm(2L, c(0, 2.5)[classes[j]] + c(0, 1), 0.8)
    data.frame(g = j, y1 = centre[1] + stats::rnorm(sizes[j]),
               y2 = centre[2] + stats::rnorm(sizes[j]))
  }))
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (!inherits(w, "latents_failed_replicates")) invokeRestart("muffleWarning")
  })
  fit <- quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2,
                        family = "additive", n_starts = 2, seed = 1))
  null_model <- quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1,
                               family = "additive", n_starts = 1, seed = 1))
  testthat::local_mocked_bindings(.latents_blrt_replicates = r2_replicates_one_invalid,
                                  .package = "latents")
  result <- NULL
  expect_warning(result <- quiet(bootstrap_lrt(null_model, fit, iter = 3L, seed = 1L)),
                 class = "latents_failed_replicates")
  r2_expect_bootstrap_rule(result)
})

test_that("bootstrap_lrt() for transition fits takes the p-value over valid replicates", {
  activity <- c("browse", "lectures")
  null <- lta(course_engagement, activity, "student", n_profiles = 2,
              time = "sequence", n_starts = 1, seed = 1)
  alternative <- lta(course_engagement, activity, "student", n_profiles = 2,
                     time = "sequence", transition_covariates = "previous_grade",
                     n_starts = 1, seed = 1)
  testthat::local_mocked_bindings(.latents_blrt_replicates = r2_replicates_one_invalid,
                                  .package = "latents")
  result <- NULL
  expect_warning(result <- bootstrap_lrt(null, alternative, data = course_engagement,
                                         iter = 3L, seed = 2),
                 class = "latents_failed_replicates")
  r2_expect_bootstrap_rule(result)
})
