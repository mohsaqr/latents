scaling_data <- function() {
  set.seed(42)
  group <- rep(seq_len(30L), each = 6L)
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = if (j %% 2L == 0L) c(0.8, 0.2) else c(0.25, 0.75))
  }, integer(1))
  # Two indicators on very different scales, so a scaling error cannot hide.
  data.frame(g = group,
             small = stats::rnorm(length(group), c(0, 2)[profile], 1),
             large = stats::rnorm(length(group), c(500, 900)[profile], 120))
}

fit_scaled <- function(data, ...) {
  multilpa(data, c("small", "large"), "g", n_profiles = 2, n_group_classes = 2,
           n_starts = 3, seed = 7, tol = 1e-10, ...)
}

test_that("grand scaling is the fit to grand z-scores, and differs only in units", {
  data <- scaling_data()
  raw <- fit_scaled(data, model = "VVI")
  scaled <- fit_scaled(data, model = "VVI", centering = "grand", scaling = "grand")
  # The same model fitted to columns standardized by hand.
  by_hand <- data
  by_hand[c("small", "large")] <- lapply(data[c("small", "large")],
                                         \(v) (v - mean(v)) / stats::sd(v))
  manual <- fit_scaled(by_hand, model = "VVI")
  expect_equal(scaled$log_likelihood, manual$log_likelihood, tolerance = 1e-8)
  expect_equal(unname(scaled$means), unname(manual$means), tolerance = 1e-6)
  # Change of variables: the scaled likelihood plus n * sum(log sd) is the raw.
  sds <- vapply(data[c("small", "large")], stats::sd, numeric(1))
  expect_equal(scaled$log_likelihood - nrow(data) * sum(log(sds)),
               raw$log_likelihood, tolerance = 1e-8)
  expect_identical(max.col(scaled$subject_posteriors),
                   max.col(raw$subject_posteriors))
  expected <- sweep(sweep(raw$means, 2L, colMeans(data[c("small", "large")])),
                    2L, sds, "/")
  expect_equal(unname(scaled$means), unname(expected), tolerance = 1e-5)
  # The fit keeps what it divided by, and gives back the columns supplied.
  expect_identical(scaled$scaling, "grand")
  expect_equal(unname(scaled$scaling_factors[1L, ]), unname(sds), tolerance = 1e-12)
  expect_equal(get_results(scaled, "data")[c("small", "large")],
               data[c("small", "large")], tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("person scaling is the fit to within-person z-scores", {
  data <- scaling_data()
  scaled <- fit_scaled(data, model = "VVI", centering = "person",
                       scaling = "person")
  by_hand <- data
  by_hand[c("small", "large")] <- lapply(data[c("small", "large")], \(v)
    stats::ave(v, data$g, FUN = \(u) (u - mean(u)) / stats::sd(u)))
  manual <- fit_scaled(by_hand, model = "VVI")
  expect_equal(scaled$log_likelihood, manual$log_likelihood, tolerance = 1e-8)
  expect_equal(unname(scaled$means), unname(manual$means), tolerance = 1e-6)
  # Every unit's transformed values have mean zero and spread one.
  stored <- as.data.frame(scaled$indicator_data)
  spread <- aggregate(cbind(small, large) ~ g,
                      data = cbind(stored, g = data$g), FUN = stats::sd)
  expect_equal(unname(as.matrix(spread[c("small", "large")])),
               matrix(1, 30L, 2L), tolerance = 1e-12)
  expect_equal(get_results(scaled, "data")[c("small", "large")],
               data[c("small", "large")], tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("scaling changes a spherical solution, as a hand-standardized fit does", {
  data <- scaling_data()
  raw <- fit_scaled(data, model = "EII")
  scaled <- fit_scaled(data, model = "EII", scaling = "grand")
  by_hand <- data
  by_hand[c("small", "large")] <- lapply(data[c("small", "large")],
                                         \(v) v / stats::sd(v))
  manual <- fit_scaled(by_hand, model = "EII")
  expect_equal(scaled$log_likelihood, manual$log_likelihood, tolerance = 1e-8)
  # One variance for both indicators: the raw fit is dominated by `large`.
  sds <- vapply(data[c("small", "large")], stats::sd, numeric(1))
  expect_false(isTRUE(all.equal(scaled$log_likelihood - nrow(data) * sum(log(sds)),
                                raw$log_likelihood, tolerance = 1e-6)))
})

test_that("predict() scales new data as the fit was scaled", {
  data <- scaling_data()
  grand <- fit_scaled(data, model = "VVI", centering = "grand", scaling = "grand")
  person <- fit_scaled(data, model = "VVI", centering = "person",
                       scaling = "person")
  # The training rows passed as new data come back with the fit's posteriors.
  lapply(list(grand, person), \(fit) {
    again <- predict(fit, newdata = data, type = "posterior")
    expect_equal(again$probability_profile_1,
                 unname(fit$subject_posteriors[, 1L]), tolerance = 1e-8)
  })
})

test_that("an enumeration passes scaling to every candidate", {
  data <- scaling_data()
  grid <- enumerate_classes(data, c("small", "large"), "g", n_profiles = 2,
                            n_group_classes = 2, model = "VVI", n_starts = 3,
                            seed = 7, tol = 1e-10, centering = "grand",
                            scaling = "grand")
  direct <- fit_scaled(data, model = "VVI", centering = "grand", scaling = "grand")
  candidate <- get_candidate(grid, "profiles2_groups2")
  expect_identical(candidate$scaling, "grand")
  expect_equal(candidate$log_likelihood, direct$log_likelihood, tolerance = 1e-8)
})

test_that("an undefined standard deviation is refused by class", {
  data <- scaling_data()
  # One group observed once: no within-person spread to divide by.
  single <- rbind(data, data.frame(g = 31L, small = 1, large = 700))
  expect_error(fit_scaled(single, model = "VVI", scaling = "person"),
               class = "latents_bad_data")
  # A group whose values never change on one indicator.
  flat <- data
  flat$small[flat$g == 3L] <- 0.5
  expect_error(fit_scaled(flat, model = "VVI", scaling = "person"),
               class = "latents_bad_data")
  constant <- data
  constant$small <- 2
  expect_error(fit_scaled(constant, model = "VVI", scaling = "grand"),
               class = "latents_bad_data")
  expect_error(fit_scaled(data, model = "VVI", scaling = "within"))
})

test_that("bootstrap tests refuse person scaling and mixed scaling", {
  data <- scaling_data()
  one <- multilpa(data, c("small", "large"), "g", n_profiles = 1,
                  n_group_classes = 1, n_starts = 1, seed = 1, scaling = "person",
                  centering = "person")
  two <- fit_scaled(data, model = "VVI", scaling = "person",
                    centering = "person")
  expect_error(bootstrap_lrt(one, two, data = data, iter = 2),
               class = "latents_unsupported_bootstrap")
  plain <- multilpa(data, c("small", "large"), "g", n_profiles = 1,
                    n_group_classes = 1, n_starts = 1, seed = 1)
  scaled <- fit_scaled(data, model = "VVI", scaling = "grand")
  expect_error(bootstrap_lrt(plain, scaled, data = data, iter = 2),
               class = "latents_incomparable_models")
})

test_that("scales do not depend on the order of the rows", {
  data <- scaling_data()
  shuffled <- data[sample(nrow(data)), ]
  x <- as.matrix(data[c("small", "large")])
  y <- as.matrix(shuffled[c("small", "large")])
  index_x <- match(data$g, unique(sort(data$g)))
  index_y <- match(shuffled$g, unique(sort(data$g)))
  expect_equal(.multilpa_indicator_scales(x, index_x, 30L, "person"),
               .multilpa_indicator_scales(y, index_y, 30L, "person"),
               tolerance = 1e-12)
})

test_that("a group's spread is judged on its own scale, and an empty group is left alone", {
  # Regression: the threshold once used the largest value of every group, so a
  # small but real spread beside a group near 1e9 was refused.
  x <- matrix(c(0:3, 1e9 + 100 * (0:3)), ncol = 1L, dimnames = list(NULL, "y"))
  groups <- rep(1:2, each = 4L)
  expect_equal(.multilpa_indicator_scales(x, groups, 2L, "person")[, "y"],
               c(stats::sd(0:3), stats::sd(100 * (0:3))), tolerance = 1e-10)
  x[1:4, 1L] <- NA
  expect_equal(unname(.multilpa_indicator_scales(x, groups, 2L, "person")[1L, "y"]), 1)
  flat <- matrix(c(5, 5, 5, 1, 2, 3), ncol = 1L, dimnames = list(NULL, "y"))
  expect_error(.multilpa_indicator_scales(flat, rep(1:2, each = 3L), 2L, "person"),
               class = "latents_bad_data")
})

test_that("a fit with no continuous indicator predicts under any centring and scaling", {
  # Regression: predict() applied grand offsets and scales the fit never had.
  set.seed(9)
  data <- data.frame(g = rep(1:10, each = 5), c = sample(c("a", "b", "c"), 50, TRUE))
  invisible(lapply(c("none", "grand", "person"), function(level) {
    fit <- multilpa(data, "c", "g", n_profiles = 2, n_group_classes = 1,
                    categorical = "c", centering = if (level == "person") "grand" else level,
                    scaling = level, n_starts = 2, seed = 1)
    expect_identical(nrow(predict(fit, data)), nrow(data))
  }))
})
