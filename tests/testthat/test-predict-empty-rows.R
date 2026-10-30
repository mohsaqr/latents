test_that("profile predictions leave empty rows unclassified and keep zero-weight rows", {
  data <- .mixture_with_seed(702, {
    location <- rep(c(-3, 3), c(50, 100))
    data.frame(a = rnorm(150, location), b = rnorm(150, location), w = 1)
  })
  data[9, c("a", "b")] <- NA_real_
  data$a[10] <- NA_real_
  data$w[11] <- 0
  fit <- lpa(data, c("a", "b"), 2, weights = "w", missing = "fiml",
             n_starts = 1, seed = 1)
  raw <- fit$subject_posteriors
  fitted <- predict(fit, type = "posterior")
  fields <- setdiff(names(fitted), "row")
  expect_identical(fitted$row, seq_len(nrow(data)))
  expect_true(all(is.na(fitted[9, fields])))
  # A zero-weight unit is in no count but keeps its posterior in output (ML-1).
  expect_false(anyNA(fitted[c(10, 11), fields]))
  expect_identical(fit$subject_posteriors, raw)
  expect_equal(predict(fit, type = "density")$log_density[9], 0,
               tolerance = 1e-12)

  newcomers <- predict(fit, data, type = "posterior")
  expect_true(all(is.na(newcomers[9, fields])))
  expect_false(anyNA(newcomers[c(10, 11), fields]))
  expect_equal(newcomers[-9, ], fitted[-9, ], tolerance = 1e-12)

  fit$drop_empty <- FALSE
  legacy <- predict(fit, type = "posterior")
  expect_false(anyNA(legacy[c(9, 11), fields]))
  expect_false(anyNA(predict(fit, data, type = "posterior")))
})

test_that("group predictions distinguish empty observations from empty groups", {
  data <- .mixture_with_seed(703, {
    id <- rep(seq_len(30), each = 5)
    location <- ifelse(id <= 10, -3, 3)
    data.frame(id = id, a = rnorm(150, location), b = rnorm(150, location), w = 1)
  })
  data[1, c("a", "b")] <- NA_real_
  data[6:10, c("a", "b")] <- NA_real_
  data$w[11:15] <- 0
  fit <- multilpa(data, c("a", "b"), "id", 2, 2, weights = "w",
                  missing = "fiml", n_starts = 1, seed = 1)
  fitted <- predict(fit, type = "posterior")
  individual <- c("profile", "posterior", "probability_profile_1",
                  "probability_profile_2")
  group <- c("group_class", "group_posterior")
  expect_true(all(is.na(fitted[c(1, 6:10), individual])))
  expect_false(anyNA(fitted[1, group]))
  expect_true(all(is.na(fitted[6:10, group])))
  # The zero-weight group keeps its rows' and its own classification (ML-1).
  expect_false(anyNA(fitted[11:15, c(individual, group)]))
  expect_false(anyNA(fitted[2:5, c(individual, group)]))
  new <- predict(fit, data, type = "posterior")
  expect_true(all(is.na(new[6:10, c(individual, group)])))
  expect_false(anyNA(new[11:15, c(individual, group)]))

  fit$drop_empty <- FALSE
  legacy <- predict(fit, type = "posterior")
  expect_false(anyNA(legacy[, c(individual, group)]))
})

test_that("prediction recognizes observations of every indicator family", {
  data <- .mixture_with_seed(704, {
    location <- rep(c(-2, 2), each = 60)
    data.frame(y = rnorm(120, location), c = sample(c("a", "b"), 120, TRUE),
                o = sample(1:3, 120, TRUE), k = rpois(120, exp(location / 2)))
  })
  fit <- lpa(data, names(data), 2, categorical = "c", ordinal = "o",
             count = "k", missing = "fiml", n_starts = 1, seed = 1)
  new <- data[1:5, ]
  new[] <- NA
  new$y[2] <- data$y[2]
  new$c[3] <- "a"
  new$o[4] <- 2L
  new$k[5] <- 0L
  predicted <- predict(fit, new, type = "posterior")
  expect_true(all(is.na(predicted[1, -1])))
  expect_false(anyNA(predicted[2:5, ]))
})

test_that("new predictions classify a half-observed joined categorical pair", {
  data <- .mixture_with_seed(705, {
    data.frame(a = rbinom(120, 1, .3), b = rbinom(120, 1, .6),
                c = rbinom(120, 1, .5))
  })
  fit <- lca(data, names(data), 1, local_dependence = "a:b",
             missing = "fiml", n_starts = 1, seed = 1)
  new <- data.frame(a = c(1, NA, NA), b = c(NA, 0, NA), c = NA_real_)
  predicted <- predict(fit, new, type = "posterior")
  expect_false(anyNA(predicted[1:2, ]))
  expect_true(all(is.na(predicted[3, -1])))
})

test_that("empty rows also have missing noise membership predictions", {
  data <- .mixture_with_seed(706, data.frame(a = rnorm(100), b = rnorm(100)))
  fit <- lpa(data, names(data), 1, noise = TRUE, missing = "fiml",
             n_starts = 1, seed = 1)
  new <- data[1:2, ]
  new[1, ] <- NA_real_
  predicted <- predict(fit, new, type = "posterior")
  expect_true("probability_noise" %in% names(predicted))
  expect_true(all(is.na(predicted[1, -1])))
  expect_false(anyNA(predicted[2, ]))
  expect_equal(predict(fit, new, type = "density")$log_density[1], 0,
               tolerance = 1e-12)
  fit$drop_empty <- FALSE
  legacy <- predict(fit, new, type = "posterior")
  expect_equal(legacy$probability_noise[1], fit$noise_probability,
               tolerance = 1e-12)
  expect_equal(legacy$probability_profile_1[1],
               unname(fit$profile_probabilities[1, 1]), tolerance = 1e-12)
  new$a[1] <- 0
  expect_error(predict(fit, new, type = "posterior"),
                class = "latents_unsupported_noise")
  expect_error(predict(fit, new, type = "density"),
                class = "latents_unsupported_noise")
  training <- data
  training$a[1] <- NA_real_
  expect_error(lpa(training, names(training), 1, noise = TRUE,
                   missing = "fiml", n_starts = 1, seed = 1),
                class = "latents_unsupported_noise")
})
