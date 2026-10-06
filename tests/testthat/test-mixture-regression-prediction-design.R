test_that("membership predictions retain training factor coding", {
  set.seed(83)
  data <- data.frame(x = rnorm(240), z = rep(c("a", "b", "c"), 80))
  state <- rbinom(240, 1, c(0.2, 0.5, 0.8)[match(data$z, c("a", "b", "c"))])
  data$y <- 4 * state + data$x + rnorm(240)
  fit <- mixture_regression(y ~ x, data, 2, membership = ~z,
                            n_starts = 1, seed = 1, vcov_type = "none")
  expected <- predict(fit, type = "class_response")
  reversed <- data
  reversed$z <- factor(reversed$z, levels = c("c", "b", "a"))
  expect_equal(predict(fit, reversed, type = "class_response"), expected)
  single <- predict(fit, data[1, , drop = FALSE], type = "class_response")
  expect_equal(single$prior, expected$prior[expected$row == 1])
  expect_equal(single$fitted, expected$fitted[expected$row == 1])
  options_before <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(options_before), add = TRUE)
  expect_equal(predict(fit, data, type = "class_response"), expected)
  bad <- data[1, , drop = FALSE]
  bad$z <- "unseen"
  expect_error(predict(fit, bad), "new level")
  bad$z <- NA_character_
  expect_error(predict(fit, bad), "missing values")
})

test_that("two-level membership designs retain factors and transformed terms", {
  set.seed(84)
  data <- data.frame(g = rep(seq_len(30), each = 6), x = rnorm(180),
                     z = rep(c("a", "b", "c"), 60))
  data$w <- rep(rep(c("u", "v", "w"), 10), each = 6)
  data$y <- rnorm(180)
  spec <- .mixture_spec(y ~ x, data, 2, "gaussian", "g", "observation", 2,
                        NULL, ~z + poly(x, 2), ~w, "varying", 1e-6, "error")
  object <- list(spec = spec)
  one <- data[1, , drop = FALSE]
  new <- .mixture_new_spec(object, one, FALSE)
  expect_equal(unname(new$w), unname(spec$w[1, , drop = FALSE]))
  expect_equal(unname(new$v), unname(spec$v[1, , drop = FALSE]))
  one$g <- NA_integer_
  expect_error(.mixture_new_spec(object, one, FALSE), class = "latents_bad_data")
})

test_that("regression plots carry membership columns and raw transformed predictors", {
  skip_on_cran()
  set.seed(85)
  data <- data.frame(g = rep(1:30, each = 4), x = runif(120, 1, 5),
                     z = rep(c("a", "b"), each = 60))
  data$y <- log(data$x) + rnorm(120)
  fit <- mixture_regression(y ~ log(x), data, 1, id = "g", class_level = "group",
                            membership = ~z, n_starts = 1, seed = 1,
                            vcov_type = "none")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, predictor = "x"))
  row_fit <- mixture_regression(y ~ x, data, 1, id = "g",
                                n_starts = 1, seed = 1, vcov_type = "none")
  subset <- data[c(9, 1), ]
  expect_identical(predict(row_fit, subset, type = "posterior")$g, subset$g)
})
