# Public reproductions and independent references for the extension audit.
# The comparisons with independent references are in tests/equivalence/test-audit-extension-repairs.R.

test_that("extreme-response regression posteriors retain the prior at crossing means", {
  set.seed(551)
  data <- data.frame(x = runif(240, -2, 2))
  classes <- rbinom(240, 1, .35)
  data$y <- ifelse(classes == 1, 3 + data$x, -3 - data$x) + rnorm(240, sd = .7)
  fit <- mixture_regression(y ~ x, data, 2, variance = "equal",
                            n_starts = 3, seed = 3, vcov_type = "none")
  beta <- fit$params$beta
  newdata <- data.frame(x = -diff(beta[1L, ]) / diff(beta[2L, ]), y = 1e8)
  posterior <- predict(fit, newdata, type = "posterior")
  probabilities <- as.matrix(posterior[c("probability_class_1", "probability_class_2")])
  prior <- plogis(fit$params$gamma[1L, 2L])
  expect_true(fit$converged)
  expect_equal(as.numeric(probabilities), c(1 - prior, prior), tolerance = 1e-12)
  expect_equal(unname(rowSums(probabilities)), 1, tolerance = 1e-15)
})

test_that("growth covariance caching does not merge distinct Gram matrices", {
  spec <- list(n = 4L, x = matrix(1, 4L, 1L), z = matrix(numeric(), 4L, 0L),
    random_design = cbind(1, c(1, 1, 1 - 1e-6, 1 + 1e-6)),
    y = c(-1, 1, -1, 1), offset = rep(0, 4L), group_index = rep(1:2, each = 2))
  statistics <- .growth_statistics(spec)
  expect_identical(statistics$pattern, 1:2)
  together <- .growth_class_terms(statistics, 0, diag(1e6, 2L), 1)$log_density
  # Computing each person alone supplies no opportunity for cache aliasing.
  separate <- vapply(1:2, function(person) {
    at <- which(spec$group_index == person)
    one <- spec
    one$n <- length(at)
    one$x <- spec$x[at, , drop = FALSE]
    one$z <- spec$z[at, , drop = FALSE]
    one$random_design <- spec$random_design[at, , drop = FALSE]
    one$y <- spec$y[at]
    one$offset <- spec$offset[at]
    one$group_index <- rep(1L, length(at))
    .growth_class_terms(.growth_statistics(one), 0, diag(1e6, 2L), 1)$log_density
  }, numeric(1))
  expect_equal(together, separate, tolerance = 1e-12)
  expect_gt(abs(diff(together)), .1)
})

test_that("extreme-response growth posteriors retain the prior at crossing means", {
  set.seed(512)
  data <- data.frame(id = rep(1:65, each = 4), time = rep(0:3, 65))
  classes <- rbinom(65, 1, .3)
  data$y <- ifelse(classes[data$id] == 1, 3 + data$time, -3 - data$time) +
    rep(rnorm(65, sd = .7), each = 4) + rnorm(260, sd = .5)
  fit <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
    random = ~ 1, variance = "equal", random_covariance = "equal",
    n_starts = 2, seed = 3, vcov_type = "none")
  beta <- fit$params$beta
  newdata <- data.frame(id = 100, time = -diff(beta[1L, ]) / diff(beta[2L, ]), y = 1e8)
  posterior <- predict(fit, newdata, type = "posterior")
  prior <- plogis(fit$params$gamma[1L, 2L])
  expect_true(fit$converged)
  expect_equal(unname(unlist(posterior[c("posterior_class_1", "posterior_class_2")])),
    c(1 - prior, prior), tolerance = 1e-12)
})

test_that("Brownian predictions refuse negative process times", {
  set.seed(161)
  data <- data.frame(id = rep(1:40, each = 5), time = rep(0:4, 40))
  process <- unlist(lapply(1:40, function(i) c(0, cumsum(rnorm(4, sd = .6)))))
  data$y <- 1 + .4 * data$time + rep(rnorm(40), each = 5) + process + rnorm(200, sd = .4)
  fit <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
    random = ~ 1, residual_correlation = "brownian", n_starts = 1, seed = 1,
    vcov_type = "none")
  newdata <- data[1:5, ]
  expect_true(fit$converged)
  expect_true(all(is.finite(predict(fit, newdata, type = "conditional")$fitted)))
  newdata$time <- -1:3
  expect_error(predict(fit, newdata, type = "conditional"), class = "latents_bad_data")
  expect_error(predict(fit, newdata), class = "latents_bad_data")
})

test_that("fitted class curves can vary a person-level membership predictor", {
  skip_if_not_installed("ggplot2")
  set.seed(91)
  data <- data.frame(id = rep(1:30, each = 4), time = rep(0:3, 30))
  data$z <- rep(rnorm(30), each = 4)
  data$y <- 1 + .4 * data$time + 2 * data$z + rep(rnorm(30), each = 4) + rnorm(120)
  fit <- mixture_regression(y ~ z + time, data, 1, id = "id", class_level = "group",
    membership = ~ z, n_starts = 0, vcov_type = "none")
  figure <- plot(fit, what = "fitted", predictor = "z")
  curves <- figure$layers[[2L]]$data
  expected <- fit$params$beta[1L, 1L] + fit$params$beta[2L, 1L] * curves$predictor_value +
    fit$params$beta[3L, 1L] * mean(data$time)
  expect_s3_class(figure, "ggplot")
  expect_equal(curves$outcome, as.vector(expected), tolerance = 1e-12)
  expect_equal(nrow(ggplot2::ggplot_build(figure)$data[[2L]]), 100L)
})

test_that("ordered random-effect contrasts are retained in spread and new predictions", {
  skip_if_not_installed("ggplot2")
  set.seed(541)
  data <- data.frame(id = rep(1:70, each = 6), time = rep(0:5, 70),
    condition = ordered(rep(letters[1:3], 140)))
  design <- model.matrix(~ condition, data)
  effects <- matrix(rnorm(70 * 3), 70, 3) * rep(c(1, .7, .5), each = 70)
  data$y <- 1 + .4 * data$time + rowSums(design * effects[data$id, ]) + rnorm(420, sd = .5)
  fit <- mixture_regression(y ~ time + condition, data, 1, id = "id",
    class_level = "group", random = ~ condition, n_starts = 1, seed = 1,
    vcov_type = "none")
  expect_true(fit$converged)
  expect_equal(.growth_random_design(fit$spec, data), fit$spec$random_design)
  typical <- .mixture_typical_rows(data[c("time", "condition")], "time", 0:5)
  expect_true(is.ordered(typical$condition))
  random <- model.matrix(~ condition, typical,
    contrasts.arg = attr(fit$spec$random_design, "contrasts"))
  variance <- rowSums((random %*% fit$params$random_covariance[[1L]]) * random)
  bands <- data.frame(class = "class_1", trajectory_time = 0:5, estimate = 0)
  spread <- .growth_spread(fit, bands, "time", .8)
  expect_equal(spread$spread_high, unname(qnorm(.9) * sqrt(variance)), tolerance = 1e-12)
  expect_equal(predict(fit, data, type = "conditional"), predict(fit, type = "conditional"),
    tolerance = 1e-10)
})

test_that("native numeric person IDs with identical display remain separate facets", {
  skip_if_not_installed("ggplot2")
  set.seed(42)
  ids <- c(.3, .1 + .2, 1:18)
  data <- data.frame(id = rep(ids, each = 4), time = rep(0:3, 20))
  data$y <- 1 + .5 * data$time + rep(rnorm(20), each = 4) + rnorm(80, sd = .4)
  fit <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
    random = ~ 1, n_starts = 1, seed = 1, vcov_type = "none")
  figure <- plot(fit, what = "individuals", persons = ids[1:2], time = "time")
  expect_equal(length(unique(figure$data$person)), 2L)
  expect_equal(nlevels(figure$data$person_label), 2L)
  expect_equal(nrow(ggplot2::ggplot_build(figure)$layout$layout), 2L)
})

test_that("trajectory data validation checks fitted values and permits extra distal columns", {
  set.seed(76)
  data <- data.frame(id = rep(1:25, each = 4), time = rep(0:3, 25))
  data$y <- 1 + .2 * data$time + rep(rnorm(25), each = 4) + rnorm(100, sd = .5)
  fit <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
    random = ~ 1, n_starts = 1, seed = 1, vcov_type = "none")
  expect_null(.trajectory_check_data(fit, transform(data, distal = 1)))
  changed <- data
  changed$y[1L] <- changed$y[1L] + .01
  expect_error(.trajectory_check_data(fit, changed), class = "latents_bad_inference_data")
  expect_error(parameter_inference(fit, changed), class = "latents_bad_inference_data")
  data$weight <- rep(1:25, each = 4)
  weighted <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
    random = ~ 1, weights = "weight", n_starts = 1, seed = 1, vcov_type = "none")
  expect_null(.trajectory_check_data(weighted, transform(data, weight = weight * 2)))
  changed <- data
  changed$weight[1:4] <- changed$weight[1:4] * 3
  expect_error(.trajectory_check_data(weighted, changed), class = "latents_bad_inference_data")
})

test_that("transformed random designs use their training prediction basis", {
  training <- data.frame(time = 0:9)
  random <- ~ poly(time, 2)
  frame <- model.frame(random, training)
  spec <- list(random = random, random_terms = terms(frame),
    random_design = model.matrix(random, frame), model_data = training)
  newdata <- data.frame(time = c(.5, 2.5, 4.5, 6.5))
  expected <- cbind(1, predict(poly(training$time, 2), newdata$time))
  expect_equal(unname(.growth_random_design(spec, newdata)), unname(expected),
    ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("plot interval matching separates term identity from profile identity", {
  errors <- data.frame(parameter = "response", outcome = c("profile_1", "profile_1"),
    term = c("x\ry:z", "x:y\rz"), conf_low = c(.1, .3), conf_high = c(.2, .4))
  actual <- .gg_match_bounds(errors, "response", c(1, 1), rev(errors$term))
  expect_equal(actual$lower, c(.3, .1))
  expect_equal(actual$upper, c(.4, .2))
})

test_that("categorical response plots use collision-safe inference intervals", {
  skip_if_not_installed("ggplot2")
  set.seed(482)
  data <- data.frame(check.names = FALSE,
    `a:b` = sample(c("c", "other"), 240, TRUE),
    a = sample(c("b:c", "other"), 240, TRUE))
  fit <- lca(data, c("a:b", "a"), 1, n_starts = 1, seed = 1)
  inference <- parameter_inference(fit)
  response <- subset(inference, parameter == "response")
  expect_false(anyDuplicated(response$term) > 0L)
  figure <- plot(fit, what = "responses", category = "first")
  terms <- c("a%3Ab:c", "a:b%3Ac")
  expected <- response[match(terms, response$term), ]
  expect_equal(figure$data$lower, expected$conf_low, tolerance = 1e-12)
  expect_equal(figure$data$upper, expected$conf_high, tolerance = 1e-12)
  expect_true(all(is.finite(figure$data$lower)))
})

test_that("seed scoping preserves the public RNG state without changing the normal kind", {
  data <- data.frame(x = 1:8, y = c(2, 1, 4, 3, 7, 5, 8, 9))
  fit <- mixture_regression(y ~ x, data, 1, n_starts = 0, vcov_type = "none")
  old_kind <- RNGkind()
  old_seed <- .Random.seed
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  RNGkind(normal.kind = "Box-Muller")
  set.seed(183)
  invisible(rnorm(1)) # Leave a cached normal variate that .Random.seed omits.
  expected_cached <- rnorm(1)
  set.seed(183)
  invisible(rnorm(1))
  seed <- .Random.seed
  kind <- RNGkind()
  scoped <- function() simulate(fit, nsim = 1, seed = 19)
  first <- scoped()
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  expect_identical(scoped(), first)
  # The uniform seed is restored; R's inaccessible cached normal is not.
  expect_false(identical(rnorm(1), expected_cached))
})

test_that("ordinal separation diagnostics use logits rather than raw slopes", {
  set.seed(819)
  data <- data.frame(x = rnorm(240))
  latent <- .3 + .8 * data$x + rlogis(240)
  data$y <- cut(latent, c(-Inf, -.6, .6, Inf), ordered_result = TRUE)
  ordinary <- mixture_regression(y ~ x, data, 1, family = "ordinal",
    n_starts = 0, vcov_type = "none")
  data$x <- data$x * .01
  expect_no_warning(scaled <- mixture_regression(y ~ x, data, 1, family = "ordinal",
    n_starts = 0, vcov_type = "none"))
  expect_gt(abs(scaled$params$beta[1L, 1L]), 30)
  expect_false(scaled$params$separation)
  expect_equal(scaled$log_likelihood, ordinary$log_likelihood, tolerance = 1e-8)
})
