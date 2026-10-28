# Release audit 2026-10-28: stable shared-dispersion names (REL-1), and
# zero-weight units absent from every inference design and score (REL-3).

test_that("shared regression dispersions retain public names and finite errors", {
  .latents_local_seed(2810L)
  groups <- 70L
  rows <- rep(seq_len(groups), each = 5L)
  class <- rep(rep_len(1:2, groups), each = 5L)
  data <- data.frame(id = rows, x = rnorm(length(rows)), r = rnorm(length(rows)))
  data$gaussian <- c(-3, 3)[class] + 0.4 * data$x + 0.3 * data$r +
    rnorm(length(rows), sd = 0.7)
  data$negative_binomial <- rnbinom(length(rows), size = 4,
    mu = exp(c(0.2, 2.5)[class] + 0.2 * data$x + 0.15 * data$r))
  invisible(lapply(c("gaussian", "negative_binomial"), function(family) {
    invisible(lapply(1:2, function(k) {
      fit <- mixture_regression(reformulate(c("x", "r"), family), data, k,
        family = family, id = "id", class_level = "group", common = ~r,
        variance = "equal", n_starts = 2L, seed = 1L, max_iter = 2000L,
        tol = 1e-10)
      key <- if (family == "gaussian") "log_sigma.all" else "log_dispersion.all"
      quantity <- if (family == "gaussian") "sigma" else "dispersion"
      classes <- get_results(fit, "classes")
      expect_true(fit$converged)
      expect_true(key %in% names(coef(fit)))
      expect_identical(names(coef(fit)), rownames(vcov(fit)))
      expect_false(any(grepl("[.]all[.]class_", names(coef(fit)))))
      expected <- if (family == "gaussian") sqrt(fit$params$sigma2) else fit$params$sigma2
      expect_equal(classes[[quantity]], unname(expected), tolerance = 1e-10)
      expect_true(all(is.finite(classes[[paste0(quantity, "_std_error")]])))
    }))
  }))
})

test_that("an excluded regression row cannot remove inference", {
  .latents_local_seed(72L)
  data <- data.frame(x = rnorm(100L), r = rnorm(100L), w = 1)
  data$y <- 2 + 1.5 * data$x + 0.4 * data$r + rnorm(100L)
  padded <- rbind(data, data.frame(x = 1e8, r = -1e8, w = 0, y = 1e8))
  fit_one <- function(frame) mixture_regression(y ~ x + r, frame, 1L,
    common = ~r, weights = "w", n_starts = 1L, seed = 2L)
  baseline <- fit_one(data)
  added <- fit_one(padded)
  expect_equal(coef(added), coef(baseline), tolerance = 1e-12)
  expect_equal(as.numeric(logLik(added)), as.numeric(logLik(baseline)), tolerance = 1e-12)
  expect_true(all(is.finite(vcov(added))))
  expect_equal(vcov(added), vcov(baseline), tolerance = 1e-10)
  expect_equal(get_results(added, "coefficients"), get_results(baseline, "coefficients"),
               tolerance = 1e-10)
  expect_identical(added$inference$clustered_on, "row")
  expect_equal(nrow(added$expectation$tau), nrow(padded))
})

test_that("excluded groups cannot choose regression inference coordinates", {
  .latents_local_seed(2813L)
  n_groups <- 90L
  index <- rep(seq_len(n_groups), each = 6L)
  g <- rnorm(n_groups)
  m <- rnorm(n_groups)
  group_class <- rbinom(n_groups, 1L, plogis(0.5 * g))
  class <- rbinom(length(index), 1L,
    plogis(ifelse(group_class[index] == 0L, -1.5, 1.5) + 0.3 * m[index]))
  data <- data.frame(id = index, x = rnorm(length(index)), r = rnorm(length(index)),
                     m = m[index], g = g[index],
                     weight = rep_len(c(0.6, 1, 1.4), n_groups)[index])
  data$y <- ifelse(class == 0L, -3, 3) + 0.4 * data$x + 0.2 * data$r +
    rnorm(nrow(data), sd = 0.7)
  sentinel <- data[seq_len(6L), ]
  # The excluded group precedes the carrying groups in their sorted index.
  sentinel$id <- 0L
  sentinel[c("x", "r", "m", "g", "y")] <- 1e8
  sentinel$weight <- 0
  padded <- rbind(data, sentinel)
  cases <- list(
    clustered = list(n_classes = 1L),
    group = list(n_classes = 2L, class_level = "group", membership = ~m),
    two_level = list(n_classes = 2L, n_group_classes = 2L,
                      membership = ~m, group_membership = ~g))
  invisible(lapply(cases, function(arguments) {
    fit_one <- function(frame) do.call(mixture_regression, c(list(
      formula = y ~ x + r, data = frame, common = ~r, id = "id", weights = "weight",
      n_starts = 2L, seed = 1L, tol = 1e-9, max_iter = 3000L), arguments))
    baseline <- fit_one(data)
    added <- fit_one(padded)
    expect_true(baseline$converged)
    expect_true(all(is.finite(vcov(baseline))))
    expect_equal(coef(added), coef(baseline), tolerance = 1e-12)
    expect_equal(as.numeric(logLik(added)), as.numeric(logLik(baseline)), tolerance = 1e-12)
    expect_true(all(is.finite(vcov(added))))
    expect_equal(vcov(added), vcov(baseline), tolerance = 1e-10)
    expect_identical(added$inference$clustered_on, "id")
    expect_equal(nrow(added$expectation$tau), nrow(padded))
    expect_equal(added$spec$n_groups, n_groups + 1L)
  }))
})
