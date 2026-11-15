# parameter_inference() for growth mixture models: the Wald table agrees with
# the fit's own tables, class matching undoes a relabelling, and the person
# bootstrap is reproducible and refuses what it does not cover.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

inference_fixture <- function(n = 160L, seed = 3L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:4
  x <- stats::rnorm(n)
  class <- 1L + (stats::runif(n) < stats::plogis(0.2 + 0.8 * x))
  out <- data.frame(id = rep(seq_len(n), each = length(times)), time = rep(times, n),
                    x = rep(x, each = length(times)))
  out$y <- c(2, 5)[class[out$id]] + c(0.6, -0.4)[class[out$id]] * out$time +
    rep(stats::rnorm(n), each = length(times)) +
    rep(stats::rnorm(n, 0, 0.3), each = length(times)) * out$time +
    stats::rnorm(nrow(out), 0, 0.7)
  out
}

inference_fit <- function(data, ...) {
  mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                     random = ~ 1 + time, membership = ~ x, n_starts = 2, seed = 1, ...)
}

test_that("the Wald table reproduces the fit's coefficient, random and class tables", {
  data <- inference_fixture()
  fit <- inference_fit(data)
  table <- parameter_inference(fit)
  expect_named(table, c("block", "class", "term", "estimate", "standard_error",
                        "statistic", "p_value", "p_adjusted", "conf_low", "conf_high"))
  coefficients <- as.data.frame(get_results(fit, "coefficients"))
  rows <- subset(table, block == "coefficient")
  expect_equal(rows$estimate, coefficients$estimate, tolerance = 1e-12)
  # Both are delta-method errors of the same linear map taken by central
  # differences with h = 1e-5 of a standard error; their arithmetic differs
  # only by rounding, eps / h, about 1e-10.
  expect_equal(rows$standard_error, coefficients$std_error, tolerance = 1e-8)
  expect_equal(rows$p_value, coefficients$p_value, tolerance = 1e-10)
  random <- as.data.frame(get_results(fit, "random"))
  for (kind in c("variance", "sd", "correlation")) {
    ours <- subset(table, block == paste0("random_", kind))
    theirs <- subset(random, parameter == kind)
    expect_equal(ours$estimate, theirs$estimate, tolerance = 1e-12)
    expect_equal(ours$standard_error, theirs$std_error, tolerance = 1e-10)
    expect_equal(ours$conf_low, theirs$conf_low, tolerance = 1e-10)
  }
  classes <- as.data.frame(get_results(fit, "classes"))
  expect_equal(subset(table, block == "share")$estimate, classes$share, tolerance = 1e-12)
  expect_equal(subset(table, block == "share")$standard_error, classes$share_std_error,
               tolerance = 1e-10)
  expect_equal(subset(table, block == "residual_sd")$estimate, classes$residual_sd,
               tolerance = 1e-12)
  membership <- as.data.frame(get_results(fit, "membership"))
  expect_equal(subset(table, block == "membership")$estimate, membership$estimate,
               tolerance = 1e-12)
  # Only coefficients and logits are tested; the family is those tests.
  expect_true(all(is.na(subset(table, !block %in% c("coefficient", "membership"))$p_value)))
  expect_equal(table$p_adjusted, table$p_value)
  adjusted <- parameter_inference(fit, adjust = "BH")
  expect_equal(adjusted$p_adjusted, stats::p.adjust(table$p_value, "BH"))
  expect_identical(attr(adjusted, "adjust"), "BH")
  expect_identical(attr(table, "method"), "wald")
})

test_that("relabelled classes are matched back to the original", {
  data <- inference_fixture()
  fit <- inference_fit(data, vcov_type = "none")
  natural <- latents:::.growth_natural_parameters(fit)
  swapped <- fit
  swapped$params <- latents:::.growth_permute_params(fit$spec, fit$params, c(2L, 1L))
  # A relabelling is the same model: the same likelihood, other names.
  expect_equal(latents:::.growth_expectation(fit$spec, fit$stats, swapped$params)$log_likelihood,
               fit$log_likelihood, tolerance = 1e-9)
  expect_false(isTRUE(all.equal(natural$values(swapped$params, fit$spec),
                                natural$values(fit$params, fit$spec))))
  aligned <- latents:::.growth_align_refit(fit, swapped)
  expect_equal(natural$values(aligned, fit$spec), natural$values(fit$params, fit$spec),
               tolerance = 1e-10)
})

test_that("the person bootstrap is reproducible and keeps the caller's stream", {
  skip_on_cran()
  data <- inference_fixture(n = 100L)
  fit <- inference_fit(data, vcov_type = "none")
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  one <- parameter_inference(fit, method = "bootstrap", iter = 4, n_starts = 1, seed = 5)
  expect_identical(stats::runif(1), before)
  two <- parameter_inference(fit, method = "bootstrap", iter = 4, n_starts = 1, seed = 5)
  expect_equal(one, two)
  expect_identical(attr(one, "method"), "bootstrap")
  expect_identical(attr(one, "n_valid"), 4L)
  expect_identical(dim(attr(one, "replicates")), c(4L, nrow(one)))
  expect_true(all(is.na(one$p_value)))
  expect_true(all(one$conf_low <= one$estimate + 1e-8 | one$conf_high >= one$estimate - 1e-8))
})

test_that("growth inference refuses what it does not cover", {
  data <- inference_fixture(n = 80L)
  fit <- inference_fit(data, vcov_type = "none")
  expect_error(parameter_inference(fit, data = subset(data, id > 2)),
               class = "latents_bad_inference_data")
  unconverged <- fit
  unconverged$converged <- FALSE
  expect_error(parameter_inference(unconverged), class = "latents_no_converge")
  degenerate <- fit
  degenerate$degenerate_classes <- "class 1"
  expect_error(parameter_inference(degenerate), class = "latents_boundary_fit")
  fixed <- parameter_inference(degenerate, boundary = "fix")
  held <- subset(fixed, startsWith(block, "random_") & class == "class_1")
  expect_true(all(is.na(held$standard_error)))
  expect_true(all(is.finite(subset(fixed, block == "coefficient")$standard_error)))
})

multilevel_inference_fixture <- function(n_schools = 14L, per_school = 8L, seed = 3L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_schools), function(s) {
    improving <- if (s %% 2L == 1L) 0.85 else 0.2
    do.call(rbind, lapply(seq_len(per_school), function(i) {
      k <- if (stats::runif(1) < improving) 1L else 2L
      data.frame(school = s, student = (s - 1) * per_school + i, wave = 0:3,
                 score = (if (k == 1L) 50 + 3 * 0:3 else rep(55, 4)) +
                   stats::rnorm(1, 0, 2) + stats::rnorm(4, 0, 1.5))
    }))
  }))
}

multilevel_inference_fit <- function(data, ...) {
  mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                     class_level = "group", random = ~ 1, n_group_classes = 2,
                     cluster = "school", n_starts = 2, seed = 1, ...)
}

test_that("a multilevel growth mixture reports its group-class shares", {
  skip_on_cran()
  fit <- multilevel_inference_fit(multilevel_inference_fixture())
  table <- parameter_inference(fit, boundary = "fix")
  expect_true("group_membership" %in% table$block)
  shares <- subset(table, block == "group_share")
  expect_identical(shares$class, c("group_class_1", "group_class_2"))
  # The group-class shares are the fit's own (its group_classes table), and
  # their interval is formed on the logit scale, so it stays inside (0, 1).
  reported <- get_results(fit, "group_classes")
  expect_equal(shares$estimate, as.vector(tapply(reported$group_share, reported$group_class,
                                                 `[`, 1L)), tolerance = 1e-10)
  expect_equal(sum(shares$estimate), 1, tolerance = 1e-12)
  expect_true(all(shares$conf_low > 0 & shares$conf_high < 1 &
                    is.finite(shares$standard_error)))
  expect_true(all(is.na(shares$p_value)))
  # One group class has no share to estimate: no group_share rows.
  single <- mixture_regression(score ~ wave, multilevel_inference_fixture(),
                               n_classes = 2, id = "student", class_level = "group",
                               random = ~ 1, n_group_classes = 1, cluster = "school",
                               n_starts = 1, seed = 1)
  expect_false("group_share" %in% parameter_inference(single, boundary = "fix")$block)
})

test_that("relabelled classes and group classes are matched back to the original", {
  skip_on_cran()
  fit <- multilevel_inference_fit(multilevel_inference_fixture(), vcov_type = "none")
  natural <- latents:::.growth_natural_parameters(fit)
  original <- natural$values(fit$params, fit$spec)
  # Every combination of a class and a group-class relabelling, including
  # one of each alone: a matching that confused the two levels, or matched
  # group classes through the trajectory labels, fails one of them.
  for (order in list(1:2, 2:1)) for (group_order in list(1:2, 2:1)) {
    swapped <- fit
    swapped$params <- latents:::.growth_permute_params(fit$spec, fit$params, order,
                                                       group_order)
    # A relabelling is the same model: the same likelihood, other names.
    expect_equal(latents:::.growth_expectation(fit$spec, fit$stats,
                                               swapped$params)$log_likelihood,
                 fit$log_likelihood, tolerance = 1e-9)
    matched <- latents:::.growth_match_classes(fit, swapped)
    expect_identical(matched$order, order)
    expect_identical(matched$group_order, group_order)
    if (!identical(order, 1:2) || !identical(group_order, 1:2)) {
      expect_false(isTRUE(all.equal(natural$values(swapped$params, fit$spec), original)))
    }
    aligned <- latents:::.growth_align_refit(fit, swapped)
    expect_equal(natural$values(aligned, fit$spec), original, tolerance = 1e-10)
  }
  # A group-class relabelling moves the group-class logits and the
  # trajectory-class intercepts of each group class, nothing else.
  moved <- latents:::.growth_permute_params(fit$spec, fit$params, 1:2, 2:1)
  expect_equal(moved$beta, fit$params$beta)
  expect_equal(moved$class_logits[1:2, ], fit$params$class_logits[2:1, ])
  expect_equal(moved$delta[, 2L], -fit$params$delta[, 2L])
  expect_error(latents:::.growth_permute_params(fit$spec, fit$params, 1:2, c(1L, 1L)),
               "permutation of the group classes")
})

test_that("a cluster resample keeps every drawn cluster's persons apart", {
  data <- data.frame(school = c(1, 1, 1, 2, 2, 2), student = c(1, 1, 2, 1, 1, 3),
                     score = 1:6)
  rows <- split(seq_len(nrow(data)), data$school)
  # School 1 drawn twice: its persons must become four persons, not two.
  resampled <- latents:::.multilpa_resample_groups(data, rows, "school", c(1L, 1L, 2L))
  renumbered <- latents:::.growth_renumber_persons(resampled, "school", "student")
  expect_identical(renumbered$school, rep(1:3, c(3L, 3L, 3L)))
  expect_identical(renumbered$student, c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 5L, 6L))
  expect_identical(renumbered$score, c(1:3, 1:3, 4:6))
  persons <- tapply(renumbered$school, renumbered$student, function(s) length(unique(s)))
  expect_true(all(persons == 1L))
})

test_that("the multilevel bootstrap resamples clusters and is reproducible", {
  skip_on_cran()
  fit <- multilevel_inference_fit(multilevel_inference_fixture(), vcov_type = "none")
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  one <- suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 4,
                                              n_starts = 1, seed = 5, boundary = "fix"))
  expect_identical(stats::runif(1), before)
  two <- suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 4,
                                              n_starts = 1, seed = 5, boundary = "fix"))
  expect_equal(one, two)
  wald <- parameter_inference(fit, boundary = "fix")
  expect_identical(as.data.frame(one)[c("block", "class", "term")],
                   as.data.frame(wald)[c("block", "class", "term")])
  expect_equal(one$estimate, wald$estimate, tolerance = 1e-10)
  expect_identical(attr(one, "method"), "bootstrap")
  expect_identical(dim(attr(one, "replicates")), c(attr(one, "n_valid"), nrow(one)))
  expect_gte(attr(one, "n_valid"), 2L)
  # Group-class shares of every replicate are shares.
  replicates <- attr(one, "replicates")
  shares <- replicates[, one$block == "group_share", drop = FALSE]
  expect_equal(rowSums(shares), rep(1, nrow(shares)), tolerance = 1e-10)
  # A replicate's units are clusters: the persons of a resample are those of
  # its clusters, so a refit sees 14 clusters whatever persons it draws.
  record <- new.env()
  record$seen <- integer()
  trace(latents:::.growth_renumber_persons, exit = bquote(
    assign("seen", c(get("seen", envir = .(record)), length(unique(returnValue()$school))),
           envir = .(record))), print = FALSE, where = asNamespace("latents"))
  on.exit(untrace(latents:::.growth_renumber_persons, where = asNamespace("latents")),
          add = TRUE)
  suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 2,
                                       n_starts = 1, seed = 6, boundary = "fix"))
  expect_identical(record$seen, c(14L, 14L))
})
