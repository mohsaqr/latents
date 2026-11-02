skip_on_cran()

# Review of 2026-11-01, findings F2 and F9.
# F2: a bootstrap table's replicate columns are paired with its rows by the
# parameters' labels and estimates, in order, never by row names. A reordered
# copy with ordinary row names (base `[` then reset row names, or
# dplyr::arrange()) used to pass and attach replicates to the wrong labels.
# F9: whether a parameter varies across the replicates is decided relative to
# its own values, so data in small units, or far from the origin, keep every
# parameter that varies.

review_bootstrap <- local({
  made <- list()
  function(shift = 0, scale = 1) {
    key <- sprintf("%g/%g", shift, scale)
    if (!is.null(made[[key]])) return(made[[key]])
    data <- datasets::faithful
    data$eruptions <- shift + scale * data$eruptions
    data$waiting <- shift + scale * data$waiting
    fit <- suppressWarnings(lpa(data, c("eruptions", "waiting"), 2, model = "VVI",
                                n_starts = 2, seed = 1))
    boot <- suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 6,
                                                 n_starts = 1, seed = 1))
    made[[key]] <<- boot
    boot
  }
})

# The tidy replicates of `table`, one parameter per label key, must give back
# the table's standard errors and percentile bounds.
expect_reproduces_table <- function(table, labels) {
  replicates <- as.data.frame(table, what = "replicates")
  key <- function(frame) do.call(paste, c(lapply(labels, function(name) frame[[name]]),
                                          sep = "|"))
  rows <- key(table)
  expect_false(anyDuplicated(rows) > 0L)
  estimated <- !is.na(table$standard_error)
  spread <- vapply(rows, function(row) stats::sd(replicates$estimate[key(replicates) == row]),
                   numeric(1))
  expect_equal(unname(spread[estimated]), table$standard_error[estimated], tolerance = 1e-8)
  level <- attr(table, "level")
  bounds <- vapply(rows, function(row) {
    stats::quantile(replicates$estimate[key(replicates) == row],
                    c((1 - level) / 2, (1 + level) / 2), names = FALSE)
  }, numeric(2))
  expect_equal(unname(bounds[1L, estimated]), table$conf_low[estimated], tolerance = 1e-8)
  expect_equal(unname(bounds[2L, estimated]), table$conf_high[estimated], tolerance = 1e-8)
}

profile_labels <- c("level", "outcome", "term", "parameter")

test_that("an unmodified table, whatever its row names, reproduces its errors", {
  boot <- review_bootstrap()
  expect_reproduces_table(boot, profile_labels)
  renamed <- boot
  row.names(renamed) <- sprintf("p%02d", seq_len(nrow(boot)))
  expect_reproduces_table(renamed, profile_labels)
  same_order <- boot[seq_len(nrow(boot)), , drop = FALSE]
  row.names(same_order) <- NULL
  expect_reproduces_table(same_order, profile_labels)
  # A column added after the statistics changes no parameter's identity.
  annotated <- boot
  annotated$note <- "checked"
  expect_reproduces_table(annotated, profile_labels)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(renamed), "ggplot")
})

test_that("a reordered copy is refused even when its row names look original", {
  boot <- review_bootstrap()
  reversed <- boot[rev(seq_len(nrow(boot))), , drop = FALSE]
  expect_error(as.data.frame(reversed, what = "replicates"), class = "latents_bad_argument")
  row.names(reversed) <- NULL
  expect_identical(row.names(reversed), row.names(boot))
  expect_error(as.data.frame(reversed, what = "replicates"), class = "latents_bad_argument")
  skip_if_not_installed("ggplot2")
  expect_error(plot(reversed), class = "latents_bad_argument")
})

test_that("dplyr::arrange() keeps the row names but not the pairing, and is refused", {
  skip_if_not_installed("dplyr")
  boot <- review_bootstrap()
  arranged <- dplyr::arrange(boot, dplyr::desc(estimate))
  expect_identical(row.names(arranged), row.names(boot))
  expect_error(as.data.frame(arranged, what = "replicates"), class = "latents_bad_argument")
  skip_if_not_installed("ggplot2")
  expect_error(plot(arranged, what = "mean"), class = "latents_bad_argument")
})

test_that("repeated rows, changed labels and lost label columns are refused", {
  boot <- review_bootstrap()
  n <- nrow(boot)
  # Same number of rows, one parameter twice.
  repeated <- boot[c(1L, 1L, seq_len(n)[-c(1L, 2L)]), , drop = FALSE]
  row.names(repeated) <- NULL
  expect_error(as.data.frame(repeated, what = "replicates"), class = "latents_bad_argument")
  appended <- boot[c(seq_len(n), 1L), , drop = FALSE]
  expect_error(as.data.frame(appended, what = "replicates"), class = "latents_bad_argument")
  relabelled <- boot
  relabelled$outcome[1L] <- "profile_9"
  expect_error(as.data.frame(relabelled, what = "replicates"), class = "latents_bad_argument")
  restated <- boot
  restated$estimate[2L] <- restated$estimate[2L] + 1
  expect_error(as.data.frame(restated, what = "replicates"), class = "latents_bad_argument")
  dropped <- boot
  dropped$term <- NULL
  expect_s3_class(dropped, "latents_bootstrap")
  expect_error(as.data.frame(dropped, what = "replicates"), class = "latents_bad_argument")
  selected <- boot[, setdiff(names(boot), "term")]
  expect_error(as.data.frame(selected, what = "replicates"), class = "latents_bad_argument")
  no_estimate <- boot
  no_estimate$estimate <- NULL
  expect_error(as.data.frame(no_estimate, what = "replicates"), class = "latents_bad_argument")
  skip_if_not_installed("ggplot2")
  expect_error(plot(relabelled), class = "latents_bad_argument")
  expect_error(plot(dropped), class = "latents_bad_argument")
})

test_that("rescaled or shifted data keep every varying parameter on the plot", {
  skip_if_not_installed("ggplot2")
  original <- review_bootstrap()
  small <- review_bootstrap(scale = 1e-12)
  shifted <- review_bootstrap(shift = 1e6)
  # The rescaled means vary by about 1e-13 to 1e-12 across the replicates,
  # far below any fixed threshold in the data's units.
  ranges <- vapply(seq_len(4L), function(column) {
    diff(range(attr(small, "replicates")[, column]))
  }, numeric(1))
  expect_true(all(ranges > 0 & ranges < 1e-11))
  reference <- plot(original)$data
  lapply(list(small, shifted), function(boot) {
    drawn <- plot(boot)$data
    expect_identical(drawn$panel, reference$panel)
    expect_identical(drawn$at, reference$at)
    expect_identical(drawn$series, reference$series)
    means <- plot(boot, what = "mean")$data
    expect_identical(nrow(means), 4L)
  })
  # The estimates and intervals carry the rescaling: means by the scale,
  # variances by its square, probabilities unchanged.
  unit_scale <- unname(c(mean = 1e-12, variance = 1e-24, probability = 1)[original$parameter])
  expect_equal(small$estimate / unit_scale, original$estimate, tolerance = 1e-5)
  expect_equal(small$conf_low / unit_scale, original$conf_low, tolerance = 1e-5)
  expect_equal(small$conf_high / unit_scale, original$conf_high, tolerance = 1e-5)
  offset <- ifelse(original$parameter == "mean", 1e6, 0)
  expect_equal(shifted$estimate - offset, original$estimate, tolerance = 1e-5)
  expect_equal(shifted$conf_high - offset, original$conf_high, tolerance = 1e-5)
})

test_that("a constant parameter is still classified constant", {
  skip_if_not_installed("ggplot2")
  # The single group class's probability is one in every replicate.
  lapply(list(review_bootstrap(), review_bootstrap(scale = 1e-12)), function(boot) {
    shares <- plot(boot, what = "probability")
    expect_identical(nrow(shares$data), 2L)
    expect_match(shares$labels$caption, "group_class_1", fixed = TRUE)
  })
  set.seed(7)
  kept <- cbind(fixed = rep(0.3, 5L),
                zero = rep(0, 5L),
                rounding = 1 - c(0, 1, 2, 1, 0) * .Machine$double.eps,
                tiny = 1e-200 * (1 + stats::runif(5L)),
                far = 1e8 + stats::rnorm(5L, sd = 1e-3),
                ordinary = stats::rnorm(5L))
  expect_identical(.latents_bootstrap_spread(kept),
                   c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  # Scale equivariance: the classification does not depend on the units.
  expect_identical(.latents_bootstrap_spread(kept * 1e-30), .latents_bootstrap_spread(kept))
  expect_identical(.latents_bootstrap_spread(kept * 1e30), .latents_bootstrap_spread(kept))
  expect_error(.latents_bootstrap_spread(as.data.frame(kept)))
})

test_that("the transition and growth bootstraps carry the same identity guard", {
  set.seed(8)
  n <- 60L
  first <- sample.int(2L, n, replace = TRUE)
  second <- ifelse(stats::runif(n) < 0.8, first, 3L - first)
  long <- data.frame(id = rep(seq_len(n), each = 2L), time = rep(1:2, n),
                     y = stats::rnorm(2L * n, c(-1.5, 1.5)[c(rbind(first, second))]))
  transition <- suppressWarnings(lta(long, "y", id = "id", time = "time", n_profiles = 2,
                                     n_starts = 2, seed = 1))
  lta_boot <- suppressWarnings(parameter_inference(transition, data = long,
                                                   method = "bootstrap", iter = 8,
                                                   n_starts = 1, seed = 1))
  lta_labels <- .latents_bootstrap_columns(lta_boot)$labels
  expect_reproduces_table(lta_boot, lta_labels)
  reversed <- lta_boot[rev(seq_len(nrow(lta_boot))), , drop = FALSE]
  row.names(reversed) <- NULL
  expect_error(as.data.frame(reversed, what = "replicates"), class = "latents_bad_argument")

  set.seed(3)
  people <- 80L
  times <- 0:4
  x <- stats::rnorm(people)
  class <- 1L + (stats::runif(people) < stats::plogis(0.2 + 0.8 * x))
  growth <- data.frame(id = rep(seq_len(people), each = length(times)),
                       time = rep(times, people), x = rep(x, each = length(times)))
  growth$y <- c(2, 5)[class[growth$id]] + c(0.6, -0.4)[class[growth$id]] * growth$time +
    rep(stats::rnorm(people), each = length(times)) + stats::rnorm(nrow(growth), 0, 0.7)
  fit <- suppressWarnings(mixture_regression(y ~ time, growth, n_classes = 2, id = "id",
                                             class_level = "group", random = ~ 1,
                                             n_starts = 2, seed = 1, vcov_type = "none"))
  growth_boot <- suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 4,
                                                      n_starts = 1, seed = 5))
  expect_s3_class(growth_boot, "latents_bootstrap")
  expect_reproduces_table(growth_boot, .latents_bootstrap_columns(growth_boot)$labels)
  reordered <- growth_boot[rev(seq_len(nrow(growth_boot))), , drop = FALSE]
  row.names(reordered) <- NULL
  expect_error(as.data.frame(reordered, what = "replicates"), class = "latents_bad_argument")
})
