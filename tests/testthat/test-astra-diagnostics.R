# Astra A09/A10/A16/A17: classification eligibility, nesting and weights must
# agree in the fitted-model tables, diagnostics and enumeration plot data.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

astra_diagnostic_fits <- local({
  saved <- NULL
  function() {
    if (!is.null(saved)) return(saved)
    withr::local_seed(912)
    groups <- 60L
    n <- groups * 6L
    type <- sample(1:2, groups, TRUE)
    state <- vapply(rep(type, each = 6L), function(h) {
      sample(1:2, 1L, prob = if (h == 1L) c(.8, .2) else c(.2, .8))
    }, integer(1))
    data <- data.frame(id = rep(seq_len(groups), each = 6L), time = rep(1:6, groups),
      x = stats::rnorm(n), y = stats::rnorm(n, c(-1.2, 1.2)[state]),
      z = stats::rnorm(n, c(-1, 1)[state]),
      w = rep(c(stats::runif(groups - 5L, .1, 4), rep(0, 5L)), each = 6L))
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    grid <- quiet(enumerate_regressions(y ~ x, data, n_classes = 2,
      n_group_classes = 2, id = "id", weights = "w", n_starts = 2, seed = 3))
    single_grid <- quiet(enumerate_regressions(y ~ x, data, n_classes = 2,
      weights = "w", n_starts = 2, seed = 3))
    group_grid <- quiet(enumerate_regressions(y ~ x, data, n_classes = 2,
      id = "id", class_level = "group", weights = "w", n_starts = 2, seed = 3))
    trailing <- data[data$time == 6L & data$id <= 20L, ]
    trailing$time <- 7L
    trailing[c("y", "z")] <- NA_real_
    transition_data <- rbind(data, trailing)
    transition_fit <- function(drop, transitions) quiet(lta(transition_data,
      c("y", "z"), "id", 2, "time", transitions = transitions,
      missing = "fiml", drop_empty = drop, weights = "w", n_starts = 1, seed = 2))
    saved <<- list(data = data, transition_data = transition_data,
      grid = grid, single_grid = single_grid, group_grid = group_grid,
      regression = get_results(grid, "model", n_classes = 2, n_group_classes = 2),
      single = get_results(single_grid, "model", n_classes = 2),
      group = get_results(group_grid, "model", n_classes = 2),
      one = quiet(mixture_regression(y ~ x, data, n_classes = 1, n_starts = 1,
        seed = 3, vcov_type = "none")),
      lta_drop = transition_fit(TRUE, "occasion"),
      lta_keep = transition_fit(FALSE, "occasion"),
      homogeneous = transition_fit(TRUE, "homogeneous"),
      cross = quiet(multilpa(data, c("y", "z"), "id", 2, 2,
        family = "restricted_cross_level", weights = "w", n_starts = 1, seed = 2)),
      additive = quiet(multilpa(data, c("y", "z"), "id", family = "additive",
        n_group_classes = 2, weights = "w", n_starts = 1, seed = 2)))
    saved
  }
})

astra_hand_quality <- function(p, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(p))
  kept <- w > 0
  p <- p[kept, , drop = FALSE]
  w <- w[kept]
  assigned <- max.col(p, ties.method = "first")
  average <- vapply(seq_len(ncol(p)), function(k) {
    mine <- assigned == k
    if (any(mine)) sum(w[mine] * p[mine, k]) / sum(w[mine]) else NA_real_
  }, numeric(1))
  list(counts = colSums(p * w), n = nrow(p), modal = tabulate(assigned, ncol(p)),
       average = average, entropy = if (ncol(p) < 2L) NA_real_ else
         1 + sum(w * rowSums(ifelse(p > 0, p * log(p), 0))) / (sum(w) * log(ncol(p))))
}

test_that("cross-level diagnostics exclude zero weights at both mapped levels", {
  fits <- astra_diagnostic_fits()
  fit <- fits$cross
  dg <- diagnostics(fit, plots = FALSE)
  expected <- astra_hand_quality(fit$group_posteriors, fit$sampling_weights)
  direct <- .latents_level_classification(fit, "groups")
  expect_equal(direct$estimated_n, unname(expected$counts))
  expect_identical(direct$n_assigned, expected$modal)
  public <- get_results(fit, "group_classes")
  expect_identical(public$n_assigned, expected$modal)
  expect_equal(sum(public$n_assigned), nobs(fit))
  table <- get_results(dg, "classification")
  groups <- table[table$level == "groups", ]
  expect_equal(groups$estimated_n, unname(expected$counts))
  expect_identical(groups$n_modal, expected$modal)
  expect_equal(sum(groups$estimated_n), nobs(fit))
  row_weights <- fit$sampling_weights[fit$sufficient_statistics$index]
  expected_rows <- astra_hand_quality(fit$subject_posteriors, row_weights)
  rows <- table[table$level == "individuals", ]
  expect_equal(rows$estimated_n, unname(expected_rows$counts))
  expect_identical(rows$n_modal, expected_rows$modal)
  expect_equal(rows$average_posterior, expected_rows$average)
})

test_that("both LTA engines obey empty-row policy and zero-weight eligibility", {
  fits <- astra_diagnostic_fits()
  invisible(lapply(c("lta_drop", "lta_keep", "homogeneous"), function(name) {
    fit <- fits[[name]]
    observed <- rowSums(!is.na(fits$transition_data[c("y", "z")])) > 0
    keep <- (isFALSE(fit$drop_empty) | observed) & fits$transition_data$w > 0
    w <- fit$sampling_weights[fit$group_index][keep]
    expected <- astra_hand_quality(fit$subject_posteriors[keep, , drop = FALSE], w)
    dg <- diagnostics(fit, plots = FALSE)
    checks <- get_results(dg, "checks")
    expect_equal(checks$value[checks$check == "relative entropy" &
                                checks$level == "individuals"], expected$entropy, label = name)
    entropy <- get_results(dg, "entropy")
    rows <- entropy$level == "individuals"
    expect_equal(entropy$n_units[rows], sum(keep), label = name)
    expect_equal(entropy$relative_entropy[rows], expected$entropy, label = name)
    classification <- get_results(dg, "classification")
    classification <- classification[classification$level == "individuals", ]
    expect_equal(classification$estimated_n, unname(expected$counts), label = name)
    expect_identical(classification$n_modal, expected$modal, label = name)
  }))
})

test_that("regression diagnostics use the real classification levels and counts", {
  fits <- astra_diagnostic_fits()
  two <- get_results(diagnostics(fits$regression, plots = FALSE), "checks")
  expect_true(all(c("individuals", "groups") %in% two$level))
  groups <- get_results(diagnostics(fits$regression, plots = FALSE), "classification")
  groups <- groups[groups$level == "groups", ]
  reference <- astra_hand_quality(fits$regression$expectation$rho,
                                  fits$regression$spec$sampling_weights)
  expect_equal(groups$estimated_n, unname(reference$counts))
  expect_identical(groups$n_modal, reference$modal)
  grouped <- get_results(diagnostics(fits$group, plots = FALSE), "checks")
  expect_true("groups" %in% grouped$level)
  expect_false("individuals" %in% grouped$level)
  one <- get_results(diagnostics(fits$one, plots = FALSE), "checks")
  expect_true(all(one$status[one$domain == "Classification"] == "n/a"))
  expect_true(all(grepl("one class", one$detail[one$domain == "Classification"])))
})

test_that("every newly supported family has entropy as the primary diagnostic table", {
  fits <- astra_diagnostic_fits()
  invisible(lapply(c("lta_drop", "cross", "additive", "regression", "single", "group", "one"),
    function(name) {
      dg <- diagnostics(fits[[name]], plots = FALSE)
      entropy <- get_results(dg, "entropy")
      expect_equal(as.data.frame(dg), entropy, label = name)
      expect_true(all(c("entropy", "classification", "average_posteriors") %in%
                        names(get_results(dg, "all"))), label = name)
      expect_named(entropy, c("level", "n_classes", "n_units", "entropy_sum", "relative_entropy"))
    }))
})

test_that("weighted regression quality views agree with their candidate summaries", {
  fits <- astra_diagnostic_fits()
  invisible(lapply(c("grid", "single_grid", "group_grid"), function(name) {
    grid <- fits[[name]]
    fit <- grid$fits[[1L]]
    table <- .latents_grid_quality_table(grid)
    grouped <- identical(fit$spec$nesting, "group")
    p <- if (grouped) fit$expectation$group_tau else fit$expectation$tau
    w <- if (grouped) fit$spec$sampling_weights else fit$spec$row_weights
    expected <- astra_hand_quality(p, w)
    expect_equal(table$avepp_min, min(expected$average, na.rm = TRUE), label = name)
    classes <- get_results(fit, "classification")
    expect_equal(table$avepp_min,
      min(classes$mean_posterior[classes$assigned == classes$class], na.rm = TRUE), label = name)
    if (identical(fit$spec$nesting, "two-level")) {
      groups <- astra_hand_quality(fit$expectation$rho, fit$spec$sampling_weights)
      expect_equal(table$group_entropy, groups$entropy)
      expect_equal(table$group_avepp_min, min(groups$average, na.rm = TRUE))
    } else expect_true(is.na(table$group_entropy))
  }))
})

test_that("hierarchical growth diagnostics retain person and cluster classification", {
  withr::local_seed(3)
  schools <- do.call(rbind, lapply(seq_len(14), function(s) {
    improving <- if (s %% 2L == 1L) .85 else .2
    do.call(rbind, lapply(seq_len(8), function(i) {
      k <- if (stats::runif(1) < improving) 1L else 2L
      data.frame(school = s, student = (s - 1) * 8 + i, wave = 0:3,
        score = (if (k == 1L) 50 + 3 * 0:3 else rep(55, 4)) +
          stats::rnorm(1, 0, 2) + stats::rnorm(4, 0, 1.5))
    }))
  }))
  fit <- suppressWarnings(mixture_regression(score ~ wave, schools, n_classes = 2,
    id = "student", class_level = "group", random = ~ 1, n_group_classes = 2,
    cluster = "school", n_starts = 1, seed = 1, vcov_type = "none"))
  dg <- diagnostics(fit, plots = FALSE)
  entropy <- get_results(dg, "entropy")
  expect_identical(entropy$level, c("individuals", "groups"))
  expect_equal(entropy$n_units, c(112, 14))
  expected <- astra_hand_quality(fit$expectation$structure$rho)
  expect_equal(entropy$relative_entropy[entropy$level == "groups"], expected$entropy)
  expect_equal(as.data.frame(dg), entropy)
})


test_that("weighted growth diagnostics and grid views use eligible persons", {
  data <- astra_diagnostic_fits()$data
  grid <- suppressWarnings(enumerate_regressions(y ~ time, data, n_classes = 2,
    id = "id", class_level = "group", random = ~ 1, weights = "w",
    n_starts = 1, seed = 3))
  fit <- get_results(grid, "model", n_classes = 2)
  expected <- astra_hand_quality(fit$expectation$posterior, fit$spec$sampling_weights)
  dg <- diagnostics(fit, plots = FALSE)
  entropy <- get_results(dg, "entropy")
  expect_equal(entropy$n_units, expected$n)
  expect_equal(entropy$relative_entropy, expected$entropy)
  expect_equal(as.data.frame(dg), entropy)
  quality <- .latents_grid_quality_table(grid)
  expect_equal(quality$avepp_min, min(expected$average, na.rm = TRUE))
  classes <- get_results(dg, "classification")
  expect_equal(classes$estimated_n, unname(expected$counts))
  expect_identical(classes$n_modal, expected$modal)
})

test_that("older weighted cross-level objects explicitly mark an unavailable row map", {
  fit <- astra_diagnostic_fits()$cross
  fit$sufficient_statistics$index <- NULL
  fit$group_index <- NULL
  dg <- diagnostics(fit, plots = FALSE)
  entropy <- get_results(dg, "entropy")
  expect_true(all(is.na(entropy$relative_entropy[entropy$level == "individuals"])))
  expect_equal(as.data.frame(dg), entropy)
  checks <- get_results(dg, "checks")
  rows <- checks$level == "individuals"
  expect_true(all(checks$status[rows] == "n/a"))
  expect_match(checks$detail[rows], "no map")
  # A saved one-group-class object can have no other classifying level.
  fit$group_posteriors <- fit$group_posteriors[, 1L, drop = FALSE]
  dg <- diagnostics(fit, plots = FALSE)
  expect_equal(nrow(get_results(dg, "classification")), 0L)
  expect_equal(nrow(get_results(dg, "average_posteriors")), 0L)
  expect_equal(as.data.frame(dg), get_results(dg, "entropy"))
})
