# Review 2026-11-01, F8: align_labels() matches category probabilities by
# position, so two fits must code each ordinal and categorical indicator with
# the same categories in the same order. Equal category counts are not enough.

skip_on_cran()

quiet_fit <- function(expr) {
  withCallingHandlers(expr,
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
}

# Two separated profiles; the ordinal `u` and the categorical `k` both carry
# the profile, `a` is continuous and `n_events` a count.
make_data <- function(seed, n = 400) {
  set.seed(seed)
  state <- sample.int(2L, n, replace = TRUE)
  pick <- function(choices, probabilities) {
    vapply(state, function(s) sample(choices, 1L, prob = probabilities[[s]]), character(1))
  }
  data.frame(
    a = stats::rnorm(n, c(-2, 2)[state]),
    u = ordered(pick(c("low", "medium", "high"), list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))),
                levels = c("low", "medium", "high")),
    v = ordered(pick(c("never", "sometimes", "often"), list(c(0.5, 0.4, 0.1), c(0.1, 0.4, 0.5))),
                levels = c("never", "sometimes", "often")),
    k = factor(pick(c("x", "y", "z"), list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))),
               levels = c("x", "y", "z")),
    n_events = stats::rpois(n, c(1, 4)[state])
  )
}

expect_refused_for <- function(fit, reference, indicator) {
  expect_error(align_labels(fit, reference), class = "latents_bad_argument")
  condition <- tryCatch(align_labels(fit, reference), latents_bad_argument = identity)
  expect_true(grepl(sprintf("`%s`", indicator), conditionMessage(condition), fixed = TRUE),
              label = conditionMessage(condition))
}

recode_ordinal <- function(data, name, levels) {
  data[[name]] <- ordered(levels[as.integer(data[[name]])], levels = levels)
  data
}

base_data <- make_data(12)
mixed <- function(data) {
  quiet_fit(lpa(data, c("a", "u", "k", "n_events"), 2, ordinal = "u", categorical = "k",
                count = "n_events", n_starts = 1, seed = 1))
}
all_ordinal <- function(data) {
  quiet_fit(lpa(data, c("u", "v"), 2, ordinal = c("u", "v"), n_starts = 1, seed = 1))
}
reference_mixed <- mixed(base_data)
reference_ordinal <- all_ordinal(base_data)

test_that("the reviewed pair with relabelled ordinal categories is refused", {
  data <- base_data[c("a", "u")]
  recoded <- recode_ordinal(data, "u", c("red", "green", "blue"))
  first <- quiet_fit(lpa(data, c("a", "u"), 2, ordinal = "u", n_starts = 1, seed = 1))
  second <- quiet_fit(lpa(recoded, c("a", "u"), 2, ordinal = "u", n_starts = 1, seed = 1))
  expect_refused_for(second, first, "u")
})

test_that("ordinal codings that differ are refused in mixed and all-ordinal models", {
  # Different labels with the same count, reversed order, a nonmonotonic
  # order of the same labels, and an extra category.
  codings <- list(
    labels = c("red", "green", "blue"),
    reversed = c("high", "medium", "low"),
    shuffled = c("medium", "low", "high")
  )
  invisible(lapply(codings, function(levels) {
    expect_refused_for(mixed(recode_ordinal(base_data, "u", levels)), reference_mixed, "u")
    expect_refused_for(all_ordinal(recode_ordinal(base_data, "u", levels)),
                       reference_ordinal, "u")
  }))
  # A fourth category: rows of the top category split into two.
  extended <- base_data
  extended$u <- ordered(ifelse(extended$u == "high" & seq_len(nrow(extended)) %% 2L == 0L,
                               "very high", as.character(extended$u)),
                        levels = c("low", "medium", "high", "very high"))
  expect_refused_for(mixed(extended), reference_mixed, "u")
  expect_refused_for(all_ordinal(extended), reference_ordinal, "u")
})

test_that("every indicator coded differently is named", {
  both <- recode_ordinal(recode_ordinal(base_data, "u", c("high", "medium", "low")),
                         "v", c("often", "sometimes", "never"))
  condition <- tryCatch(align_labels(all_ordinal(both), reference_ordinal),
                        latents_bad_argument = identity)
  expect_s3_class(condition, "latents_bad_argument")
  expect_true(grepl("`u`", conditionMessage(condition), fixed = TRUE))
  expect_true(grepl("`v`", conditionMessage(condition), fixed = TRUE))
})

test_that("categorical indicators with different category sets are refused", {
  relabelled <- base_data
  levels(relabelled$k) <- c("x", "y", "w")
  expect_refused_for(mixed(relabelled), reference_mixed, "k")
  extra <- base_data
  extra$k <- factor(ifelse(extra$k == "z" & seq_len(nrow(extra)) %% 2L == 0L, "q",
                           as.character(extra$k)), levels = c("x", "y", "z", "q"))
  expect_refused_for(mixed(extra), reference_mixed, "k")
})

test_that("fits coded the same way align exactly as before", {
  # Another sample with the same coding is accepted and aligned by the
  # unchanged alignment service.
  other <- make_data(34)
  invisible(lapply(list(list(mixed, reference_mixed), list(all_ordinal, reference_ordinal)),
                   function(pair) {
    fit <- pair[[1L]](other)
    expect_identical(align_labels(fit, pair[[2L]]), .multilpa_align_labels(fit, pair[[2L]]))
    # A fit with its labels swapped is put back on the reference's labels.
    swapped <- .multilpa_permute_profiles(pair[[2L]], 2:1)
    expect_false(identical(swapped$subject_profiles, pair[[2L]]$subject_profiles))
    aligned <- align_labels(swapped, pair[[2L]])
    expect_identical(aligned$subject_profiles, pair[[2L]]$subject_profiles)
    expect_equal(get_results(aligned, "profiles"), get_results(pair[[2L]], "profiles"),
                 tolerance = 1e-12)
  }))
  # A fit aligned to itself is unchanged.
  expect_identical(align_labels(reference_mixed, reference_mixed), reference_mixed)
})
