# A held measurement block belongs to an indicator and a category, not to a
# position in a list. These tests fix that identity: naming the same indicators
# in a different order, or declaring the same factor's levels in a different
# order, must describe the same model and report the same numbers. The old
# behaviour consumed the blocks positionally, so a permutation silently rebuilt
# a different model and returned a different likelihood.

# Every ordering of a set, so the invariance is asserted over all of them rather
# than over one hand-picked swap. Brackets belong inside a function body.
.staged_orderings <- function(items) {
  if (length(items) <= 1L) return(list(items))
  unlist(lapply(seq_along(items), function(index) {
    lapply(.staged_orderings(items[-index]),
           function(rest) c(items[[index]], rest))
  }), recursive = FALSE)
}

# Two kinds of school, differing in how often a pupil answers "yes", so the
# profiles are well separated and every stage converges without qualification.
.staged_frame <- function(seed = 127L, n_groups = 30L, size = 8L) {
  set.seed(seed)
  n <- n_groups * size
  state <- rep(c(0, 1), each = n / 2L)
  answer <- function(high) {
    factor(ifelse(stats::runif(n) < ifelse(state == 1, high, 1 - high),
                  "yes", "no"))
  }
  data.frame(school = rep(seq_len(n_groups), each = size),
             a = answer(0.9), b = answer(0.85), c = answer(0.8),
             score = stats::rnorm(n, state * 4))
}

test_that("the staged fit does not depend on the order `categorical` names", {
  frame <- .staged_frame()
  items <- c("a", "b", "c")
  indicators <- c(items, "score")
  measurement <- multilpa(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 1, categorical = items,
                          n_starts = 2, seed = 1)
  reference <- fit_staged(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 2, categorical = items,
                          measurement = measurement, n_starts = 2, seed = 1)

  invisible(lapply(.staged_orderings(items), function(ordering) {
    label <- paste(ordering, collapse = ", ")
    permuted <- fit_staged(frame, indicators, "school", n_profiles = 2,
                           n_group_classes = 2, categorical = ordering,
                           measurement = measurement, n_starts = 2, seed = 1)
    expect_equal(permuted$log_likelihood, reference$log_likelihood,
                 tolerance = 1e-8, label = label)
    expect_identical(permuted$n_parameters, reference$n_parameters)
    # Each held distribution must arrive at the indicator whose name it carries,
    # and be the one the first stage estimated for that indicator.
    expect_named(permuted$response_probabilities, ordering)
    invisible(lapply(items, function(item) {
      expect_equal(permuted$response_probabilities[[item]],
                   measurement$response_probabilities[[item]],
                   label = paste(label, item))
    }))
  }))
})

test_that("the staged fit does not depend on the order the factor levels are declared", {
  frame <- .staged_frame()
  items <- c("a", "b", "c")
  indicators <- c(items, "score")
  measurement <- multilpa(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 1, categorical = items,
                          n_starts = 2, seed = 1)
  reference <- fit_staged(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 2, categorical = items,
                          measurement = measurement, n_starts = 2, seed = 1)

  # The same responses, with one item's levels declared the other way round.
  # Nothing about the data has changed, so nothing about the fit may change.
  reversed <- frame
  reversed$a <- factor(reversed$a, levels = rev(levels(reversed$a)))
  flipped <- fit_staged(reversed, indicators, "school", n_profiles = 2,
                        n_group_classes = 2, categorical = items,
                        measurement = measurement, n_starts = 2, seed = 1)
  expect_equal(flipped$log_likelihood, reference$log_likelihood,
               tolerance = 1e-8)
  expect_identical(colnames(flipped$response_probabilities$a),
                   rev(levels(frame$a)))
  # The probability of "yes" is the probability of "yes" whichever column it is
  # printed in, so read both back by their own labels.
  invisible(lapply(levels(frame$a), function(category) {
    expect_equal(flipped$response_probabilities$a[, category],
                 measurement$response_probabilities$a[, category],
                 label = category)
  }))

  # And with every level order reversed at once.
  all_reversed <- frame
  all_reversed[items] <- lapply(items, function(item) {
    factor(frame[[item]], levels = rev(levels(frame[[item]])))
  })
  every <- fit_staged(all_reversed, indicators, "school", n_profiles = 2,
                      n_group_classes = 2, categorical = items,
                      measurement = measurement, n_starts = 2, seed = 1)
  expect_equal(every$log_likelihood, reference$log_likelihood, tolerance = 1e-8)
})

test_that("a first stage whose categories are not these categories is refused", {
  frame <- .staged_frame()
  items <- c("a", "b", "c")
  indicators <- c(items, "score")
  measurement <- multilpa(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 1, categorical = items,
                          n_starts = 2, seed = 1)
  # The same item, recorded on a different scale. There is no alignment to
  # make, so the request must be refused rather than answered with a number.
  relabelled <- frame
  relabelled$a <- factor(ifelse(frame$a == "yes", "agree", "disagree"))
  expect_error(
    fit_staged(relabelled, indicators, "school", n_profiles = 2,
               n_group_classes = 2, categorical = items,
               measurement = measurement, n_starts = 2, seed = 1),
    class = "multilpa_bad_stage")

  # An extra category the first stage never saw cannot be held either.
  extended <- frame
  extended$b <- factor(ifelse(seq_len(nrow(frame)) %% 40L == 0L, "maybe",
                              as.character(frame$b)))
  expect_error(
    fit_staged(extended, indicators, "school", n_profiles = 2,
               n_group_classes = 2, categorical = items,
               measurement = measurement, n_starts = 2, seed = 1),
    class = "multilpa_bad_stage")
})

test_that("a labelled start is aligned by label, and a label that cannot be aligned is refused", {
  frame <- .staged_frame()
  items <- c("a", "b", "c")
  indicators <- c(items, "score")
  measurement <- multilpa(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 1, categorical = items,
                          n_starts = 2, seed = 1)
  start <- starting_values(measurement, what = "measurement")

  # The labelled blocks handed over in the first stage's own order, held by a
  # second stage that names the indicators in the opposite order.
  labelled <- start
  labelled$response_probabilities <- measurement$response_probabilities
  reordered <- multilpa(frame, indicators, "school", n_profiles = 2,
                        n_group_classes = 2, categorical = rev(items),
                        n_starts = 2, seed = 1, start = labelled,
                        fixed = "measurement")
  expect_named(reordered$response_probabilities, rev(items))
  invisible(lapply(items, function(item) {
    expect_equal(reordered$response_probabilities[[item]],
                 measurement$response_probabilities[[item]], label = item)
  }))

  # A name this fit has no indicator for is refused, by class.
  wrong_item <- labelled
  names(wrong_item$response_probabilities) <- c("a", "b", "elsewhere")
  expect_error(
    multilpa(frame, indicators, "school", n_profiles = 2, n_group_classes = 2,
             categorical = items, n_starts = 1, seed = 1, start = wrong_item,
             fixed = "measurement"),
    class = "multilpa_bad_start")

  # A category label this indicator does not have is refused too.
  wrong_category <- labelled
  colnames(wrong_category$response_probabilities$a) <- c("no", "perhaps")
  expect_error(
    multilpa(frame, indicators, "school", n_profiles = 2, n_group_classes = 2,
             categorical = items, n_starts = 1, seed = 1,
             start = wrong_category, fixed = "measurement"),
    class = "multilpa_bad_start")
})

test_that("an unlabelled start keeps the documented positional contract", {
  frame <- .staged_frame()
  items <- c("a", "b", "c")
  indicators <- c(items, "score")
  measurement <- multilpa(frame, indicators, "school", n_profiles = 2,
                          n_group_classes = 1, categorical = items,
                          n_starts = 2, seed = 1)
  # `starting_values()` now labels the response blocks, so an unlabelled start
  # has to be built deliberately. It is still a supported input: a start stored
  # by an earlier version, or assembled by hand, carries no labels and must keep
  # being read positionally rather than rejected.
  start <- starting_values(measurement, what = "measurement")
  expect_identical(names(start$response_probabilities), items)
  start$response_probabilities <- unname(lapply(start$response_probabilities,
                                                unname))
  expect_null(names(start$response_probabilities))
  positional <- multilpa(frame, indicators, "school", n_profiles = 2,
                         n_group_classes = 2, categorical = items,
                         n_starts = 2, seed = 1, start = start,
                         fixed = "measurement")
  invisible(lapply(items, function(item) {
    expect_equal(positional$response_probabilities[[item]],
                 measurement$response_probabilities[[item]], label = item)
  }))
})
