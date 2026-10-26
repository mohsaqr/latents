# Adversarial latent class data against poLCA (Linzer & Lewis 2011).
#
# The poLCA suite compares the package on poLCA's five bundled datasets, which
# are well behaved: large classes, interior probabilities, comfortable degrees
# of freedom. This suite simulates the data on which latent class software
# goes wrong: a category seen in three rows, classes that never use some
# categories (probabilities on the bound in the first, a middle and the last
# position), a 5% class, weakly separated classes, more items than rows can
# support, items with up to eight categories, a model three degrees of freedom
# short of saturation, the boundary data with every item's categories in
# reverse order, a class with deterministic answers, and four classes whose
# likelihood has many local maxima.
#
# Every case makes two separate claims, so that a disagreement says whether
# the package computes wrongly or searches badly:
#   (A) computation - the package is started from poLCA's estimates, refined to
#       poLCA's fixed point, and iterates normally from there. At the same mode
#       every output is compared through `compare_fit()`, against tables built
#       from poLCA's estimates and its own E step, with the same independent
#       likelihood, observed-information and OPG references as the poLCA suite
#       (whose helpers are sourced below, not copied).
#   (B) search - the package's own random starts against poLCA's best over
#       many starts. A shortfall of the package's maximum is a finding about
#       its optimiser; a maximum above poLCA's is recorded as a row that is not
#       an output of the fit, and is a finding about poLCA's search.
# The reversed case is also compared with the boundary case through every
# output table: category order is a labelling, so nothing may depend on it.

sys.source(file.path("equivalence", "harness", "suites", "polca.R"), envir = environment())
sys.source(file.path("equivalence", "harness", "suites", "invariance.R"), envir = environment())

suite_hard_polca <- function() {
  require_suite_packages(
    "poLCA", "numDeriv",
    reason = "the single-level categorical reference and the reference observed information")
  fixtures <- .hard_fixtures()
  results <- lapply(fixtures, .hard_case)
  names(results) <- vapply(fixtures, `[[`, character(1), "name")
  reversal <- .hard_reversal(results$boundary, results$reversed)
  stack_cases(c(lapply(results, `[[`, "case"), list(reversal)))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' Draw latent class data with exactly fixed class sizes
#'
#' Class sizes are fixed rather than drawn so that "a 5% class" is 25 rows
#' under every seed, not whatever a multinomial draw gives.
#'
#' @param sizes Rows per class.
#' @param probabilities Named list, one classes-by-categories matrix per item,
#'   columns named by category label in level order.
#' @param seed The seed.
#' @return A `data.frame` of factors, one per item, rows shuffled.
.hard_simulate <- function(sizes, probabilities, seed) {
  set.seed(seed)
  class <- sample(rep(seq_along(sizes), sizes))
  frame <- as.data.frame(lapply(probabilities, function(block) {
    stopifnot("each class's probabilities must sum to one" =
                isTRUE(all.equal(unname(rowSums(block)), rep(1, nrow(block)))))
    codes <- vapply(class, function(k) sample.int(ncol(block), 1L, prob = block[k, ]),
                    integer(1))
    factor(colnames(block)[codes], levels = colnames(block))
  }), stringsAsFactors = FALSE)
  names(frame) <- names(probabilities)
  frame
}

#' A classes-by-categories matrix from its rows
#' @param ... One probability vector per class.
#' @param labels Category labels.
#' @return The matrix.
.hard_block <- function(..., labels = letters[seq_along(list(...)[[1L]])]) {
  block <- do.call(rbind, list(...))
  colnames(block) <- labels
  block / rowSums(block)
}

#' The adversarial fixtures
#'
#' Each fixture records the data, the class count and why it is hard; the
#' class probabilities that generated it are in the code beside it.
#'
#' @return A list of fixtures.
.hard_fixtures <- function() {
  binary <- function(...) lapply(list(...), function(p) .hard_block(c(p[[1L]], 1 - p[[1L]]),
                                                                   c(p[[2L]], 1 - p[[2L]]),
                                                                   if (length(p) > 2L)
                                                                     c(p[[3L]], 1 - p[[3L]])))

  # 1. One category of a four-category item seen in exactly three rows. The
  # draw never gives it, and three rows are then recoded to it, so its
  # estimated probability is about 3 / (class size): interior but carried by
  # three observations.
  sparse_probabilities <- c(
    list(S = .hard_block(c(0.6, 0.3, 0.1, 0), c(0.1, 0.3, 0.6, 0))),
    stats::setNames(binary(c(0.85, 0.2), c(0.8, 0.15), c(0.75, 0.25), c(0.9, 0.3)),
                    c("B1", "B2", "B3", "B4")))
  sparse <- .hard_simulate(c(220L, 180L), sparse_probabilities, seed = 20261026L)
  sparse$S[c(7L, 111L, 305L)] <- "d"

  # 2. Classes that never use some categories. Class 1 never answers `a` on
  # I1 (first position), `b` on I2 (middle), `d` on I3 (last); class 2 never
  # answers `d` on I1 and `a` on I4; class 3 never answers `c` or `d` on I4.
  # The maximum puts six of these seven cells on the bound, in first, middle
  # and last positions; class 3's `d` on I4 stays interior, because rows
  # answering it keep some posterior weight in class 3 (measured 2026-10-26).
  boundary_probabilities <- list(
    I1 = .hard_block(c(0, 0.2, 0.3, 0.5), c(0.6, 0.3, 0.1, 0), c(0.25, 0.25, 0.25, 0.25)),
    I2 = .hard_block(c(0.7, 0, 0.2, 0.1), c(0.1, 0.2, 0.3, 0.4), c(0.1, 0.6, 0.2, 0.1)),
    I3 = .hard_block(c(0.5, 0.4, 0.1, 0), c(0.05, 0.15, 0.3, 0.5), c(0.3, 0.1, 0.3, 0.3)),
    I4 = .hard_block(c(0.1, 0.1, 0.2, 0.6), c(0, 0.2, 0.4, 0.4), c(0.5, 0.5, 0, 0)),
    I5 = .hard_block(c(0.9, 0.1), c(0.2, 0.8), c(0.5, 0.5)))
  boundary <- .hard_simulate(c(250L, 200L, 150L), boundary_probabilities, seed = 20261027L)

  # 3. A class of 25 of 500 rows (5%). The two large classes answer the
  # first four items and the last four in opposite directions; the small one
  # says yes to everything, a pattern no mixture of the large classes makes,
  # so the maximum keeps it as a class of its own rather than splitting a
  # large one (fitted share 0.063, measured 2026-10-26).
  small_probabilities <- stats::setNames(
    binary(c(0.95, 0.85, 0.15), c(0.9, 0.9, 0.1), c(0.95, 0.8, 0.2), c(0.9, 0.85, 0.15),
           c(0.95, 0.15, 0.85), c(0.9, 0.1, 0.9), c(0.95, 0.2, 0.8), c(0.9, 0.15, 0.85)),
    sprintf("Q%d", seq_len(8L)))
  small <- .hard_simulate(c(25L, 240L, 235L), small_probabilities, seed = 20261028L)

  # 4. Weakly separated classes: every item differs by 0.3 between the
  # classes. The fitted relative entropy is 0.409 (measured 2026-10-26).
  weak_probabilities <- stats::setNames(
    binary(c(0.65, 0.35), c(0.35, 0.65), c(0.7, 0.4), c(0.6, 0.3), c(0.3, 0.6),
           c(0.65, 0.35)),
    sprintf("W%d", seq_len(6L)))
  weak <- .hard_simulate(c(330L, 270L), weak_probabilities, seed = 20261029L)

  # 5. Ten binary items, 80 rows, three classes: 32 parameters from 80 rows.
  many_probabilities <- stats::setNames(
    binary(c(0.9, 0.1, 0.5), c(0.85, 0.2, 0.8), c(0.8, 0.1, 0.2), c(0.9, 0.3, 0.9),
           c(0.2, 0.15, 0.85), c(0.85, 0.15, 0.15), c(0.75, 0.8, 0.1), c(0.9, 0.2, 0.6),
           c(0.1, 0.25, 0.8), c(0.8, 0.85, 0.15)),
    sprintf("M%02d", seq_len(10L)))
  many <- .hard_simulate(c(30L, 28L, 22L), many_probabilities, seed = 20261030L)

  # 6. Items with six to eight categories, three classes, 600 rows: 71
  # parameters, many cells with expected counts below five.
  wide_row <- function(width, peak) {
    p <- exp(-0.6 * abs(seq_len(width) - peak))
    p / sum(p)
  }
  wide_probabilities <- list(
    V6 = .hard_block(wide_row(6, 1), wide_row(6, 3.5), wide_row(6, 6)),
    V7 = .hard_block(wide_row(7, 7), wide_row(7, 2), wide_row(7, 4)),
    V8 = .hard_block(wide_row(8, 2), wide_row(8, 8), wide_row(8, 5)),
    X6 = .hard_block(wide_row(6, 4), wide_row(6, 1), wide_row(6, 6)))
  wide <- .hard_simulate(c(240L, 200L, 160L), wide_probabilities, seed = 20261031L)

  # 7. Items with two, three and four categories and three classes: 23
  # degrees of freedom in the response table against 20 parameters, three to
  # spare. Kruskal's condition (k-ranks 2 + 3 + 3 >= 2 x 3 + 2) holds, so the
  # model is generically identified, with nothing to spare.
  saturated_probabilities <- list(
    T2 = .hard_block(c(0.9, 0.1), c(0.15, 0.85), c(0.6, 0.4)),
    T3 = .hard_block(c(0.8, 0.15, 0.05), c(0.1, 0.2, 0.7), c(0.1, 0.8, 0.1)),
    T4 = .hard_block(c(0.7, 0.2, 0.05, 0.05), c(0.05, 0.1, 0.25, 0.6),
                     c(0.1, 0.6, 0.2, 0.1)))
  saturated <- .hard_simulate(c(400L, 350L, 250L), saturated_probabilities, seed = 20261101L)

  # 8. The boundary data with every item's categories in reverse order: the
  # same answers, the same labels, the opposite level order, so the floored
  # cells move from first to last position and back.
  reversed <- boundary
  reversed[] <- lapply(boundary, function(column) factor(column, levels = rev(levels(column))))

  # 9. A deterministic class: class 1 always answers `a` on D1 to D3 and
  # always `b` on the three-category E1, so its row of each of those items is
  # one cell at one minus the floor and the rest on the floor, the case where
  # a category is pinned at the top by its held siblings.
  deterministic_probabilities <- list(
    D1 = .hard_block(c(1, 0), c(0.2, 0.8), c(0.3, 0.7)),
    D2 = .hard_block(c(1, 0), c(0.25, 0.75), c(0.2, 0.8)),
    D3 = .hard_block(c(1, 0), c(0.3, 0.7), c(0.25, 0.75)),
    E1 = .hard_block(c(0, 1, 0), c(0.6, 0.1, 0.3), c(0.3, 0.1, 0.6)),
    E2 = .hard_block(c(0.4, 0.3, 0.3), c(0.7, 0.2, 0.1), c(0.1, 0.2, 0.7)),
    D4 = .hard_block(c(0.5, 0.5), c(0.85, 0.15), c(0.15, 0.85)))
  deterministic <- .hard_simulate(c(120L, 180L, 150L), deterministic_probabilities,
                                  seed = 20261102L)

  # 10. Four moderately separated classes on eight binary items and 250 rows,
  # a likelihood with several local maxima for either program's search: the
  # package's 50 starts reached 12 distinct maxima, 29 of them the best
  # (measured 2026-10-26).
  multimodal_probabilities <- stats::setNames(
    lapply(seq_len(8L), function(item) {
      high <- c(0.72, 0.7, 0.68, 0.72)
      pattern <- list(c(1, 1, 0, 0), c(1, 0, 1, 0), c(0, 1, 1, 0), c(1, 0, 0, 1),
                      c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 1, 1, 0), c(0, 0, 0, 1))[[item]]
      p <- ifelse(pattern == 1, high, 1 - high)
      .hard_block(c(p[[1L]], 1 - p[[1L]]), c(p[[2L]], 1 - p[[2L]]), c(p[[3L]], 1 - p[[3L]]),
                  c(p[[4L]], 1 - p[[4L]]))
    }),
    sprintf("H%d", seq_len(8L)))
  multimodal <- .hard_simulate(c(75L, 65L, 60L, 50L), multimodal_probabilities,
                               seed = 20261103L)

  list(
    list(name = "sparse", label = "a category seen in three rows", data = sparse, k = 2L),
    list(name = "boundary", label = "classes that never use some categories",
         data = boundary, k = 3L),
    list(name = "small", label = "a class of 5% of 500 rows", data = small, k = 3L),
    list(name = "weak", label = "weakly separated classes", data = weak, k = 2L),
    list(name = "many", label = "ten items, 80 rows, three classes", data = many, k = 3L),
    list(name = "wide", label = "items with six to eight categories", data = wide, k = 3L),
    list(name = "saturated", label = "three degrees of freedom short of saturation",
         data = saturated, k = 3L),
    list(name = "reversed", label = "the boundary data, categories reversed",
         data = reversed, k = 3L),
    list(name = "deterministic", label = "a class with deterministic answers",
         data = deterministic, k = 3L),
    list(name = "multimodal", label = "four classes, eight items, 250 rows",
         data = multimodal, k = 4L))
}

# ---------------------------------------------------------------------------
# One case: the computation claim and the search claim
# ---------------------------------------------------------------------------

#' poLCA's best solution over many random starts
#' @param observed The data.
#' @param k Class count.
#' @param nrep Random starts.
#' @param seed The seed.
#' @return A `poLCA` fit.
.hard_polca_search <- function(observed, k, nrep, seed) {
  formula <- stats::as.formula(paste("cbind(", paste(names(observed), collapse = ","), ") ~ 1"))
  set.seed(seed)
  # poLCA prints its iteration alerts with cat(); the refinement below takes
  # the winner to its fixed point whatever its own stopping rule said.
  utils::capture.output(searched <- poLCA::poLCA(
    formula, data = observed, nclass = k, nrep = nrep, maxiter = 10000L, tol = 1e-10,
    verbose = FALSE, na.rm = TRUE, calc.se = FALSE))
  searched
}

#' Run both claims on one fixture
#' @param fixture A fixture from [.hard_fixtures()].
#' @param nrep poLCA's random starts.
#' @param n_starts The package's random starts.
#' @return A list with `case` (a `harness_case`), `fit` (the computation
#'   fit), `searched` (the package's own search) and `frame` (the data fitted).
.hard_case <- function(fixture, nrep = 60L, n_starts = 50L) {
  observed <- fixture$data
  items <- names(observed)
  k <- fixture$k
  frame <- observed
  frame$unit <- factor(seq_len(nrow(frame)))
  label <- sprintf("hard: %s (%s), %d classes", fixture$name, fixture$label, k)

  searched_polca <- .hard_polca_search(observed, k, nrep, seed = 20261026L)
  reference <- .polca_refine(searched_polca)
  floor <- 1e-10
  bounded <- lapply(reference$probs, .polca_bound, floor = floor)
  point <- .polca_evaluate(bounded, reference$P, reference$codes)
  stopifnot("bounding must not move the maximum by a comparable amount" =
              abs(point$llik - reference$llik) <= 1e-7)
  levels_of <- lapply(observed, levels)
  labelled <- .polca_label(bounded, reference$codes, observed, levels_of)
  start <- starting_values(list(profile_probabilities = matrix(reference$P, 1L, k),
                                group_probabilities = 1,
                                response_probabilities = labelled))

  # (A) Started at poLCA's mode, iterating normally: at a fixed point the
  # package must stay where it starts, so every output describes poLCA's
  # point. The tolerance is the one the poLCA suite justifies (stop when the
  # likelihood stops rising).
  fit <- multilpa(frame, vars = items, id = "unit", n_profiles = k, n_group_classes = 1L,
                  categorical = items, n_starts = 1L, tol = 1e-16, max_iter = 20000L,
                  start = start)
  stopifnot("the package must converge from poLCA's estimates" = isTRUE(fit$converged))
  # The same estimates scored without a step, which tests the likelihood alone.
  at_reference <- multilpa(frame, vars = items, id = "unit", n_profiles = k,
                           n_group_classes = 1L, categorical = items, n_starts = 1L,
                           max_iter = 0L, start = start)
  textbook <- .textbook_lca_log_likelihood(reference$codes, reference$P, bounded)
  n_categories <- vapply(fit$response_probabilities, ncol, integer(1))
  residual_df <- min(fit$n_informative, prod(n_categories) - 1L) - fit$n_parameters
  # The reference observed information needs an identified maximum (the
  # poLCA suite's rule); every fixture is built to have one.
  stopifnot("the fixture must leave at least three residual degrees of freedom" =
              residual_df >= 3L)
  errors <- .polca_standard_errors(reference$codes, reference$P, bounded)

  # (B) The package's own search, against poLCA's best.
  searched <- multilpa(frame, vars = items, id = "unit", n_profiles = k, n_group_classes = 1L,
                       categorical = items, n_starts = n_starts, seed = 20261026L,
                       tol = 1e-12, max_iter = 20000L)
  package_ll <- as.numeric(logLik(searched))
  search <- rbind(
    compare_values(sprintf("%s: shortfall of the package's own search against poLCA", label),
                   0, max(point$llik - package_ll, 0), tolerance = 1e-6),
    # A maximum above poLCA's best is poLCA's search falling short; it is
    # recorded, as a non-output row, so that the reference itself is seen to
    # be a local maximum, and it says nothing against the package.
    if (package_ll > point$llik + 1e-6) {
      compare_values(sprintf(paste("%s: poLCA's best over %d starts against the package's",
                                   "higher maximum (poLCA's search, not a package defect)"),
                             label, nrep), package_ll, point$llik, tolerance = 1e-6)
    })

  extra <- rbind(
    compare_values(sprintf("%s: log-likelihood after iterating from poLCA's estimates", label),
                   point$llik, as.numeric(logLik(fit)), tolerance = 1e-6),
    compare_values(sprintf("%s: log-likelihood at poLCA's estimates", label),
                   point$llik, as.numeric(logLik(at_reference)), tolerance = 1e-6),
    compare_values(sprintf("%s: textbook log-likelihood at poLCA's estimates", label),
                   point$llik, textbook, tolerance = 1e-6),
    compare_values(sprintf("%s: poLCA's BIC", label), searched_polca$bic, fit$bic_individual,
                   tolerance = 1e-5),
    compare_values(sprintf("%s: poLCA's AIC", label), searched_polca$aic, fit$aic,
                   tolerance = 1e-5),
    compare_values(sprintf("%s: residual degrees of freedom", label), searched_polca$resid.df,
                   residual_df, tolerance = 0.5),
    .polca_inference_rows(label, fit, reference$codes, reference$P, bounded, point, errors,
                          labelled, observed),
    search)

  harness <- harness_reference(
    log_likelihood = point$llik, n_parameters = searched_polca$npar,
    posteriors = point$posterior, responses = labelled,
    standard_errors = list(probability = .polca_label(errors$probability, reference$codes,
                                                      observed, levels_of)),
    profile_probabilities = matrix(reference$P, 1L), group_probabilities = 1,
    group_log_likelihood = stats::setNames(point$by_row, seq_len(nrow(observed))),
    data = observed)
  case <- compare_fit(label, fit, harness,
                      tolerance = c(`residuals$p_value` = 1e-5, `residuals$p_adjusted` = 1e-5),
                      extra = extra)
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("poLCA %s", utils::packageVersion("poLCA")),
                                dataset = sprintf("simulated: %s", fixture$name),
                                precision = "machine")
  list(case = case, fit = fit, searched = searched, frame = frame)
}

# ---------------------------------------------------------------------------
# Category reversal
# ---------------------------------------------------------------------------

#' Every output of the boundary fit against the same fit with categories reversed
#'
#' The reversed data are fitted from the boundary fit's own solution with each
#' item's columns reversed, so that both fits describe one point and one class
#' order; the reversed case's own comparison against poLCA has already shown
#' that this point is the reversed data's maximum. Every table must then be the
#' same with the same category labels, except the thresholds, which are the
#' logits of cumulative probabilities and so depend on the order: the
#' cumulative probability up to label `j` in reversed order is one minus the
#' original cumulative probability up to the label before `j`, so the
#' reversed threshold at `j` is minus the original threshold one label down,
#' and the label that comes last in reversed order (the original first) has
#' none.
#'
#' @param original,reversed Results of [.hard_case()] for the two fixtures.
#' @return A `harness_case`.
.hard_reversal <- function(original, reversed) {
  fit <- original$fit
  start <- unclass(starting_values(fit))
  start$response_probabilities <- lapply(start$response_probabilities, function(block) {
    block[, rev(colnames(block)), drop = FALSE]
  })
  frame <- original$frame
  frame[fit$vars] <- lapply(frame[fit$vars], function(column) {
    factor(column, levels = rev(levels(column)))
  })
  stopifnot("the reversed case must fit the same answers in reverse order" =
              identical(lapply(frame[fit$vars], levels), lapply(reversed$frame[fit$vars], levels)),
            "the reversed case must fit the same answers" =
              identical(lapply(frame[fit$vars], as.character),
                        lapply(reversed$frame[fit$vars], as.character)))
  refit <- multilpa(frame, vars = fit$vars, id = "unit", n_profiles = fit$n_profiles,
                    n_group_classes = 1L, categorical = fit$vars, n_starts = 1L, tol = 1e-16,
                    max_iter = 20000L, start = start)
  stopifnot("the reversed fit must converge" = isTRUE(refit$converged))
  expected <- fit_tables(fit)
  expected$responses <- .hard_reverse_thresholds(expected$responses, fit$categorical_levels)
  label <- "hard: boundary data, categories reversed, against the original order"
  case <- .invariance_compare(label, refit, expected)
  # The point the refit reached is the reversed case's own fit from poLCA's
  # estimates on the reversed data, so it is the maximum compared there.
  same_point <- compare_values(
    sprintf("%s: log-likelihood against the reversed case's fit from poLCA's estimates", label),
    as.numeric(logLik(reversed$fit)), as.numeric(logLik(refit)), tolerance = 1e-9)
  case$comparisons <- rbind(.with_output_columns(case$comparisons),
                            data.frame(same_point, source = "invariance"))
  case$comparisons$source <- "invariance"
  case$comparisons$dataset <- "simulated: boundary"
  case$comparisons$precision <- "machine"
  case
}

#' Carry the thresholds of a responses table through a category reversal
#' @param responses The responses table of the original fit.
#' @param levels The fit's category levels, in original order.
#' @return The table with `threshold` as the reversed order defines it.
.hard_reverse_thresholds <- function(responses, levels) {
  # The original threshold of the label one place down in the original order;
  # the first original label has none below it.
  below <- vapply(seq_len(nrow(responses)), function(row) {
    item_levels <- levels[[responses$indicator[[row]]]]
    position <- match(responses$category[[row]], item_levels)
    if (position == 1L) return(NA_real_)
    hit <- responses$profile == responses$profile[[row]] &
      responses$indicator == responses$indicator[[row]] &
      responses$category == item_levels[[position - 1L]]
    -responses$threshold[hit]
  }, numeric(1))
  responses$threshold <- below
  responses
}
