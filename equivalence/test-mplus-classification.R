# The package's own entropy and AvePP tables against the values Mplus prints.
# The posteriors themselves are compared to 1e-5 elsewhere; these tests check
# what is computed from them: the formula, the rows, and which level is which.

test_that("single-level entropy and AvePP reproduce published Mplus examples", {
  invisible(lapply(c("7.9", "7.10"), function(example) {
    fixture <- readRDS(equivalence_fixture("mplus", paste0("public-", example, ".rds")))
    fit <- multilpa(fixture$data, c("y1", "y2", "y3", "y4"), "cluster", 2L, 1L,
                    variance_model = fixture$variance_model, n_starts = 20L,
                    max_iter = 1000L, tol = 1e-12, seed = 912)
    theirs <- mplus_classification_single(
      file.path("mplus", "public", paste0("ex", example, ".html")), 2L)
    # Classes are matched by their posterior-based counts, which Mplus prints
    # in its own label order.
    order <- order(fit$effective_profile_counts)[rank(theirs$counts)]
    expect_lt(max(abs(fit$effective_profile_counts[order] - theirs$counts)), 1e-3)
    mine <- latents_classification(fit, "individuals")
    discrepancy <- classification_discrepancy(theirs, mine, order)
    expect_true(all(discrepancy < mplus_printed_tolerance),
                info = paste(example, paste(signif(discrepancy, 3), collapse = ", ")))
    # A shifted value must fail: the comparison can detect a real error.
    mine$entropy <- mine$entropy + 0.002
    expect_false(all(classification_discrepancy(theirs, mine, order) <
                       mplus_printed_tolerance))
  }))
})

test_that("two-level FIML entropy and AvePP reproduce genuine Mplus at each level", {
  cases <- c(equal = "equal", varying = "full")
  invisible(lapply(names(cases), function(variance_model) {
    reference <- readRDS(equivalence_fixture("mplus",
      paste0("twolevel-missing-full-", variance_model, ".rds")))
    fit <- multilpa(reference$data, c("y1", "y2"), "clus", 2, 2,
                    covariance_model = "full", missing = "fiml",
                    variance_model = variance_model, n_starts = 4, seed = 20260917,
                    tol = 1e-13, max_iter = 10000)
    theirs <- mplus_classification_twolevel(
      file.path("mplus", "extensions-missing-full", paste0(cases[[variance_model]], ".out")),
      2L, 2L)
    # The fixture holds Mplus's own estimates in its label order.
    profile_order <- order(fit$means[, "y1"])[rank(reference$means[, 1])]
    group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
      rank(reference$profile_probabilities[, 1])]
    individuals <- classification_discrepancy(
      theirs$individuals, latents_classification(fit, "individuals"), profile_order)
    groups <- classification_discrepancy(
      theirs$groups, latents_classification(fit, "groups"), group_order)
    expect_true(all(individuals < mplus_printed_tolerance),
                info = paste(variance_model, "individuals", paste(signif(individuals, 3), collapse = ", ")))
    expect_true(all(groups < mplus_printed_tolerance),
                info = paste(variance_model, "groups", paste(signif(groups, 3), collapse = ", ")))
    # The levels are not interchangeable: Mplus's within-level AvePP does not
    # match this package's group level.
    expect_false(isTRUE(all.equal(theirs$individuals$avepp,
                                  latents_classification(fit, "groups")$avepp[group_order],
                                  tolerance = 1e-6)))
  }))
})

test_that("the parser reads a two-level output's levels by name", {
  theirs <- mplus_classification_twolevel(
    file.path("mplus", "extensions-missing-full", "equal.out"), 2L, 2L)
  expect_equal(theirs$individuals$entropy, 0.999)
  expect_equal(theirs$groups$entropy, 1.000)
  expect_equal(theirs$individuals$avepp, c(1.000, 0.999))
})
