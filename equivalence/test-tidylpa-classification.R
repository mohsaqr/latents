# The package's entropy and AvePP against tidyLPA's own formulas, applied to
# the same posteriors. Same rows and posteriors in, so this isolates the
# formulas. An independent tidyLPA fit is not compared: mclust's single EM run
# from its hierarchical start can stop at a different maximum, and its
# classification statistics then describe a different solution.

test_that("relative entropy and AvePP equal tidyLPA's on the same posteriors", {
  skip_if_not_installed("tidyLPA")
  vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy",
            "self_regulation", "test_anxiety")
  invisible(lapply(3:5, function(k) {
    fit <- lpa(srl, vars, n_profiles = k, model = "VVI", n_starts = 10, seed = 1)
    z <- fit$subject_posteriors
    most_likely <- max.col(z, ties.method = "first")
    # tidyLPA (calc_fitindices for an Mclust model) adds 1e-12 inside the log.
    tidy_entropy <- 1 + sum(rowSums(z * log(z + 1e-12))) / (nrow(z) * log(ncol(z)))
    tidy_avepp <- diag(tidyLPA:::avgprobs_mostlikely(z, most_likely))
    mine <- latents_classification(fit, "individuals")
    expect_equal(mine$entropy, tidy_entropy, tolerance = 1e-10)
    expect_equal(mine$avepp, unname(tidy_avepp), tolerance = 1e-12)
    # tidyLPA's prob_min and prob_max are the range of the same diagonal.
    expect_equal(range(mine$avepp), range(tidy_avepp), tolerance = 1e-12)
  }))
})

test_that("tidyLPA's formulas take every row, so they agree on complete data only", {
  skip_if_not_installed("tidyLPA")
  # With rows that observe nothing, this package leaves them out of the
  # classification statistics (they were never classified); tidyLPA's mclust
  # engine cannot fit missing data at all, so its formula on all rows is the
  # reference only when every row is observed.
  data <- srl[c("intrinsic_value", "test_anxiety")]
  fit <- lpa(data, names(data), n_profiles = 2, model = "VVI", n_starts = 5, seed = 1)
  expect_true(all(fit$informative))
  z <- fit$subject_posteriors
  expect_equal(latents_classification(fit, "individuals")$entropy,
               1 + sum(rowSums(z * log(z + 1e-12))) / (nrow(z) * log(ncol(z))),
               tolerance = 1e-10)
})
