# Attribution is part of the evidence contract: component formulas and local
# oracles must stay visible without becoming independent full-model fits.

test_that("component and local-oracle checks cannot become model sources", {
  rules <- utils::read.csv(file.path("verification-matrix", "test-attribution.csv"))
  selected <- function(file, title) {
    own <- rules[rules$file == file, ]
    hit <- which(vapply(own$pattern, grepl, logical(1), x = title, perl = TRUE))
    expect_gt(length(hit), 0L)
    own[hit[[1L]], ]
  }
  expect_identical(selected("test-weighted-logits.R",
    "binary soft counts give the weighted logistic regression")$kind, "component check")
  expect_identical(selected("test-tidylpa-classification.R",
    "relative entropy and AvePP equal tidyLPA's on the same posteriors")$kind,
    "component check")
  expect_identical(selected("test-features-pair-covariance.R",
    "standard errors and the score agree with numerical derivatives")$kind,
    "independent code")
  expect_identical(selected("test-features-vlmr.R",
    "kit p-values do not move with the indicators' units and origin")$kind, "internal")
  expect_identical(selected("test-features-pair-covariance.R",
    "case_01: latents agrees with latentgold")$reference, "Latent GOLD")
  expect_identical(selected("test-features-pair-covariance.R",
    "case_01: latents agrees with mplus_fixed")$reference, "Mplus")
  expect_identical(selected("test-features-direct-effects.R",
    "log likelihood, parameter count and BIC agree with both programs")$fn, "lpa")
  expect_identical(selected("test-features-joint-items.R",
    "the two-level model with a joined pair is compared in full")$fn, "multilca")
  references <- utils::read.csv(file.path("verification-matrix", "harness-references.csv"))
  local_em <- startsWith(references$pattern, "^EM with ")
  expect_equal(sum(local_em), 2L)
  expect_true(all(references$kind[local_em] == "independent code"))
})

test_that("Latent GOLD checks credit the public function that actually runs", {
  cases <- readRDS(file.path("latentgold-coverage", "kit", "cases.rds"))
  comparison <- utils::read.csv(file.path("latentgold-coverage", "comparison.csv"))
  kind <- vapply(cases, `[[`, character(1), "kind")
  lca_cases <- names(kind)[kind == "lca_enumeration"]
  lpa_cases <- names(kind)[kind == "lpa_enumeration"]
  expressions <- parse(file.path("latentgold-coverage", "compare-functions.R"))
  mapping <- Filter(function(expression) is.call(expression) &&
    identical(expression[[1L]], as.name("<-")) &&
    identical(expression[[2L]], as.name(".case_function")), as.list(expressions))
  expect_length(mapping, 1L)
  functions <- eval(mapping[[1L]][[3L]], envir = baseenv())
  expect_identical(unname(functions[c("lca_enumeration", "lpa_enumeration")]),
                   c("lca", "lpa"))
  expect_length(lca_cases, 4L)
  expect_length(lpa_cases, 4L)
  expect_true(all(comparison$fn[comparison$case %in% lca_cases] == "lca"))
  expect_true(all(comparison$fn[comparison$case %in% lpa_cases] == "lpa"))
  expect_false(any(comparison$fn == "enumeration"))
  rules <- utils::read.csv(file.path("verification-matrix", "test-attribution.csv"))
  markov <- rules[rules$file == "test-latentgold-step3-markov.R", ]
  expect_identical(markov$fn, "three_step")
  datasets <- utils::read.csv(file.path("verification-matrix", "test-datasets.csv"))
  markov_data <- datasets[grepl("test-latentgold-step3-markov.R", datasets$location,
                               fixed = TRUE), ]
  expect_gt(nrow(markov_data), 0L)
  expect_true(all(markov_data$fn == "three_step"))
})
