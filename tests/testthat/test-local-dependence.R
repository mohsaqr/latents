# Local dependence between two categorical indicators: the pair is modelled
# jointly within each class (Latent GOLD's `a <-> b | Cluster`).
# The comparisons with independent references are in tests/equivalence/test-local-dependence.R.

dependent_items <- function(seed = 5, n = 1500) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE)
  # Two items driven by one latent response within each class, so they stay
  # associated given the class; two more that are independent given it.
  latent <- stats::rnorm(n)
  first <- ifelse(class == 1, latent > -0.5, latent > 0.8)
  second <- ifelse(class == 1, latent + stats::rnorm(n, 0, 0.6) > -0.5,
                   latent + stats::rnorm(n, 0, 0.6) > 0.8)
  yes_no <- function(value) factor(value, c(FALSE, TRUE), c("no", "yes"))
  data.frame(a = yes_no(first), b = yes_no(second),
             c = factor(stats::rbinom(n, 1, c(0.8, 0.2)[class]), 0:1, c("no", "yes")),
             e = factor(stats::rbinom(n, 1, c(0.75, 0.3)[class]), 0:1, c("no", "yes")))
}

joint_fit <- function(data, pairs = "a:b") {
  quietly(lca(data, c("a", "b", "c", "e"), 2, n_starts = 10, seed = 1, tol = 1e-12,
              local_dependence = pairs))
}

test_that("the joint model does not depend on the order of the pair or how it is joined", {
  skip_on_cran()
  data <- dependent_items()
  forward <- joint_fit(data, "a:b")
  backward <- joint_fit(data, "b:a")
  expect_equal(backward$log_likelihood, forward$log_likelihood, tolerance = 1e-8)
  # Joining the two items by hand into one categorical indicator is the model.
  by_hand <- transform(data, ab = interaction(a, b, sep = ":", lex.order = TRUE))
  manual <- quietly(lca(by_hand, c("ab", "c", "e"), 2, n_starts = 10, seed = 1,
                        tol = 1e-12))
  expect_equal(manual$log_likelihood, forward$log_likelihood, tolerance = 1e-8)
  expect_identical(manual$n_parameters, forward$n_parameters)
})

test_that("the dependence table holds each class's joint table and its independence", {
  skip_on_cran()
  fit <- joint_fit(dependent_items())
  table <- get_results(fit, "local_dependence")
  expect_named(table, c("pair", "kind", "profile", "item_1", "category_1", "item_2",
                        "category_2", "probability", "independent", "ratio",
                        "covariance", "correlation"))
  expect_identical(nrow(table), 8L)
  totals <- aggregate(cbind(probability, independent) ~ profile, table, sum)
  expect_equal(totals$probability, c(1, 1), tolerance = 1e-10)
  expect_equal(totals$independent, c(1, 1), tolerance = 1e-10)
  expect_equal(table$ratio, table$probability / table$independent, tolerance = 1e-12)
  # The simulated items agree more often than independence implies.
  same <- subset(table, category_1 == category_2)
  expect_true(all(same$ratio > 1))
  expect_false("local_dependence" %in% names(get_results(
    quietly(lca(dependent_items(), c("a", "b", "c", "e"), 2, n_starts = 2, seed = 1)),
    "all")))
  expect_true("local_dependence" %in% names(get_results(fit, "all")))
})

test_that("a joint model predicts from the caller's own columns", {
  skip_on_cran()
  data <- dependent_items()
  fit <- joint_fit(data)
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-10)
  expect_identical(fit$vars, c("a:b", "c", "e"))
  expect_identical(as.character(fit$call[[1L]]), "lca")
})

test_that("pairs that cannot be modelled are refused by class", {
  data <- dependent_items(n = 300)
  fit_with <- function(pairs, ...) {
    lca(data, c("a", "b", "c", "e"), 2, n_starts = 1, seed = 1,
        local_dependence = pairs, ...)
  }
  expect_error(fit_with("ab"), class = "latents_bad_argument")
  expect_error(fit_with("a:a"), class = "latents_bad_argument")
  expect_error(fit_with(c("a:b", "b:c")), class = "latents_bad_argument")
  expect_error(fit_with("a:z"), class = "latents_bad_argument")
  continuous <- transform(data, score = stats::rnorm(nrow(data)))
  expect_error(multilpa(continuous, c("a", "score", "c"), NULL, 2, 1,
                        categorical = c("a", "c"), local_dependence = "a:score",
                        n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  plain <- quietly(lca(data, c("a", "b", "c", "e"), 2, n_starts = 1, seed = 1))
  expect_error(get_results(plain, "local_dependence"),
               class = "latents_no_local_dependence")
  expect_error(.latents_join_pairs(data["a"], .latents_parse_pairs("a:b")),
               class = "latents_bad_data")
  # A colon in a category would merge ("x:y", "z") with ("x", "y:z").
  colons <- data.frame(a = factor(rep(c("x:y", "x", "u"), each = 20)),
                       b = factor(rep(c("z", "y:z", "v"), each = 20)),
                       c = factor(rep(0:1, 30)))
  expect_error(lca(colons, c("a", "b", "c"), 1, local_dependence = "a:b", n_starts = 1,
                   seed = 1), class = "latents_bad_data")
  # Rows observing one item of the pair keep that item (review 2026-10-26
  # round 2, C2-2); the old expectation here, a latents_partial_pair warning
  # that the pair was dropped, was the bias itself.
  gappy <- data
  gappy$b[1:5] <- NA
  expect_no_warning(gappy_fit <- lca(gappy, c("a", "b", "c", "e"), 2, n_starts = 1, seed = 1,
                                     missing = "fiml", local_dependence = "a:b"),
                    class = "latents_partial_pair")
  expect_false(is.null(attr(gappy_fit$categorical_data, "latents_partial")))
})

correlated_pair <- function(seed = 9, n = 1200) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE)
  correlation <- c(0.6, -0.5)[class]
  first <- stats::rnorm(n)
  second <- correlation * first + sqrt(1 - correlation^2) * stats::rnorm(n)
  data.frame(y1 = c(0, 2.5)[class] + first, y2 = c(0, 2)[class] + second,
             y3 = stats::rnorm(n, c(0, 1.5)[class]), y4 = stats::rnorm(n, c(0, 1)[class]))
}

test_that("one pair of two indicators frees every covariance, as the full model does", {
  skip_on_cran()
  data <- correlated_pair()[c("y1", "y2")]
  pattern <- quietly(lpa(data, c("y1", "y2"), 2, n_starts = 10, seed = 1, tol = 1e-12,
                         local_dependence = "y1:y2"))
  full <- quietly(lpa(data, c("y1", "y2"), 2, n_starts = 10, seed = 1, tol = 1e-12,
                      model = "VVV"))
  expect_equal(pattern$log_likelihood, full$log_likelihood, tolerance = 1e-8)
  expect_equal(pattern$n_parameters, full$n_parameters)
  reversed <- quietly(lpa(correlated_pair(), paste0("y", 1:4), 2, n_starts = 10, seed = 1,
                          tol = 1e-12, local_dependence = "y2:y1"))
  forward <- quietly(lpa(correlated_pair(), paste0("y", 1:4), 2, n_starts = 10, seed = 1,
                         tol = 1e-12, local_dependence = "y1:y2"))
  expect_equal(reversed$log_likelihood, forward$log_likelihood, tolerance = 1e-8)
  # With one covariance matrix shared by the classes, the pair adds one.
  equal <- quietly(lpa(correlated_pair(), paste0("y", 1:4), 2, n_starts = 5, seed = 1,
                       variance_model = "equal", local_dependence = "y1:y2"))
  independent <- quietly(lpa(correlated_pair(), paste0("y", 1:4), 2, n_starts = 5,
                             seed = 1, model = "EEI"))
  expect_equal(equal$n_parameters, independent$n_parameters + 1)
})

test_that("continuous pairs that the model cannot hold are refused by class", {
  data <- transform(correlated_pair(n = 300), group = rep(1:30, each = 10),
                    item = factor(stats::rbinom(300, 1, 0.5)))
  vars <- paste0("y", 1:4)
  expect_error(lpa(data, vars, 2, model = "VVV", local_dependence = "y1:y2",
                   n_starts = 1, seed = 1), class = "latents_bad_argument")
  expect_error(lpa(data, vars, 2, covariance_model = "full", local_dependence = "y1:y2",
                   n_starts = 1, seed = 1), class = "latents_bad_argument")
  expect_error(multilpa(data, c(vars, "item"), NULL, 2, 1, categorical = "item",
                        local_dependence = "y1:item", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars[2:4], "group", 2, 2, profile_covariates = "y1",
                        local_dependence = "y2:y3", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
})

test_that("refits of a fit with chosen covariances keep the pattern", {
  skip_on_cran()
  data <- correlated_pair(n = 300)
  fit <- quietly(lpa(data, c("y1", "y2", "y3"), 2, n_starts = 3, seed = 1,
                     local_dependence = "y1:y2"))
  arguments <- .multilpa_refit_arguments(fit)
  expect_identical(arguments$local_dependence, "y1:y2")
  # A single-level fit is refitted on its rebuilt frame, which carries the
  # row identifier the refit groups by, as the bootstrap does.
  refit <- quietly(do.call(multilpa, c(list(data = .multilpa_resolve_data(fit, NULL),
                                            n_starts = 3, seed = 1), arguments)))
  expect_equal(refit$n_parameters, fit$n_parameters)
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-6)
  expect_identical(refit$covariance_pattern, fit$covariance_pattern)
})

test_that("the items of a joined pair count as indicators for the three-step guards", {
  set.seed(9)
  data <- data.frame(a = stats::rbinom(300, 1, 0.5), b = stats::rbinom(300, 1, 0.5),
                     c = stats::rbinom(300, 1, 0.5), e = stats::rbinom(300, 1, 0.5),
                     g = rep(seq_len(30), each = 10))
  fit <- quietly(lca(data, c("a", "b", "c", "e"), 2, local_dependence = "a:b",
                     n_starts = 2, seed = 1))
  expect_error(three_step(fit, data, "a", method = "modal"), class = "multilpa_indicator_reused")
  expect_error(r3step(fit, data, "b"), class = "multilpa_indicator_reused")
  expect_error(multilpa(data, c("a", "b", "c", "e"), "g", 2, 2,
                        categorical = c("a", "b", "c", "e"), local_dependence = "a:b",
                        profile_covariates = "a", n_starts = 1, seed = 1),
               class = "multilpa_indicator_reused")
})

test_that("enumeration with chosen covariances crosses the diagonal they extend", {
  skip_on_cran()
  data <- correlated_pair(n = 300)
  vars <- c("y1", "y2", "y3")
  grid <- quietly(enumerate_lpa(data, vars, n_profiles = 1:2, local_dependence = "y1:y2",
                                n_starts = 2, seed = 1))
  table <- as.data.frame(grid)
  expect_identical(nrow(table), 4L)
  expect_true(all(is.na(table$error)))
  expect_setequal(unique(table$model), c("EEI + y1:y2", "VVI + y1:y2"))
  direct <- quietly(lpa(data, vars, 2, local_dependence = "y1:y2", n_starts = 2, seed = 1))
  expect_equal(subset(table, n_profiles == 2 & model == "VVI + y1:y2")$log_likelihood,
               direct$log_likelihood, tolerance = 1e-8)
  expect_error(enumerate_lpa(data, vars, n_profiles = 1:2, local_dependence = "y1:y2",
                             model = "VVV", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(enumerate_lpa(data, vars, n_profiles = 2, local_dependence = "y1:y2",
                             model = "EEI", variance_model = "varying", n_starts = 1,
                             seed = 1), class = "latents_bad_argument")
  # `model = NULL` crosses no structure, and `variance_model` alone is the one.
  single <- as.data.frame(quietly(enumerate_lpa(data, vars, n_profiles = 2, model = NULL,
                                                local_dependence = "y1:y2", n_starts = 1,
                                                seed = 1)))
  expect_identical(single$model, "VVI + y1:y2")
  equal <- as.data.frame(quietly(enumerate_lpa(data, vars, n_profiles = 2,
                                               variance_model = "equal",
                                               local_dependence = "y1:y2", n_starts = 1,
                                               seed = 1)))
  expect_identical(equal$model, "EEI + y1:y2")
  # A pattern fit is labelled by the diagonal it extends and its pairs, not
  # as the full covariance it is estimated through.
  full <- quietly(lpa(data, vars, 2, model = "VVV", n_starts = 2, seed = 1))
  compared <- compare_models(paired = direct, full = full)
  expect_identical(compared$structure, c("VVI + y1:y2", "VVV"))
  expect_output(print(direct), "VVI + y1:y2", fixed = TRUE)
  expect_identical(get_results(direct, "model")$covariance_structure, "VVI + y1:y2")
})
