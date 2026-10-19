# Local dependence between two categorical indicators: the pair is modelled
# jointly within each class (Latent GOLD's `a <-> b | Cluster`).

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

test_that("a jointly modelled pair agrees with Latent GOLD's class-specific association", {
  skip_on_cran()
  # Latent GOLD 6.1 on these data, `a <-> b | Cluster`: log likelihood
  # -3614.1388 with 11 parameters, class sizes 0.5093 and 0.4907, and these
  # marginal probabilities of `a` and `b` being "yes" in each class.
  fit <- joint_fit(dependent_items())
  within <- function(actual, printed) {
    expect_true(all(abs(actual - printed) <= 5e-5),
                info = paste(format(actual, digits = 8), collapse = " "))
  }
  within(fit$log_likelihood, -3614.1388)
  expect_equal(fit$n_parameters, 11)
  sizes <- colMeans(fit$subject_posteriors)
  order <- order(-sizes)
  within(sizes[order], c(0.5093, 0.4907))
  table <- get_results(fit, "local_dependence")
  yes_a <- aggregate(probability ~ profile, subset(table, category_1 == "yes"), sum)
  yes_b <- aggregate(probability ~ profile, subset(table, category_2 == "yes"), sum)
  within(yes_a$probability[order], c(0.6914, 0.2184))
  within(yes_b$probability[order], c(0.6827, 0.2750))
})

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
  expect_named(table, c("pair", "profile", "item_1", "category_1", "item_2",
                        "category_2", "probability", "independent", "ratio"))
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
  # Rows observing one item of the pair lose the pair, and say so.
  gappy <- data
  gappy$b[1:5] <- NA
  expect_warning(lca(gappy, c("a", "b", "c", "e"), 2, n_starts = 1, seed = 1,
                     missing = "fiml", local_dependence = "a:b"),
                 class = "latents_partial_pair")
})
