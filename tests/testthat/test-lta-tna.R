skip_on_cran()
skip_if_not_installed("tna")

# A general lta() fit with homogeneous first-order transitions has one
# transition matrix per group class, so it hands its estimates to tna.

.lta_tna_data <- function(seed = 5L, n_groups = 120L, per = 6L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = per),
                      t = rep(seq_len(per), n_groups))
  # Half the units never move; the other half switch state often, so a
  # second group class and a stayer class are both in the data.
  state <- unlist(lapply(seq_len(n_groups), \(i) {
    if (i <= n_groups / 2) rep(i %% 2 + 1, per) else
      cumsum(stats::runif(per) > 0.4) %% 2 + 1
  }))
  rating <- function(shift) {
    as.integer(cut(stats::rnorm(length(state), ifelse(state == 2, shift, -shift)),
                   c(-Inf, -1, 0, 1, Inf)))
  }
  frame$a <- rating(1.3)
  frame$b <- rating(1.3)
  frame$c <- rating(0.9)
  frame
}

.lta_tna_fit <- function(...) {
  lta(.lta_tna_data(), c("a", "b", "c"), id = "g", time = "t", n_profiles = 2,
      ordinal = c("a", "b", "c"), n_starts = 3, seed = 1, ...)
}

.table_matrix <- function(fit, class) {
  rows <- subset(get_results(fit, "transitions"), group_class == class)
  states <- sprintf("profile_%d", seq_len(fit$n_profiles))
  matrix(rows$probability[order(match(rows$from, states), match(rows$to, states))],
         length(states), byrow = TRUE, dimnames = list(states, states))
}

test_that("get_tna() of an ordinal lta carries its estimated matrix and starts", {
  fit <- .lta_tna_fit()
  expect_s3_class(fit, "multilpa_lta")
  network <- get_tna(fit)
  expect_s3_class(network, "tna")
  expect_equal(unclass(network$weights), .table_matrix(fit, "group_class_1"),
               tolerance = 1e-12, ignore_attr = TRUE)
  starts <- get_results(fit, "initial")
  expect_equal(unname(network$inits), starts$probability, tolerance = 1e-12)
  # Rows of a transition matrix sum to one.
  expect_equal(unname(rowSums(network$weights)), rep(1, 2), tolerance = 1e-12)
})

test_that("get_tna() pools two classes by their expected moves", {
  fit <- .lta_tna_fit(n_group_classes = 2)
  network <- get_tna(fit)
  rows <- get_results(fit, "transitions")
  rows$moves <- rows$probability * rows$origin_expected_count
  pooled <- aggregate(cbind(moves, origin_expected_count) ~ from + to, data = rows,
                      FUN = sum)
  pooled$probability <- pooled$moves / pooled$origin_expected_count
  pooled <- pooled[order(pooled$from, pooled$to), ]
  expect_equal(as.vector(t(unclass(network$weights))), pooled$probability,
               tolerance = 1e-10)
})

test_that("get_group_tna() of a mover-stayer fit adds an identity network", {
  fit <- .lta_tna_fit(mover_stayer = TRUE)
  networks <- get_group_tna(fit)
  expect_s3_class(networks, "group_tna")
  expect_identical(names(networks), c("Group class 1", "Stayers"))
  expect_equal(unclass(networks[["Stayers"]]$weights), diag(2), ignore_attr = TRUE)
  expect_equal(unclass(networks[["Group class 1"]]$weights),
               .table_matrix(fit, "group_class_1"), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_error(get_tna(fit), class = "latents_unsupported_tna")
})

test_that("transitions that vary by occasion or covariate, or are second order, are refused", {
  frame <- .lta_tna_data()
  frame$shift <- as.integer(frame$t > 3)
  by_covariate <- lta(frame, c("a", "b", "c"), id = "g", time = "t",
                      n_profiles = 2, ordinal = c("a", "b", "c"),
                      transition_covariates = "shift", n_starts = 2, seed = 1)
  expect_error(get_tna(by_covariate), class = "latents_unsupported_tna")
  expect_error(get_group_tna(by_covariate), class = "latents_unsupported_tna")
  by_occasion <- lta(frame, c("a", "b", "c"), id = "g", time = "t",
                     n_profiles = 2, ordinal = c("a", "b", "c"),
                     transitions = "occasion", n_starts = 2, seed = 1)
  expect_error(get_tna(by_occasion), class = "latents_unsupported_tna")
  second_order <- lta(frame, c("a", "b", "c"), id = "g", time = "t",
                      n_profiles = 2, ordinal = c("a", "b", "c"), order = 2,
                      n_starts = 2, seed = 1)
  expect_error(get_group_tna(second_order), class = "latents_unsupported_tna")
})
