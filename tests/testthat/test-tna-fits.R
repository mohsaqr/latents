# Transition networks from a multilpa() / multilca() fit made with `time =`:
# the observed transitions between modal profiles, checked against a count
# made independently from the assignments.

.sequence_fit <- function(covariates = FALSE) {
  set.seed(1)
  data <- data.frame(g = rep(seq_len(30), each = 6), t = rep(seq_len(6), 30))
  data$y1 <- sample(1:2, nrow(data), TRUE)
  data$y2 <- ifelse(stats::runif(nrow(data)) < 0.85, data$y1, 3L - data$y1)
  data$y3 <- ifelse(stats::runif(nrow(data)) < 0.85, data$y1, 3L - data$y1)
  data$z <- rep(stats::rnorm(30), each = 6)
  # Shuffle the rows so the network cannot lean on the input order.
  data <- data[sample(nrow(data)), ]
  row.names(data) <- NULL
  fit <- multilca(data, vars = c("y1", "y2", "y3"), id = "g", time = "t",
                  n_profiles = 2, n_group_classes = 2, n_starts = 2, seed = 1,
                  profile_covariates = if (covariates) "z" else character())
  list(fit = fit, data = data)
}

test_that("the network is the row-normalised count of consecutive profiles", {
  skip_if_not_installed("tna")
  pair <- .sequence_fit()
  assigned <- get_results(pair$fit, "assignments")
  assigned <- assigned[order(assigned$g, assigned$t), ]
  same_group <- head(assigned$g, -1L) == tail(assigned$g, -1L)
  from <- head(assigned$profile, -1L)[same_group]
  to <- tail(assigned$profile, -1L)[same_group]
  counts <- table(factor(from, 1:2), factor(to, 1:2))
  expected <- unclass(counts / rowSums(counts))
  network <- get_tna(pair$fit)
  expect_equal(unname(network$weights), unname(matrix(expected, 2L, 2L)))
  first <- assigned$profile[!duplicated(assigned$g)]
  expect_equal(unname(network$inits),
               as.numeric(table(factor(first, 1:2))) / length(first))
})

test_that("the grouped networks follow each group's modal group class", {
  skip_if_not_installed("tna")
  pair <- .sequence_fit()
  networks <- get_group_tna(pair$fit, label = "Trajectory")
  expect_s3_class(networks, "group_tna")
  classes <- get_results(pair$fit, "group_posteriors")
  modal <- sort(unique(classes$group_class[classes$modal]))
  expect_setequal(names(networks), sprintf("Trajectory %d", modal))
  expect_error(get_group_tna(pair$fit, label = c("a", "b")), "single string")
})

test_that("a covariate fit gives networks and draws its response probabilities", {
  skip_if_not_installed("tna")
  fit <- suppressWarnings(.sequence_fit(covariates = TRUE)$fit)
  expect_s3_class(get_tna(fit), "tna")
  expect_s3_class(get_group_tna(fit), "group_tna")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, what = "responses"))
})

test_that("a fit without an order has no network", {
  skip_if_not_installed("tna")
  data <- .sequence_fit()$data
  fit <- multilca(data, vars = c("y1", "y2", "y3"), id = "g", n_profiles = 2,
                  n_group_classes = 1, n_starts = 1, seed = 1)
  expect_error(get_tna(fit), class = "latents_no_time")
  expect_error(get_group_tna(fit), class = "latents_no_time")
})
