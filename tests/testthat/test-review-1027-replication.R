# Review 2026-10-27: how often the starts reached the best log likelihood is
# shown where a user looks, not only in summary(). print() used to show the
# winning start's index as "best start: 5/10", which reads as a count.

replication_fit <- function(n_starts = 10L) {
  withr::local_seed(20261027)
  data <- data.frame(id = rep(seq_len(100L), each = 4L),
                     a = factor(sample(3L, 400L, TRUE)), b = factor(sample(3L, 400L, TRUE)),
                     c = factor(sample(3L, 400L, TRUE)))
  suppressWarnings(multilca(data, c("a", "b", "c"), id = "id", n_profiles = 2,
                            n_group_classes = 2, n_starts = n_starts, seed = 1))
}

test_that("the replication table counts starts within tolerance of the best", {
  fit <- replication_fit()
  table <- get_results(fit, "replication")
  expect_s3_class(table, "data.frame")
  expect_identical(nrow(table), 1L)
  expect_identical(table$n_starts, 10L)
  expect_identical(table$n_best_replicated, as.integer(fit$n_best_replicated))
  expect_equal(table$share, fit$n_best_replicated / 10)
  # The count is recomputed from the stored starts by the documented rule.
  finished <- is.finite(fit$starts$log_likelihood)
  expect_identical(table$n_best_replicated,
                   sum(finished & abs(fit$starts$log_likelihood - fit$log_likelihood) <=
                         fit$replication_tolerance))
  expect_identical(table$best_start, as.integer(fit$best_start))
})

test_that("print and diagnostics state the count, not the winning start's index", {
  fit <- replication_fit()
  printed <- utils::capture.output(print(fit))
  expected <- sprintf("best log likelihood reached by %d of 10 starts", fit$n_best_replicated)
  expect_true(any(grepl(expected, printed, fixed = TRUE)))
  expect_false(any(grepl("best start: ", printed, fixed = TRUE)))
  checked <- diagnostics(fit, plots = FALSE)
  shown <- utils::capture.output(print(checked))
  expect_true(any(grepl(sprintf("%d of 10 starts", fit$n_best_replicated), shown,
                        fixed = TRUE)))
  expect_identical(get_results(checked, "replication"), get_results(fit, "replication"))
})

test_that("a best solution found by one start of several is flagged", {
  fit <- replication_fit()
  fit$n_best_replicated <- 1L
  expect_match(.multilpa_replication_phrase(fit), "found once: refit with more starts")
  fit$n_best_replicated <- NA_integer_
  expect_match(.multilpa_replication_phrase(fit), "^best start: ")
  # A diagnostics object saved before the table existed still answers.
  checked <- diagnostics(replication_fit(), plots = FALSE)
  checked$replication <- NULL
  expect_identical(nrow(get_results(checked, "replication")), 1L)
})

test_that("a clicked result shows its header before its internals", {
  fit <- replication_fit()
  # A host that renders str() and then print() (Beatrina, RStudio) shows the
  # first characters; they are the result's header, not its list plumbing.
  shown <- c(utils::capture.output(str(fit)), utils::capture.output(print(fit)))
  expect_match(shown[[1L]], "^<multilpa> Two-level latent class analysis")
  first_900 <- substr(paste(shown, collapse = "\n"), 1L, 900L)
  expect_match(first_900, "Class sizes")
  expect_match(first_900, "reached by")
  expect_lt(length(utils::capture.output(str(fit))), 3L)
  # The raw list stays reachable.
  expect_gt(length(utils::capture.output(str(unclass(fit)))), 50L)
  # A wide response table: one row per indicator category, one column per class.
  wide <- .multilpa_wide_responses(fit)
  expect_named(wide, c("indicator", "category", "profile_1", "profile_2"))
  expect_identical(nrow(wide), 9L)
  expect_equal(unname(colSums(wide[c("profile_1", "profile_2")])), c(3, 3), tolerance = 1e-12)
})

test_that("the replication table counts converged and failed starts", {
  fit <- replication_fit()
  table <- get_results(fit, "replication")
  starts <- get_results(fit, "starts")
  expect_identical(table$n_converged, sum(is.finite(starts$log_likelihood) & starts$converged))
  expect_identical(table$n_failed, sum(!is.finite(starts$log_likelihood)))
  expect_match(.multilpa_replication_phrase(fit), sprintf("\\(%d converged", table$n_converged))
  # A failed and an unconverged start are named in the phrase.
  fit$starts$log_likelihood[2L] <- NA_real_
  fit$starts$converged[3L] <- FALSE
  expect_match(.multilpa_replication_phrase(fit), "1 stopped at max_iter; 1 failed")
  # One parenthesis, however many notes.
  fit$n_best_replicated <- 1L
  expect_identical(lengths(regmatches(.multilpa_replication_phrase(fit),
                                      gregexpr("\\(", .multilpa_replication_phrase(fit)))), 1L)
})

test_that("a class with no assigned member does not turn the lowest posterior into NA", {
  fit <- replication_fit()
  checked <- diagnostics(fit, plots = FALSE)
  averages <- checked$average_posteriors
  own <- averages$assigned_class == averages$class & averages$level == "individuals"
  emptied <- which(own)[[1L]]
  checked$average_posteriors$average_posterior[emptied] <- NA_real_
  # The line-per-diagnostic print of an object saved before the checks.
  checked$checks <- NULL
  shown <- utils::capture.output(print(checked))
  line <- grep("Lowest avg posterior", shown, value = TRUE)
  expect_false(grepl(" NA", line, fixed = TRUE))
  expect_match(line, sprintf("class %s: no unit assigned", averages$class[emptied]))
})
