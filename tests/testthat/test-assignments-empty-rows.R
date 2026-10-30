assignment_empty_case <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_seed(27)
    profile <- c(rep(1L, 90L), rep(2L, 210L))
    data <- data.frame(a = rnorm(300L, c(-3, 3)[profile]),
                       b = rnorm(300L, c(-3, 3)[profile]))
    data[9L, ] <- NA_real_
    data$a[10L] <- NA_real_
    fit <- quietly(lpa(data, c("a", "b"), n_profiles = 2L,
                        missing = "fiml", n_starts = 2L, seed = 1L))
    cache <<- list(fit = fit, data = data, truth = profile)
    cache
  }
})

test_that("assignments preserve rows without classifying empty observations", {
  case <- assignment_empty_case()
  original <- case$fit$subject_posteriors
  assignments <- get_results(case$fit, "assignments", data = case$data)
  fields <- c("profile", "uncertainty", "posterior_profile_1",
              "posterior_profile_2")
  expect_identical(assignments[c("a", "b")], case$data)
  expect_true(all(is.na(assignments[9L, fields])))
  expect_false(anyNA(assignments[10L, fields]))
  expect_equal(sum(!is.na(assignments$profile)), nobs(case$fit))
  classification <- get_results(case$fit, "classification")
  expect_equal(as.integer(table(assignments$profile)), classification$n_modal)
  expect_identical(case$fit$subject_posteriors, original)
  expect_equal(unname(original[9L, ]), as.vector(case$fit$profile_probabilities))
  recovery <- get_results(case$fit, "assignments",
                          data = transform(case$data, truth = case$truth),
                          truth = "truth")
  expect_equal(sum(recovery$n[is.na(recovery$class)]), 1)
  expect_equal(sum(recovery$n[!is.na(recovery$class)]), 299)
})

test_that("both posterior formats and all-table output mark excluded cases", {
  fit <- assignment_empty_case()$fit
  wide <- get_results(fit, "posteriors", format = "wide")
  long <- get_results(fit, "posteriors")
  expect_identical(wide$row, seq_len(300L))
  expect_true(all(is.na(wide[9L, c("profile", "posterior_profile_1",
                                  "posterior_profile_2")])))
  expect_true(all(is.na(long$posterior[long$row == 9L])))
  expect_true(all(is.na(long$modal[long$row == 9L])))
  expect_false(anyNA(long$posterior[long$row == 10L]))
  expect_equal(sum(long$modal, na.rm = TRUE), 299L)
  all_tables <- get_results(fit, "all")
  expect_true(is.na(all_tables$assignments$profile[9L]))
  expect_true(all(is.na(all_tables$posteriors$modal[all_tables$posteriors$row == 9L])))
})

test_that("legacy reporting keeps every assignment; zero weights stay out of counts", {
  case <- assignment_empty_case()
  data <- transform(case$data, weight = rep(1, 300L))
  data$weight[11L] <- 0
  fit <- quietly(lpa(data, c("a", "b"), n_profiles = 2L,
                      weights = "weight", missing = "fiml", drop_empty = FALSE,
                      n_starts = 2L, seed = 1L))
  assigned <- get_results(fit, "assignments")
  expect_false(anyNA(assigned[9L, c("profile", "uncertainty",
                                  "posterior_profile_1", "posterior_profile_2")]))
  # A zero-weight unit keeps its row, posterior and assignment in output but
  # is absent from every count (ML-1).
  expect_false(anyNA(assigned[11L, c("profile", "uncertainty",
                                   "posterior_profile_1", "posterior_profile_2")]))
  expect_false(anyNA(assigned$profile))
  expect_equal(sum(get_results(fit, "classification")$n_modal), nrow(data) - 1L)
})

test_that("zero-weight rows keep their assignment while empty rows lose theirs", {
  case <- assignment_empty_case()
  data <- transform(case$data, weight = rep(1, 300L))
  data$weight[11L] <- 0
  fit <- quietly(lpa(data, c("a", "b"), n_profiles = 2L, weights = "weight",
                      missing = "fiml", n_starts = 2L, seed = 1L))
  fields <- c("profile", "uncertainty", "posterior_profile_1", "posterior_profile_2")
  assigned <- get_results(fit, "assignments")
  expect_true(all(is.na(assigned[9L, fields])))
  # A zero-weight unit is in no count but keeps its row, posterior and
  # assignment in output (ML-1).
  expect_false(anyNA(assigned[c(10L, 11L), fields]))
  posteriors <- get_results(fit, "posteriors")
  expect_false(anyNA(posteriors$posterior[posteriors$row == 11L]))
  expect_true(all(is.na(posteriors$posterior[posteriors$row == 9L])))
  # Counts leave out both the empty row and the zero-weight row.
  expect_identical(nobs(fit), 298L)
  expect_equal(sum(get_results(fit, "classification")$n_modal), 298)
})

test_that("group classification requires an observed member, not every member", {
  withr::local_seed(5)
  data <- data.frame(school = rep(seq_len(30L), each = 8L))
  shift <- c(-2, 2)[1L + (data$school %% 2L)]
  data$a <- rnorm(nrow(data), shift)
  data$b <- rnorm(nrow(data), -shift)
  data[data$school == 7L, c("a", "b")] <- NA_real_
  data[1L, c("a", "b")] <- NA_real_
  fit <- quietly(multilpa(data, c("a", "b"), id = "school",
                          n_profiles = 2L, n_group_classes = 2L,
                          n_starts = 2L, seed = 1L, missing = "fiml"))
  assigned <- get_results(fit, "assignments")
  expect_true(all(is.na(assigned$profile[data$school == 7L])))
  expect_true(all(is.na(assigned$group_class[data$school == 7L])))
  expect_true(is.na(assigned$profile[1L]))
  expect_false(is.na(assigned$group_class[1L]))
  wide <- get_results(fit, "group_posteriors", format = "wide")
  long <- get_results(fit, "group_posteriors")
  expect_true(all(is.na(wide[wide$group == 7L, c("group_class",
                            "posterior_group_class_1", "posterior_group_class_2")])))
  expect_true(all(is.na(long$modal[long$group == 7L])))
  expect_true(all(is.na(long$posterior[long$group == 7L])))
  expect_equal(sum(long$modal, na.rm = TRUE),
               sum(get_results(fit, "classification", level = "groups")$n_modal))
  legacy <- fit
  legacy$drop_empty <- FALSE
  expect_false(anyNA(get_results(legacy, "assignments")$group_class))
  expect_false(anyNA(get_results(legacy, "group_posteriors")$posterior))
})
