# Empty observations must not acquire observed class membership indirectly
# through a descriptive summary, a general transition model, or sensitivity.
empty_consumer_data <- function() {
  withr::local_seed(702)
  data <- data.frame(id = rep(seq_len(24), each = 4), time = rep(1:4, 24))
  shift <- rep(rep(c(-2, 2), each = 12), each = 4)
  data$a <- stats::rnorm(96, shift)
  data$b <- stats::rnorm(96, shift)
  data[c(1:4, 9), c("a", "b")] <- NA
  data[14, "a"] <- NA
  data
}

test_that("fitted-class descriptives retain unclassified rows as a missing stratum", {
  data <- empty_consumer_data()
  fit <- quietly(multilpa(data, c("a", "b"), "id", 2, n_group_classes = 2,
                         missing = "fiml", n_starts = 1, seed = 2))
  profiles <- descriptives(fit, by = "profile")
  excluded <- profiles[is.na(profiles$profile), ]
  expect_equal(excluded$n, c(0L, 0L))
  expect_equal(excluded$n_missing, c(5L, 5L))
  expect_equal(sum(profiles$n[profiles$variable == "a"] +
                     profiles$n_missing[profiles$variable == "a"]), 96L)
  groups <- descriptives(fit, by = "group_class")
  expect_equal(groups$n_missing[is.na(groups$group_class)], c(4L, 4L))
  expect_false(is.na(.multilpa_assignment_values(fit, "group_class")[9]))
  expect_false(is.na(.multilpa_assignment_values(fit, "profile")[14]))
  fit$drop_empty <- FALSE
  expect_false(anyNA(descriptives(fit, by = "profile")$profile))
  fit$sampling_weights <- c(0, rep(1, 23))
  expect_true(all(is.na(.multilpa_assignment_values(fit, "profile")[1:4])))
})

test_that("sensitivity ignores differences only in excluded classifications", {
  data <- empty_consumer_data()
  fit <- quietly(lpa(data, c("a", "b"), 2, missing = "fiml", n_starts = 1, seed = 2))
  changed <- fit
  changed$subject_profiles[!fit$informative] <- 3L - changed$subject_profiles[!fit$informative]
  expect_equal(.sensitivity_profiles(fit, NULL, NULL, NULL, NULL)$agreement(changed), 1)
  fit$drop_empty <- changed$drop_empty <- FALSE
  expect_equal(.sensitivity_profiles(fit, NULL, NULL, NULL, NULL)$agreement(changed), 91 / 96)
})

test_that("general transition assignments apply the shared classification policy", {
  data <- empty_consumer_data()
  fit <- quietly(lta(data, c("a", "b"), "id", 2, time = "time",
                     transitions = "occasion", missing = "fiml",
                     n_starts = 1, seed = 2, max_iter = 50))
  assigned <- get_results(fit, "assignments")
  fields <- c("profile", "uncertainty", "posterior_profile_1", "posterior_profile_2")
  expect_equal(nrow(assigned), nrow(data))
  expect_identical(assigned$id, data$id)
  expect_true(all(is.na(assigned[c(1:4, 9), fields])))
  expect_true(all(is.na(assigned$group_class[1:4])))
  expect_false(is.na(assigned$group_class[9]))
  expect_false(anyNA(assigned[14, fields]))
  expect_identical(as.data.frame(fit, what = "assignments"), assigned)
  expect_identical(get_results(fit, "all")$assignments, assigned)
  expect_identical(summary(fit)$tables$assignments, assigned)
  expect_false(anyNA(fit$subject_posteriors))
  changed <- fit
  changed$expectation$subject_posteriors[!fit$informative, ] <-
    changed$expectation$subject_posteriors[!fit$informative, 2:1]
  expect_equal(.sensitivity_transitions(fit, data, NULL, NULL, NULL)$agreement(changed), 1)
  fit$drop_empty <- FALSE
  expect_false(anyNA(get_results(fit, "assignments")[fields]))
  fit$sampling_weights <- c(0, rep(1, 23))
  expect_true(all(is.na(get_results(fit, "assignments")[1:4, c(fields, "group_class")])))
})

test_that("sequence summaries and networks use classified groups and occasions", {
  data <- empty_consumer_data()
  fit <- quietly(multilpa(data, c("a", "b"), "id", 2, n_group_classes = 2,
                         time = "time", missing = "fiml", n_starts = 1, seed = 2))
  totals <- get_results(fit, "sequence_summary")
  expect_equal(sum(totals$groups), 23L)
  # Sequence length describes scheduled records within classified groups;
  # the one wholly empty group is excluded, the empty occasion in group 3
  # remains part of that group's recorded span.
  expect_equal(sum(totals$observations), 92L)
  expect_equal(sum(totals$complete), 23L)
  expect_equal(nrow(.gg_case_values(fit)), 91L)
  expect_equal(nrow(.gg_case_values(fit, "groups")), 23L)
  sequences <- .multilpa_modal_sequences(fit)
  expect_equal(nrow(sequences$states), 23L)
  expect_equal(sum(!is.na(as.matrix(sequences$states))), 91L)
  expect_equal(sum(get_results(fit, "classification", level = "individuals")$n_modal), 91L)
  expect_equal(sum(get_results(fit, "classification", level = "groups")$n_modal), 23L)
  fit$drop_empty <- FALSE
  expect_equal(sum(get_results(fit, "sequence_summary")$groups), 24L)
  fit$sampling_weights <- c(0, rep(1, 23))
  expect_equal(sum(get_results(fit, "sequence_summary")$groups), 23L)
})

test_that("homogeneous transition sequence lengths leave empty groups unclassified", {
  data <- empty_consumer_data()
  fit <- quietly(lta(data, c("a", "b"), "id", 2, time = "time",
                     missing = "fiml", n_starts = 1, seed = 2, max_iter = 50))
  lengths <- get_results(fit, "sequence_lengths")
  expect_equal(nrow(lengths), 24L)
  expect_true(is.na(lengths$group_class[lengths$group == 1L]))
  expect_false(is.na(lengths$group_class[lengths$group == 3L]))
  expect_equal(lengths$occasions, rep(4L, 24L))
  fit$drop_empty <- FALSE
  expect_false(anyNA(get_results(fit, "sequence_lengths")$group_class))
})
