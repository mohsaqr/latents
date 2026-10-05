# The numbers stated in the prose of the workflow-lpa and workflow-lca articles
# (vignettes/articles/). The articles' own code is extracted and run, so
# the objects checked here are the ones the text describes, produced by the
# same calls in the same order under the same seed. A change to the
# estimation that moves one of the stated numbers fails here.

skip_on_cran()
skip_if_not_installed("knitr")

.run_vignette <- function(name) {
  path <- test_path("..", "..", "vignettes", "articles", paste0(name, ".Rmd"))
  skip_if_not(file.exists(path), "the vignette sources are not available")
  script <- tempfile(fileext = ".R")
  knitr::purl(path, output = script, quiet = TRUE, documentation = 0L)
  environment <- new.env(parent = globalenv())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(sys.source(script, envir = environment))
  environment
}

test_that("the numbers in the latent profile workflow hold", {
  run <- .run_vignette("workflow-lpa")
  candidates <- get_results(run$models, "candidates")
  eee <- candidates[candidates$model == "EEE", ]
  # "The EEE structure attains the lowest values"
  expect_identical(candidates$model[which.min(candidates$bic_individual)], "EEE")
  # "BIC is lowest for one profile, and the two-profile model is 4.7 points higher"
  expect_identical(which.min(eee$bic_individual), 1L)
  expect_equal(round(eee$bic_individual[2] - eee$bic_individual[1], 1), 4.7)
  # "their BIC continues to fall up to six profiles" (no covariances)
  falling <- vapply(c("EEI", "VVI"), function(structure) {
    all(diff(candidates$bic_individual[candidates$model == structure]) < 0)
  }, logical(1))
  expect_true(all(falling))
  # "The test rejects the one-profile model (p = 0.005)"
  expect_equal(as.data.frame(run$lrt)$p_value, 0.005)
  # "69% ... higher means on the first four ... 31% ... differ little in
  # test anxiety"
  fit <- run$fit
  shares <- get_results(fit, "profile_probabilities")$probability
  high <- which.max(shares)
  expect_equal(round(sort(shares), 2), c(0.31, 0.69))
  profiles <- get_results(fit, "profiles")
  high_means <- profiles$mean[profiles$profile == high]
  low_means <- profiles$mean[profiles$profile != high]
  expect_true(all(high_means[1:4] > low_means[1:4]))
  expect_lt(abs(high_means[5] - low_means[5]), 0.3)
  # The table the text shows carries its standard errors.
  expect_false(anyNA(profiles$mean_standard_error))
  expect_false(anyNA(profiles$variance_standard_error))
  # "relative entropy of 0.74 and average posterior probabilities of at least 0.88"
  expect_equal(round(get_results(fit, "entropy")$relative_entropy[1], 2), 0.74)
  averages <- get_results(fit, "average_posteriors")
  diagonal <- averages$average_posterior[averages$assigned_class == averages$class &
                                           averages$level == "individuals"]
  expect_gte(min(diagonal), 0.88)
})

test_that("the numbers in the latent class workflow hold", {
  run <- .run_vignette("workflow-lca")
  candidates <- get_results(run$models, "candidates")
  bic <- candidates$bic_individual
  # "BIC is lowest for four classes, 13.2 points below the three-class model"
  expect_identical(which.min(bic), 4L)
  expect_equal(round(bic[3] - bic[4], 1), 13.2)
  # "ICL is lower for three classes than for four ... lowest for one class ...
  # among more than one class, lowest for two ... relative entropy below 0.5"
  icl <- candidates$icl_individual
  expect_lt(icl[3], icl[4])
  expect_identical(which.min(icl), 1L)
  expect_identical(which.min(icl[-1L]) + 1L, 2L)
  expect_lt(candidates$profile_entropy[2], 0.5)
  # "The fourth class contains 1.3% of the reports"
  four <- candidate_fit(run$models, n_profiles = 4, n_group_classes = 1)
  expect_equal(round(min(get_results(four, "profile_probabilities")$probability), 3),
               0.013)

  fit <- run$fit
  responses <- get_results(fit, "responses")
  yes <- responses[responses$category == "yes", ]
  job <- yes[yes$indicator == "part_time_job", ]
  social <- yes[yes$indicator == "on_social_media", ]
  # Classes are identified by content, not by their number.
  job_class <- job$profile[which.max(job$probability)]
  screen_class <- social$profile[which.max(social$probability)]
  away_class <- setdiff(1:3, c(job_class, screen_class))
  expect_length(away_class, 1L)
  expect_gt(max(job$probability), 0.9)
  expect_lt(social$probability[social$profile == away_class], 0.05)
  # "58% ... 4% ... 39%"
  shares <- get_results(fit, "profile_probabilities")$probability
  expect_equal(round(shares[c(screen_class, job_class, away_class)], 2),
               c(0.58, 0.04, 0.39))
  # "relative entropy of 0.66"
  expect_equal(round(get_results(fit, "entropy")$relative_entropy[1], 2), 0.66)
  # "5.30 ... 4.63 ... 4.53"
  expect_equal(round(run$happiness$estimate[c(away_class, screen_class, job_class)], 2),
               c(5.30, 4.63, 4.53))
})
