# Within-person z-scores against Mplus (mplus/person-scaling/). Mplus has no
# within-person standardization, so it is given the data standardized within
# each group beforehand (sample standard deviation, denominator n - 1), and
# fits the two-level VVI model with 2 profiles and 2 group classes. The
# package, given the raw data and `centering = "person", scaling = "person"`,
# reports the log likelihood of the scaled values, so it must reproduce
# Mplus's fit to them: likelihood, parameter count and posteriors.

test_that("person centring and scaling reproduce Mplus on within-person z-scores", {
  directory <- file.path("mplus", "person-scaling")
  lines <- read_mplus_lines(file.path(directory, "scaled.out"))
  expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
  expect_true(any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(directory, "scaled-results.dat"), quiet = TRUE)
  vars <- paste0("y", 1:3)
  raw <- utils::read.table(file.path(directory, "raw.dat"), col.names = c("id", "g", vars))
  fit <- suppressWarnings(multilpa(raw, vars, "g", n_profiles = 2, n_group_classes = 2,
                                   model = "VVI", centering = "person", scaling = "person",
                                   n_starts = 30, seed = 1, tol = 1e-10, max_iter = 5000))
  expect_identical(as.integer(fit$n_parameters), npar)
  expect_lt(abs(fit$log_likelihood - values[2L * npar + 2L]), 5e-4)
  saved <- utils::read.table(file.path(directory, "scaled-posteriors.dat"), na.strings = "*")
  columns <- montecarlo_save_columns(lines)
  names(saved) <- columns
  joint <- as.matrix(saved[grep("^CPROB", columns)])
  profiles <- cbind(joint[, 1] + joint[, 3], joint[, 2] + joint[, 4])
  rows <- match(saved$ID, raw$id)
  expect_lt(montecarlo_relabel(fit$subject_posteriors[rows, , drop = FALSE], profiles)$gap,
            1e-3)
})
