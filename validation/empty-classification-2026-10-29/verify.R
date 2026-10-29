# Run from the package root with one BLAS thread.
pkgload::load_all(quiet = TRUE)
set.seed(27)
profile <- c(rep(1L, 90L), rep(2L, 210L))
data <- data.frame(a = rnorm(300L, ifelse(profile == 1L, -3, 3)),
                   b = rnorm(300L, ifelse(profile == 1L, -3, 3)))
data[9L, ] <- NA_real_
fit <- lpa(data, c("a", "b"), n_profiles = 2L, missing = "fiml",
           n_starts = 3L, seed = 1L)
assigned <- get_results(fit, "assignments")
print(assigned[8:10, ])
print(c(input_rows = nrow(data), output_rows = nrow(assigned),
        classified = sum(!is.na(assigned$profile)), nobs = nobs(fit),
        classification_total = sum(get_results(fit, "classification")$n_modal)))
stopifnot(nrow(assigned) == 300L, sum(!is.na(assigned$profile)) == 299L,
          all(is.na(assigned[9L, ])), !anyNA(fit$subject_posteriors),
          isTRUE(all.equal(unname(fit$subject_posteriors[9L, ]),
                            as.vector(fit$profile_probabilities))))

# The saved pre-fix reproduction is optional; when available, verify that
# reporting corrections did not alter the original Gaussian fit.
before_path <- "tmp/empty-assignment-review/reproduction.rds"
if (file.exists(before_path)) {
  before <- readRDS(before_path)$fit
  fields <- c("log_likelihood", "means", "variances", "subject_posteriors",
              "group_posteriors", "profile_probabilities")
  stopifnot(all(vapply(fields, function(name) {
    isTRUE(all.equal(fit[[name]], before[[name]], tolerance = 1e-12))
  }, logical(1))))
  cat("Original fit parameters, likelihood and internal posteriors preserved.\n")
}
write.csv(assigned, "validation/empty-classification-2026-10-29/assignments.csv",
          row.names = FALSE, na = "NA")
sources <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
             list.files("tests/testthat", pattern = "[.]R$", full.names = TRUE))
write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
          "validation/empty-classification-2026-10-29/source-manifest.csv",
          row.names = FALSE)
cat("Synthetic reproduction and source manifest saved.\n")
