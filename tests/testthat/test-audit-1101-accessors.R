skip_on_cran()

# Audit 2026-11-01, slice `accessors`: regression tests for the fixes and the
# independent checks of the public result surface (pooling, predict, the
# likelihood methods and the condition registry).

audit_pool_data <- function() {
  data <- subset(course_engagement, student <= 30)
  withr::with_seed(11, lapply(seq_len(3), function(i) {
    within(data, previous_grade <- previous_grade +
             stats::rnorm(length(previous_grade), sd = 0.3))
  }))
}
audit_pool_vars <- c("browse", "lectures", "forum_read")

test_that("pooled estimates follow Rubin's rules, computed by hand", {
  completed <- audit_pool_data()
  pooled <- quietly(pool_imputations(
    completed, audit_pool_vars, "student", n_profiles = 2, n_group_classes = 1,
    profile_covariates = "previous_grade", n_starts = 2, seed = 1))
  long <- get_results(pooled, "imputations")
  estimates <- get_results(pooled)
  m <- 3L
  key <- factor(paste(long$level, long$outcome, long$term, long$parameter),
                levels = unique(paste(long$level, long$outcome, long$term,
                                      long$parameter)))
  q <- as.vector(tapply(long$estimate, key, mean))
  u <- as.vector(tapply(long$standard_error^2, key, mean))
  b <- as.vector(tapply(long$estimate, key, stats::var))
  total <- u + (1 + 1 / m) * b
  r <- (1 + 1 / m) * b / u
  df <- (m - 1) * (1 + 1 / r)^2
  expect_equal(estimates$estimate, q)
  expect_equal(estimates$within, u)
  expect_equal(estimates$between, b)
  expect_equal(estimates$standard_error, sqrt(total))
  expect_equal(estimates$df, df)
  expect_equal(estimates$fmi, (r + 2 / (df + 3)) / (r + 1))
  expect_equal(estimates$conf_high, q + stats::qt(0.975, df) * sqrt(total))
  tested <- !is.na(estimates$statistic)
  expect_equal(estimates$p_value[tested],
               2 * stats::pt(-abs(q / sqrt(total)), df)[tested])
  # The pooled covariance: mean within-imputation covariance plus the
  # inflated between-imputation covariance, off the diagonal too.
  fits <- lapply(completed, function(data) {
    quietly(multilpa(data, audit_pool_vars, "student", n_profiles = 2,
                     n_group_classes = 1, profile_covariates = "previous_grade",
                     n_starts = 2, seed = 1))
  })
  expect_identical(get_results(pooled, "fits")$profile_order, rep("1,2", m))
  within <- Reduce(`+`, lapply(fits, vcov)) / m
  draws <- do.call(rbind, lapply(fits, coef))[, rownames(within)]
  expect_equal(unname(vcov(pooled)),
               unname(within + (1 + 1 / m) * stats::cov(draws)))
  expect_equal(unname(confint(pooled)[, 2L]), estimates$conf_high)
})

test_that("covariate fits with held measurement pool, with one row per estimated coefficient", {
  completed <- audit_pool_data()
  base <- quietly(multilpa(completed[[1L]], audit_pool_vars, "student",
                           n_profiles = 2, n_group_classes = 1, n_starts = 2,
                           seed = 1))
  start <- starting_values(base, what = "measurement")
  pooled <- quietly(pool_imputations(
    completed, audit_pool_vars, "student", n_profiles = 2, n_group_classes = 1,
    profile_covariates = "previous_grade", fixed = "measurement",
    start = start, n_starts = 2, seed = 1))
  estimates <- get_results(pooled)
  expect_identical(estimates$level, c("profile", "profile"))
  expect_identical(estimates$term, c("(Intercept)", "previous_grade"))
  expect_identical(dim(vcov(pooled)), c(2L, 2L))
  expect_equal(unname(sqrt(diag(vcov(pooled)))), estimates$standard_error)
  expect_identical(names(coef(pooled)), rownames(vcov(pooled)))
})

test_that("get_results() of pooled imputations offers \"all\" and refuses unused arguments", {
  completed <- audit_pool_data()
  pooled <- quietly(pool_imputations(
    completed, audit_pool_vars, "student", n_profiles = 2, n_group_classes = 1,
    n_starts = 2, seed = 1))
  all <- get_results(pooled, "all")
  expect_identical(names(all), c("estimates", "imputations", "fits"))
  expect_identical(all$estimates, get_results(pooled, "estimates"))
  expect_identical(all$fits, get_results(pooled, "fits"))
  expect_identical(get_results(pooled, "imp"), get_results(pooled, "imputations"))
  expect_error(get_results(pooled, "estimates", level = 0.9),
               class = "latents_bad_argument")
  expect_error(get_results(pooled, "nope"), class = "latents_bad_argument")
  expect_error(as.data.frame(pooled, what = "all"), class = "latents_bad_argument")
  expect_identical(as.data.frame(pooled), get_results(pooled))
  expect_warning(covariance <- vcov(pooled, type = "robust"),
                 class = "latents_unused_argument")
  expect_identical(covariance, vcov(pooled))
})

test_that("predict() of a profile model refuses a bad type and warns about a misnamed argument", {
  fit <- quietly(lpa(srl, names(srl), n_profiles = 2, n_starts = 2, seed = 1))
  expected <- predict(fit)
  expect_warning(misnamed <- predict(fit, data = srl), class = "latents_unused_argument")
  expect_identical(misnamed, expected)
  expect_error(predict(fit, type = "bogus"), class = "latents_bad_argument")
  expect_identical(predict(fit, type = "post"), predict(fit, type = "posterior"))
  # The fitted rows' predictions reproduce the stored posteriors.
  posterior <- predict(fit, type = "posterior")
  wide <- get_results(fit, "posteriors", format = "wide")
  expect_equal(unname(as.matrix(posterior[c("probability_profile_1", "probability_profile_2")])),
               unname(as.matrix(wide[c("posterior_profile_1", "posterior_profile_2")])))
  expect_equal(sum(predict(fit, type = "density")$log_density), fit$log_likelihood)
})

test_that("logLik(), nobs(), AIC() and BIC() agree with the stored fit", {
  data <- subset(course_engagement, student <= 30)
  fit <- quietly(multilpa(data, audit_pool_vars, "student", n_profiles = 2,
                          n_group_classes = 2, n_starts = 2, seed = 1))
  likelihood <- logLik(fit)
  expect_equal(as.numeric(likelihood), fit$log_likelihood)
  expect_identical(attr(likelihood, "df"), fit$n_parameters)
  expect_identical(nobs(fit), 30L)
  expect_equal(AIC(fit), fit$aic)
  expect_equal(BIC(fit), fit$bic)
  expect_equal(BIC(fit), -2 * fit$log_likelihood + fit$n_parameters * log(30))
})

test_that("vcov() and confint() cover the estimated parameters only; coef() every one", {
  data <- subset(course_engagement, student <= 30)
  base <- quietly(multilpa(data, audit_pool_vars, "student", n_profiles = 2,
                           n_group_classes = 1, n_starts = 2, seed = 1))
  held <- quietly(multilpa(data, audit_pool_vars, "student", n_profiles = 2,
                           n_group_classes = 1, profile_covariates = "previous_grade",
                           fixed = "measurement",
                           start = starting_values(base, what = "measurement"),
                           n_starts = 2, seed = 1))
  expect_length(coef(held), 14L)
  expect_identical(rownames(vcov(held)),
                   c("profile.coefficient.profile_1.(Intercept)",
                     "profile.coefficient.profile_1.previous_grade"))
  expect_identical(rownames(confint(held)), rownames(vcov(held)))
  expect_error(confint(held, parm = names(coef(held))[1L]),
               class = "latents_held_parameter")
})

test_that("units of sampling weight zero keep their assignment, as documented", {
  data <- withr::with_seed(3, data.frame(
    y1 = c(stats::rnorm(40), stats::rnorm(40, 3)),
    y2 = c(stats::rnorm(40), stats::rnorm(40, 3)), w = 1))
  data$w[c(5L, 60L)] <- 0
  fit <- quietly(lpa(data, c("y1", "y2"), n_profiles = 2, weights = "w",
                     n_starts = 2, seed = 1))
  assignments <- get_results(fit, "assignments")
  expect_false(anyNA(assignments$profile[c(5L, 60L)]))
  expect_identical(predict(fit)$profile, assignments$profile)
  expect_identical(nobs(fit), 78L)
  expect_equal(sum(get_results(fit, "classification")$n_modal), 78)
})

test_that("the replication table has the columns ?get_results lists", {
  fit <- quietly(lpa(srl, names(srl), n_profiles = 2, n_starts = 2, seed = 1))
  expect_identical(names(get_results(fit, "replication")),
                   c("n_starts", "n_screened", "n_converged", "n_failed",
                     "n_best_replicated", "share", "tolerance", "best_start",
                     "log_likelihood"))
})

test_that("every latents_ condition class raised in R/ is documented, and every documented one is raised", {
  # test-tidy-core.R's audit of the same contract still matches the
  # `multilpa_` prefix the classes carried before the rename, so it compares
  # two empty sets. This one reads the `latents_` classes.
  source_dir <- test_path("..", "..", "R")
  skip_if_not(file.exists(file.path(source_dir, "conditions.R")),
              "R/ sources are not available")
  pattern <- paste0('class\\s*=\\s*(c\\()?\\s*"latents_[A-Za-z0-9_]*"',
                    '(\\s*,\\s*"latents_[A-Za-z0-9_]*")*')
  raised <- unlist(lapply(list.files(source_dir, pattern = "[.]R$", full.names = TRUE),
                          function(path) {
    lines <- readLines(path, warn = FALSE)
    text <- paste(ifelse(grepl("^\\s*#", lines), "", lines), collapse = "\n")
    at <- gregexpr(pattern, text)[[1L]]
    if (at[1L] == -1L) return(character())
    found <- regmatches(text, list(at))[[1L]]
    last_at <- function(what, before) {
      hits <- gregexpr(what, substr(text, 1L, before), perl = TRUE)[[1L]]
      if (hits[1L] == -1L) 0L else max(hits)
    }
    # A `class =` belongs to the constructor opened most recently before it;
    # `structure(x, class = ...)` builds an object, not a condition, unless
    # it is the condition passed to message(); `.multilpa_notice()` raises
    # classed messages; `.latents_lrt_structure(` is not `structure(`.
    is_condition <- vapply(at, function(start) {
      last_at(paste0("(error|warning|message)Condition\\(|message\\(structure\\(|",
                     "\\.multilpa_notice\\("), start) >=
        last_at("(?<![A-Za-z0-9_.])structure\\(|class\\s*<-", start) - nchar("message(")
    }, logical(1))
    unlist(regmatches(found[is_condition],
                      gregexpr("latents_[A-Za-z0-9_]+", found[is_condition])))
  }))
  doc <- paste(readLines(file.path(source_dir, "conditions.R"), warn = FALSE),
               collapse = "\n")
  items <- regmatches(doc, gregexpr("\\\\item\\{[^}]*\\}", doc))[[1L]]
  documented <- unique(gsub("`", "", unlist(regmatches(
    items, gregexpr("`latents_[A-Za-z0-9_]+`", items)))))
  # Documented as "Not raised": a half-observed pair row is summed over.
  documented <- setdiff(documented, "latents_partial_pair")
  expect_gt(length(documented), 90L)
  expect_setequal(unique(raised), documented)
})
