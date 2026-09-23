# Compare the native two-level latent class model against a retained genuine
# Mplus 9 run with CATEGORICAL binary indicators.
# Run from the project root: Rscript equivalence/mplus/categorical/compare.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("equivalence", "fixtures.R"))
artifact_dir <- file.path("equivalence", "mplus", "categorical")

indicators <- paste0("u", 1:5)
raw <- read.table(file.path(artifact_dir, "categorical.dat"),
                  col.names = c(indicators, "clus", "id"))
str(raw)
print(head(raw))
print(vapply(raw, class, character(1)))
print(summary(raw))
stopifnot(
  "unexpected data shape" = nrow(raw) == 1200L && length(unique(raw$clus)) == 60L,
  "indicators must be binary 0/1" =
    all(vapply(raw[, indicators], function(value)
      identical(sort(unique(value)), c(0L, 1L)), logical(1))),
  "identifiers must be unique" = !anyDuplicated(raw$id))

output <- readLines(file.path(artifact_dir, "categorical.out"))
stopifnot(
  "Mplus did not terminate normally" =
    any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
  "Mplus did not replicate the best loglikelihood" =
    any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", output)),
  "Mplus reported an error" =
    !any(grepl("*** ERROR", output, fixed = TRUE)),
  "the run did not treat the indicators as categorical" =
    any(grepl("TAU\\(U\\) FOR LATENT CLASS PATTERN", output)))

n_parameters <- 13L
values <- scan(file.path(artifact_dir, "categorical-results.dat"), quiet = TRUE)
stopifnot("unexpected saved-results layout" =
            length(values) == 2L * n_parameters + 7L &&
            values[2L * n_parameters + 1L] == n_parameters)
# TECH1 fixes the layout: thresholds 1-5 belong to cw#1 and 6-10 to cw#2, shared
# across the cb patterns, then the cb logit, the cw#1 intercept, and the cw#1 ON
# cb slope. The MODEL cw block is what imposes measurement invariance across cb,
# matching this package's model; without it Mplus frees 20 thresholds, not 10.
mplus_thresholds <- rbind(values[1:5], values[6:10])
group_logit <- values[n_parameters - 2L]
intercept <- values[n_parameters - 1L]
slope <- values[n_parameters]
first_profile <- stats::plogis(c(intercept + slope, intercept))
mplus_profile_probabilities <- cbind(first_profile, 1 - first_profile)
mplus_group_probabilities <- c(stats::plogis(group_logit),
                               1 - stats::plogis(group_logit))
mplus_log_likelihood <- values[2L * n_parameters + 2L]

fit <- multilpa(raw, indicators, "clus", n_profiles = 2, n_group_classes = 2,
                  categorical = indicators, n_starts = 40, seed = 20260918,
                  tol = 1e-13, max_iter = 20000)
stopifnot("the native fit did not converge" = fit$converged,
          "the native fit did not replicate its best likelihood" =
            fit$n_best_replicated >= 2L,
          "parameter counts disagree" = fit$n_parameters == n_parameters)

# A threshold is qlogis(P(y = lowest category)), which is exactly the Mplus
# u$1 parameterization, so the two are directly comparable once labels align.
native_thresholds <- vapply(fit$response_probabilities, function(block) {
  as.vector(.multilpa_categorical_thresholds(block))
}, numeric(fit$n_profiles))
profile_order <- order(native_thresholds[, 1L])[rank(mplus_thresholds[, 1L])]
group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
  rank(mplus_profile_probabilities[, 1L])]

differences <- c(
  log_likelihood = abs(fit$log_likelihood - mplus_log_likelihood),
  thresholds = max(abs(native_thresholds[profile_order, ] - mplus_thresholds)),
  profile_probabilities = max(abs(
    fit$profile_probabilities[group_order, profile_order] -
      mplus_profile_probabilities)),
  group_probabilities = max(abs(
    fit$group_probabilities[group_order] - mplus_group_probabilities)))
print(differences, digits = 12)

comparison <- data.frame(
  quantity = c(sprintf("threshold_profile%d_%s",
                       rep(seq_len(2L), times = 5L),
                       rep(indicators, each = 2L)),
               "profile_probability_class1", "profile_probability_class2",
               "group_probability_class1", "log_likelihood"),
  native = c(as.vector(t(native_thresholds[profile_order, ])),
             fit$profile_probabilities[group_order, profile_order][, 1L],
             fit$group_probabilities[group_order][1L], fit$log_likelihood),
  mplus = c(as.vector(t(mplus_thresholds)),
            mplus_profile_probabilities[, 1L],
            mplus_group_probabilities[1L], mplus_log_likelihood))
comparison$absolute_difference <- abs(comparison$native - comparison$mplus)
print(comparison, digits = 10, row.names = FALSE)

# Mplus saves eight significant digits and prints the likelihood to three
# decimals; the tolerances below are that precision, not an estimate.
stopifnot(
  "thresholds disagree with Mplus" = differences[["thresholds"]] < 1e-4,
  "profile probabilities disagree with Mplus" =
    differences[["profile_probabilities"]] < 1e-5,
  "group probabilities disagree with Mplus" =
    differences[["group_probabilities"]] < 1e-5,
  "log likelihood disagrees with Mplus" =
    differences[["log_likelihood"]] < 1e-3)

saveRDS(list(data = raw, mplus_thresholds = mplus_thresholds,
             mplus_profile_probabilities = mplus_profile_probabilities,
             mplus_group_probabilities = mplus_group_probabilities,
             mplus_log_likelihood = mplus_log_likelihood,
             n_parameters = n_parameters),
        mplus_fixture("twolevel-categorical.rds"))
write.csv(comparison, file.path(artifact_dir, "comparison.csv"), row.names = FALSE)
cat("Categorical comparison complete.\n")
