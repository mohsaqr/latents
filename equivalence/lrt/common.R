# compare_models(test = "lrt") against the references of make-kit.R, as one
# table: one row per compared quantity. Sourced by compare.R (which writes
# comparison.csv) and by equivalence/test-lrt.R; the package must be loaded.
lrt_comparison <- function(kit) {
  reference <- readRDS(file.path(kit, "reference.rds"))
  row_result <- function(case, quantity, ours, theirs, tolerance) {
    data.frame(case = case, quantity = quantity, latents = ours, reference = theirs,
               difference = abs(ours - theirs), tolerance = tolerance,
               agrees = abs(ours - theirs) <= tolerance)
  }
  set.seed(1)

  # 1. Covariance structures, against mclust.
  structures <- c("EEI", "VVI", "EEE", "VVV")
  fits <- lapply(stats::setNames(structures, structures), function(structure) {
    lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 2, model = structure,
        n_starts = 20, tol = 1e-10, max_iter = 5000, seed = 1)
  })
  chains <- list(c("EEI", "VVI", "VVV"), c("EEI", "EEE", "VVV"))
  mclust_rows <- do.call(rbind, lapply(chains, function(chain) {
    tested <- do.call(compare_models, c(fits[chain], test = "lrt"))
    tested <- as.data.frame(tested)
    do.call(rbind, lapply(2:3, function(i) {
      expected <- subset(reference$mclust$lrt, null == chain[i - 1L] & alternative == chain[i])
      case <- sprintf("mclust %s in %s", chain[i - 1L], chain[i])
      likelihood <- subset(reference$mclust$fits, structure == chain[i])$log_likelihood
      rbind(row_result(case, sprintf("%s log likelihood", chain[i]),
                       fits[[chain[i]]]$log_likelihood, likelihood, 1e-6),
            row_result(case, "LR statistic", tested$lr_statistic[i], expected$statistic, 2e-6),
            row_result(case, "df", tested$lr_df[i], expected$df, 0),
            row_result(case, "p value", tested$lr_p_value[i], expected$p_value,
                       1e-5 * expected$p_value))
    }))
  }))

  # 2. Membership on a covariate, against glca::gofglca().
  gss <- reference$glca$data
  items <- c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE")
  without <- lca(gss, items, n_classes = 3, n_starts = 50, tol = 1e-10, max_iter = 5000,
                 seed = 1)
  with <- lca(gss, items, n_classes = 3, profile_covariates = "AGE", n_starts = 50,
              tol = 1e-10, max_iter = 5000, seed = 1)
  glca_test <- as.data.frame(compare_models(without = without, with = with, test = "lrt"))
  expected <- reference$glca$lrt
  glca_statistic <- 2 * (expected$alternative_log_likelihood - expected$null_log_likelihood)
  glca_rows <- rbind(
    row_result("glca ~ 1 in ~ AGE", "null log likelihood", without$log_likelihood,
               expected$null_log_likelihood, 1e-6),
    row_result("glca ~ 1 in ~ AGE", "alternative log likelihood", with$log_likelihood,
               expected$alternative_log_likelihood, 1e-6),
    row_result("glca ~ 1 in ~ AGE", "LR statistic (glca log likelihoods)",
               glca_test$lr_statistic[2L], glca_statistic, 2e-6),
    row_result("glca ~ 1 in ~ AGE", "LR statistic (gofglca deviance, 3 decimals)",
               glca_test$lr_statistic[2L], expected$statistic, 5e-4),
    row_result("glca ~ 1 in ~ AGE", "df", glca_test$lr_df[2L], expected$df, 0),
    row_result("glca ~ 1 in ~ AGE", "p value (gofglca)", glca_test$lr_p_value[2L],
               expected$p_value, 5e-4))

  # 3. Measurement invariant across occasions against specific to each, Mplus.
  raw <- utils::read.table(file.path(kit, "mplus", "lta.dat"))
  long <- do.call(rbind, lapply(1:2, function(occasion) {
    answers <- raw[, (occasion - 1L) * 3L + 1:3]
    names(answers) <- paste0("u", 1:3)
    cbind(person = seq_len(nrow(raw)), occasion = occasion, answers)
  }))
  long[paste0("u", 1:3)] <- lapply(long[paste0("u", 1:3)], factor)
  transition_fits <- lapply(c(invariant = "invariant", occasion = "occasion"), function(m) {
    lta(long, paste0("u", 1:3), "person", n_profiles = 2, time = "occasion",
        categorical = paste0("u", 1:3), measurement = m, n_starts = 10, seed = 1,
        tol = 1e-10, max_iter = 3000)
  })
  mplus_test <- as.data.frame(compare_models(invariant = transition_fits$invariant,
                                             occasion = transition_fits$occasion,
                                             test = "lrt"))
  expected <- reference$mplus
  mplus_rows <- rbind(
    row_result("Mplus LTA invariant in occasion", "invariant log likelihood",
               transition_fits$invariant$log_likelihood, expected$fits$log_likelihood[1L], 5e-4),
    row_result("Mplus LTA invariant in occasion", "occasion log likelihood",
               transition_fits$occasion$log_likelihood, expected$fits$log_likelihood[2L], 5e-4),
    row_result("Mplus LTA invariant in occasion", "LR statistic",
               mplus_test$lr_statistic[2L], expected$lrt$statistic, 2e-3),
    row_result("Mplus LTA invariant in occasion", "df", mplus_test$lr_df[2L],
               expected$lrt$df, 0),
    row_result("Mplus LTA invariant in occasion", "log10 p value",
               log10(mplus_test$lr_p_value[2L]), log10(expected$lrt$p_value), 1e-2))

  result <- rbind(mclust_rows, glca_rows, mplus_rows)
  row.names(result) <- NULL
  result
}
