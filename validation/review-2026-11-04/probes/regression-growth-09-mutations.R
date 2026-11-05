# In-memory mutations of table formulas; test files run from a scratchpad copy, so no snapshot is touched.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
ns <- asNamespace("latents")
wrap <- function(name, post) { old <- get(name, envir = ns); new <- function(...) post(old(...)); environment(new) <- environment(); unlockBinding(name, ns); assign(name, new, envir = ns); lockBinding(name, ns) }
wrap(".mixture_coefficient_table", function(t) { t$p_adjusted <- p.adjust(t$p_value, "BH"); if (!is.null(t$exp_estimate)) t$exp_estimate <- t$estimate; t })
wrap(".mixture_fit_table", function(t) { t$sabic <- t$sabic + 1; t$bic_rows <- t$bic; t })
wrap(".mixture_recovery_table", function(t) { t$share <- t$n / sum(t$n); t })
wrap(".mixture_group_class_table", function(t) { if (nrow(t)) { t$group_share <- 1 - t$group_share; if ("std_error" %in% names(t)) t$std_error <- t$std_error * 2 }; t })
wrap(".mixture_class_table", function(t) { if (!is.null(t$sigma_std_error)) t$sigma_std_error <- t$sigma_std_error * 2; if (!is.null(t$dispersion_std_error)) t$dispersion_std_error <- t$dispersion_std_error * 2; if (!is.null(t$prior_std_error)) t$prior_std_error <- t$prior_std_error * 2; t })
wrap(".growth_random_table", function(t) { t$conf_low <- t$estimate - 2 * (t$estimate - t$conf_low); t$std_error <- t$std_error * 2; t })
wrap(".growth_class_table", function(t) { t$share_conf_low <- t$share_conf_low / 2; t$share_std_error <- t$share_std_error * 2; t$residual_sd_std_error <- t$residual_sd_std_error * 2; t })
wrap(".growth_cluster_fit_table", function(t) { t$bic_persons <- t$bic; t$sabic <- t$sabic + 1; t })
wrap(".growth_fit_table", function(t) { t$bic_rows <- t$bic; t$sabic <- t$sabic + 1; t })
files <- c("test-mixture-regression.R", "test-mixture-regression-labels.R", "test-release-regression.R", "test-review-regression.R", "test-review-1026-mixreg.R", "test-review-1026-mixreg-r2.R", "test-audit-1101-mixture-regression.R", "test-growth-mixture.R", "test-multilevel-growth.R", "test-audit-1101-growth-mixture.R", "test-review-1026-growth-r2.R")
res <- do.call(rbind, lapply(files, function(f) {
  out <- as.data.frame(testthat::test_file(file.path("rg-mut", "testthat", f), reporter = "silent", load_package = "none", stop_on_failure = FALSE))
  bad <- subset(out, failed > 0 | error)
  if (nrow(bad)) { cat("FAILED in", f, ":\n"); print(bad[c("test", "failed", "error")]) }
  data.frame(file = f, tests = nrow(out), expectations = sum(out$nb), failed = sum(out$failed), errors = sum(out$error), secs = round(sum(out$real), 1))
}))
print(res)
