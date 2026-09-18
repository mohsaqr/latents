# Verify the Lo-Mendell-Rubin adjustment factor against two retained genuine
# Mplus 9 TECH11 outputs, and record that the VLMR reference distribution is
# deliberately not reproduced.
# Run from the project root: Rscript validation/mplus/lmr/compare.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
artifact_dir <- file.path("validation", "mplus", "lmr")

# Mplus prints TECH11 as a fixed block; read the labelled numbers out of it.
read_tech11 <- function(prefix) {
  output <- readLines(file.path(artifact_dir, paste0(prefix, ".out")))
  stopifnot("Mplus did not terminate normally" =
              any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", output)),
            "TECH11 was not produced" =
              any(grepl("VUONG-LO-MENDELL-RUBIN", output)))
  trailing_number <- function(line) {
    trimmed <- trimws(line)
    as.numeric(regmatches(trimmed, regexpr("-?[0-9]+\\.?[0-9]*$", trimmed)))
  }
  value <- function(pattern) {
    trailing_number(grep(pattern, output, value = TRUE)[1L])
  }
  list(statistic = value("2 Times the Loglikelihood Difference"),
       df = value("Difference in the Number of Parameters"),
       adjusted = trailing_number(
         output[grep("LO-MENDELL-RUBIN ADJUSTED", output)[1L] + 2L]),
       vlmr_p = value("^ +P-Value"))
}

references <- lapply(c(lmr2 = "lmr2", lmr3 = "lmr3"), function(prefix) {
  raw <- read_tech11(prefix)
  data.frame(comparison = prefix, statistic = raw$statistic, df = raw$df,
             mplus_adjusted = raw$adjusted)
})
reference <- do.call(rbind, references)
reference$n <- 500L

# The native adjustment is a pure function of the statistic, df and n, so it is
# applied here directly to the Mplus statistics rather than refitting in R.
reference$native_adjusted <- reference$statistic /
  (1 + 1 / (reference$df * log(reference$n)))
reference$absolute_difference <- abs(reference$native_adjusted - reference$mplus_adjusted)
print(reference, digits = 10, row.names = FALSE)

# Mplus prints TECH11 to three decimals, so the statistic it reports is itself
# rounded; the tolerance below is that printing precision, not an estimate.
stopifnot("LMR adjustment disagrees with Mplus" =
            max(reference$absolute_difference) < 5e-3)

# End-to-end: fit both models natively on the same published data and compare
# the statistic Mplus computed from its own one- and two-class fits. With one
# group class the likelihood is the pooled single-level LPA likelihood, so the
# cluster column carries no information here and only supplies the group layout.
fixture <- readRDS(file.path("tests", "fixtures", "mplus", "public-7.9.rds"))
indicators <- c("y1", "y2", "y3", "y4")
# Example 7.9 is the shared-variance specification, which is also the Mplus
# TYPE = MIXTURE default; the TECH11 run above inherited it.
stopifnot("fixture 7.9 is not the equal-variance example" =
            identical(fixture$variance_model, "equal"))
fitted_null <- fit_multilpa(fixture$data, indicators, "cluster", 1L, 1L,
                          variance_model = "equal", n_starts = 5L,
                          max_iter = 5000L, tol = 1e-12, seed = 912)
fitted_alternative <- fit_multilpa(fixture$data, indicators, "cluster", 2L, 1L,
                                 variance_model = "equal", n_starts = 30L,
                                 max_iter = 5000L, tol = 1e-12, seed = 912)
stopifnot("native parameter difference does not match Mplus" =
            fitted_alternative$n_parameters - fitted_null$n_parameters ==
            reference$df[1L])
stopifnot("native fits did not converge" =
            fitted_null$converged && fitted_alternative$converged)
native <- lmr_lrt_multilpa(fitted_null, fitted_alternative)
print(native, digits = 10, row.names = FALSE)
stopifnot(
  "native statistic disagrees with Mplus" =
    abs(native$statistic - reference$statistic[1L]) < 5e-3,
  "native adjusted statistic disagrees with Mplus" =
    abs(native$adjusted_statistic - reference$mplus_adjusted[1L]) < 5e-3,
  "a p-value must not be reported" = is.na(native$p_value))

saveRDS(reference, file.path("tests", "fixtures", "mplus", "lmr-tech11.rds"))
write.csv(reference, file.path(artifact_dir, "comparison.csv"), row.names = FALSE)
cat("\nVLMR reference distribution (Mplus Mean/Standard Deviation and p-value)\n",
    "is deliberately not reproduced; see R/lmr.R and the README.\n", sep = "")
cat("LMR comparison complete.\n")
