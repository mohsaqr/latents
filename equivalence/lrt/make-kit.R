# Reference values for the likelihood-ratio test of compare_models(test = "lrt").
# Run from the project root: Rscript equivalence/lrt/make-kit.R
# Writes equivalence/lrt/reference.rds. The Mplus outputs in mplus/ are made by
# running /Applications/MplusDemo/mpdemo on mplus/*.inp (see README.md).
set.seed(1)
suppressPackageStartupMessages(library(mclust))  # Mclust() calls mclustBIC() unqualified
kit <- file.path("equivalence", "lrt")

# 1. mclust: covariance structures nested letter by letter (volume, shape,
#    orientation), two profiles on the Old Faithful eruptions.
faithful_data <- datasets::faithful
structures <- c("EEI", "VVI", "EEE", "VVV")
mclust_fits <- lapply(stats::setNames(structures, structures), function(structure) {
  # mclust's default EM tolerance (1e-5, relative) leaves the likelihood a few
  # 1e-4 short of its maximum; the test needs it to 1e-8.
  fit <- mclust::Mclust(faithful_data, G = 2, modelNames = structure, verbose = FALSE,
                        control = mclust::emControl(tol = c(1e-12, 1e-12),
                                                    itmax = c(1e6, 1e6)))
  data.frame(structure = structure, log_likelihood = fit$loglik, n_parameters = fit$df)
})
mclust_table <- do.call(rbind, mclust_fits)
row.names(mclust_table) <- NULL
pairs <- data.frame(null = c("EEI", "VVI", "EEI", "EEE"),
                    alternative = c("VVI", "VVV", "EEE", "VVV"))
mclust_lrt <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  null <- mclust_fits[[pairs$null[i]]]
  alternative <- mclust_fits[[pairs$alternative[i]]]
  statistic <- 2 * (alternative$log_likelihood - null$log_likelihood)
  df <- alternative$n_parameters - null$n_parameters
  data.frame(null = pairs$null[i], alternative = pairs$alternative[i],
             statistic = statistic, df = df,
             p_value = stats::pchisq(statistic, df, lower.tail = FALSE))
}))

# 2. glca: latent class membership on a covariate, tested by gofglca()'s
#    deviance table. Complete cases, so both programs see the same rows.
gss <- glca::gss08
items <- c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE")
gss <- gss[stats::complete.cases(gss[c(items, "AGE")]), c(items, "AGE")]
gss$AGE <- (gss$AGE - mean(gss$AGE)) / stats::sd(gss$AGE)
glca_without <- glca::glca(glca::item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ 1,
                           data = gss, nclass = 3, n.init = 50, seed = 1,
                           verbose = FALSE)
glca_with <- glca::glca(glca::item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ AGE,
                        data = gss, nclass = 3, n.init = 50, seed = 1,
                        verbose = FALSE)
deviance <- glca::gofglca(glca_without, glca_with, test = "chisq")$dtable
glca_lrt <- data.frame(null = "~ 1", alternative = "~ AGE",
                       null_log_likelihood = deviance[1L, "loglik"],
                       alternative_log_likelihood = deviance[2L, "loglik"],
                       null_parameters = deviance[1L, "npar"],
                       alternative_parameters = deviance[2L, "npar"],
                       statistic = deviance[2L, "Deviance"], df = deviance[2L, "Df"],
                       p_value = deviance[2L, "Pr(>Chi)"])

# 3. Mplus 9 demo: two-occasion LTA, measurement invariant against specific to
#    each occasion. Mplus prints the likelihood to three decimals.
read_mplus <- function(name) {
  lines <- readLines(file.path(kit, "mplus", paste0(name, ".out")))
  stopifnot("Mplus did not terminate normally" =
              any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
            "Mplus did not replicate its best likelihood" =
              any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
  number <- function(pattern) {
    as.numeric(sub(".*\\s(-?[0-9.]+)\\s*$", "\\1", grep(pattern, lines, value = TRUE)[1L]))
  }
  data.frame(model = name, log_likelihood = number("H0 Value"),
             n_parameters = number("^Number of Free Parameters"))
}
mplus_table <- rbind(read_mplus("lta-invariant"), read_mplus("lta-occasion"))
mplus_statistic <- 2 * (mplus_table$log_likelihood[2L] - mplus_table$log_likelihood[1L])
mplus_df <- mplus_table$n_parameters[2L] - mplus_table$n_parameters[1L]
mplus_lrt <- data.frame(null = "invariant", alternative = "occasion",
                        statistic = mplus_statistic, df = mplus_df,
                        p_value = stats::pchisq(mplus_statistic, mplus_df,
                                                lower.tail = FALSE))

reference <- list(mclust = list(fits = mclust_table, lrt = mclust_lrt),
                  glca = list(data = gss, lrt = glca_lrt),
                  mplus = list(fits = mplus_table, lrt = mplus_lrt),
                  versions = c(R = as.character(getRversion()),
                               mclust = as.character(utils::packageVersion("mclust")),
                               glca = as.character(utils::packageVersion("glca")),
                               mplus = "Mplus 9 demo (mpdemo)"))
saveRDS(reference, file.path(kit, "reference.rds"))
print(mclust_table, digits = 10); print(mclust_lrt, digits = 10)
print(glca_lrt, digits = 10); print(mplus_table, digits = 10); print(mplus_lrt, digits = 10)
