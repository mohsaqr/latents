# External check of the weight convention against Mplus 9 (TYPE = COMPLEX,
# WEIGHT IS wt; y ON x). A one-profile, full-covariance weighted fit of (y, x)
# estimates the weighted mean and covariance of the pair; the regression of y
# on x and its conditional log likelihood follow from them exactly. Mplus
# rescales the weights to sum to n, as latents does, so the log likelihoods
# agree only if that scaling is the same.
# Run from the package root: Rscript equivalence/weights/compare-mplus-complex.R
devtools::load_all(".", quiet = TRUE)
reference_dir <- file.path("..", "JStats", "validation", "r-reference")
data <- utils::read.table(file.path(reference_dir, "mplus-demo-complex-weight-on.dat"),
                          col.names = c("y", "x", "wt", "clus", "strat", "sub"))
mplus <- c(slope = 0.892, intercept = -0.108, residual_variance = 3.779,
           log_likelihood = -4167.500)
regression <- function(fit) {
  covariance <- fit$covariances[, , 1L]
  means <- fit$means[1L, ]
  slope <- covariance["y", "x"] / covariance["x", "x"]
  # Conditional log likelihood = joint minus the marginal of x, which at the
  # weighted maximum is -n/2 (log(2 pi var_x) + 1).
  marginal_x <- -nrow(data) / 2 * (log(2 * pi * covariance["x", "x"]) + 1)
  c(slope = slope, intercept = means[["y"]] - slope * means[["x"]],
    residual_variance = covariance["y", "y"] - slope^2 * covariance["x", "x"],
    log_likelihood = fit$log_likelihood - marginal_x)
}
fit <- withCallingHandlers(
  lpa(data, c("y", "x"), 1, model = "VVV", weights = "wt", n_starts = 1, seed = 1),
  latents_single_level = function(notice) invokeRestart("muffleMessage"))
latents <- regression(fit)
# An independent closed form: weighted least squares with the weights scaled
# to sum to n is the weighted maximum likelihood regression.
scaled <- data$wt * nrow(data) / sum(data$wt)
wls <- stats::lm(y ~ x, data = data, weights = scaled)
residual_variance <- sum(scaled * stats::residuals(wls)^2) / nrow(data)
closed_form <- c(slope = unname(stats::coef(wls)[["x"]]),
                 intercept = unname(stats::coef(wls)[["(Intercept)"]]),
                 residual_variance = residual_variance,
                 log_likelihood = -nrow(data) / 2 *
                   (log(2 * pi * residual_variance) + 1))
# Mplus prints three decimals. Its residual variance (3.779) sits 6e-4 below
# the closed-form maximum (3.7796): the likelihood is flat there, so its
# printed log likelihood still matches to 1e-4. The Mplus gate is therefore
# 1e-3; the closed form is the tight one.
comparison <- data.frame(quantity = names(mplus), mplus = unname(mplus),
                         closed_form = unname(closed_form),
                         latents = unname(latents),
                         versus_closed_form = unname(latents - closed_form),
                         versus_mplus = unname(latents - mplus))
comparison$agrees <- abs(comparison$versus_closed_form) <= 1e-8 &
  abs(comparison$versus_mplus) <= 1e-3
print(comparison, row.names = FALSE, digits = 7)
stopifnot("latents must reproduce the weighted estimates" = all(comparison$agrees))
