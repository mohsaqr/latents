# depmixS4 1.5-1 `transition = ~ x` on a panel (N = 250, T = 4, 4 binary
# items, K = 2), cross-checked in JStats against LMest. Reference:
# JStats tests/fixtures/lta-covariate-ref.json.
source(file.path("equivalence", "lta-extensions", "common.R"))
ref <- jsonlite::fromJSON(file.path(jstats, "tests", "fixtures", "lta-covariate-ref.json"))
long <- to_long(ref$data$responses, data.frame(x = ref$data$x))
items <- paste0("y", 1:4)
long[items] <- lapply(long[items], factor)
fit <- lta(long, items, "subject", n_profiles = 2, time = "occasion",
           categorical = items, transition_covariates = "x", n_starts = 10,
           seed = 1, tol = 1e-10, max_iter = 2000)
# Conditional transition matrices at the reference grid of x.
conditional <- function(x_value) {
  vapply(1:2, function(k) {
    logit <- fit$transition_coefficients[, 1L, k, 1L] %*% c(1, x_value)
    move <- plogis(logit)
    out <- numeric(2)
    out[k] <- 1 - move
    out[3L - k] <- move
    out
  }, numeric(2)) |> t()
}
ours <- lapply(ref$expected$xGrid, conditional)
theirs <- lapply(seq_along(ref$expected$xGrid), function(g) {
  matrix(ref$expected$condTransition[g, , ], 2L)
})
# Align profiles by the response probability of category "1" per item.
ours_profile <- vapply(fit$measurement[[1L]]$response_probabilities,
                       function(p) p[, 2L], numeric(2))
gaps <- vapply(list(1:2, 2:1), function(perm) {
  max(vapply(seq_along(ours), function(g) {
    max(abs(ours[[g]][perm, perm] - theirs[[g]]))
  }, numeric(1)))
}, numeric(1))
slopes <- sort(abs(fit$transition_coefficients["x", 1L, , 1L]), decreasing = TRUE)
errors <- get_results(fit, "transition_coefficients")
slope_se <- subset(errors, term == "x")$standard_error[order(
  -abs(subset(errors, term == "x")$estimate))]
result <- rbind(
  row_result("depmix_covariate", "log likelihood", fit$log_likelihood, ref$expected$logLik, 1e-4),
  row_result("depmix_covariate", "parameter count", fit$n_parameters, ref$expected$npar, 0),
  row_result("depmix_covariate", "max |conditional transition diff|", min(gaps), 0, 1e-3),
  row_result("depmix_covariate", "max |abs slope diff|",
             max(abs(slopes - sort(abs(ref$expected$transSlopeEst), decreasing = TRUE))), 0, 1e-3),
  row_result("depmix_covariate", "max relative slope SE diff",
             max(abs(slope_se / sort(ref$expected$transSlopeSE, decreasing = TRUE)[
               order(order(-abs(ref$expected$transSlopeEst)))] - 1)), 0, 0.05))
print(result, row.names = FALSE, digits = 8)
saveRDS(result, file.path("equivalence", "lta-extensions", "result-depmix-covariate.rds"))
