devtools::load_all(quiet = TRUE)
.lca_design <- function(n_groups = 60, per = 8) {
  set.seed(2024)
  g <- rep(seq_len(n_groups), each = per)
  cluster <- rep(sample(1:2, n_groups, TRUE), each = per)
  z <- rep(stats::rnorm(n_groups), each = per)
  class <- ifelse(stats::runif(length(g)) < stats::plogis(ifelse(cluster == 1L, 1, -1) + 0.8 * z), 1L, 2L)
  draw <- function(p1, p2) ifelse(stats::runif(length(g)) < ifelse(class == 1L, p1, p2), 1L, 2L)
  data.frame(g = g, z = z, y1 = draw(0.85, 0.2), y2 = draw(0.8, 0.25), y3 = draw(0.75, 0.3), y4 = draw(0.2, 0.7), y5 = draw(0.9, 0.35))
}
data <- .lca_design(150, 10)
fit <- multilca(data, vars = paste0("y", 1:5), id = "g", n_profiles = 2, n_group_classes = 2, profile_covariates = "z",
                n_starts = 10, seed = 1, tol = 1e-12, max_iter = 100000L)
print(logLik(fit), digits = 12)
print(coef(fit))
print(head(get_results(fit, "probabilities"), 20))
set.seed(1)
g <- glca::glca(glca::item(y1, y2, y3, y4, y5) ~ z, group = g, data = data, nclass = 2, ncluster = 2,
                n.init = 20, maxiter = 100000, eps = 1e-12, verbose = FALSE)
cat("glca loglik", format(g$gof$loglik, digits = 13), "\n")
str(g$param)
str(g$std.err)
str(g$convergence)
# get_residuals on an LCA
lfit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"), n_classes = 2, n_starts = 2, seed = 1)
print(get_residuals(lfit))
str(get_results(lfit, "posteriors"))
print(get_results(lfit, "probabilities"))
str(student_esm[c("time_with_friends", "on_social_media", "sports")])
