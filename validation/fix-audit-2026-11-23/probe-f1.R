devtools::load_all(quiet = TRUE)
new_rule <- latents:::.latents_em_settled
old_rule <- function(gain, steps, threshold, log_likelihood) abs(gain) <= threshold
set_rule <- function(f) assignInNamespace(".latents_em_settled", f, ns = "latents")
read_mplus <- function(case, names) {
  d <- utils::read.table(sprintf("tests/equivalence/data/mplus-%s.dat", case), col.names = names)
  d[d == -999] <- NA; d
}
lca_data <- read_mplus("lca_k3", c("id", sprintf("u%d", 1:6)))
eei_data <- read_mplus("eei_k3poor", c("id", "y1", "y2", "y3", "y4"))
row <- function(rule, label, fit) {
  cl <- subset(get_classification(fit), level == "individuals")
  data.frame(rule = rule, case = label, loglik = sprintf("%.10f", as.numeric(logLik(fit))),
             iterations = fit$iterations, converged = fit$converged,
             counts = paste(sprintf("%.5f", sort(cl$estimated_n)), collapse = " "))
}
runs <- lapply(c("old", "new"), \(rule) {
  set_rule(if (rule == "old") old_rule else new_rule)
  t0 <- proc.time()[["elapsed"]]
  out <- rbind(
    row(rule, "lca_k3 tol1e-12", lca(lca_data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1, tol = 1e-12, max_iter = 100000L)),
    row(rule, "lca_k3 default", lca(lca_data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1)),
    row(rule, "lca_k3 tol1e-10", lca(lca_data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1, tol = 1e-10, max_iter = 100000L)),
    row(rule, "lca_k3 plain 1e-12", lca(lca_data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1, tol = 1e-12, max_iter = 100000L, acceleration = "none")),
    row(rule, "lca_k3 tol1e-15", lca(lca_data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1, tol = 1e-15, max_iter = 100000L)),
    row(rule, "eei tol1e-12", lpa(eei_data, c("y1","y2","y3","y4"), n_profiles = 3, model = "EEI", n_starts = 50, seed = 1, tol = 1e-12, max_iter = 100000L)),
    row(rule, "eei default", lpa(eei_data, c("y1","y2","y3","y4"), n_profiles = 3, model = "EEI", n_starts = 50, seed = 1)),
    row(rule, "esm lca default", lca(student_esm, c("time_with_friends", "on_social_media", "tv_video_games", "sports"), n_classes = 3, n_starts = 10, seed = 1)))
  message(rule, " elapsed ", proc.time()[["elapsed"]] - t0)
  out
})
print(do.call(rbind, runs), right = FALSE)
