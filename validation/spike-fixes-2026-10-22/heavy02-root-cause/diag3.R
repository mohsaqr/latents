# Is heavy_02's lower point a maximum or a saddle? Hessian eigenvalues of the
# independent log likelihood, and plain EM continued from it with and without
# a small nudge.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
sys.source(file.path(root, "reference.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars; data <- case$data
fits <- readRDS("/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad/heavy/fits.rds")
effects <- case$spec$effects
points <- lapply(fits, ref_point_from_fit, effects = effects)
cat("LL check: ", vapply(points, ref_log_likelihood, numeric(1), data = data, vars = vars), "\n")
curvature <- function(point) {
  values <- ref_pack(point, "VVI")
  objective <- function(v) ref_log_likelihood(ref_unpack(v, point, "VVI"), data, vars)
  centre <- objective(values)
  scale <- vapply(seq_along(values), function(i) {
    step <- 1e-3 * max(abs(values[[i]]), 0.1)
    up <- values; up[[i]] <- up[[i]] + step
    down <- values; down[[i]] <- down[[i]] - step
    1 / sqrt(abs((objective(up) - 2 * centre + objective(down)) / step^2))
  }, numeric(1))
  standardized <- function(u) objective(values + u * scale)
  settings <- list(eps = 0.1, d = 0.1, r = 4, zero.tol = 1e-8)
  zero <- numeric(length(values))
  hessian <- numDeriv::hessian(standardized, zero, method.args = settings)
  gradient <- numDeriv::grad(standardized, zero, method.args = settings)
  e <- eigen((hessian + t(hessian)) / 2, symmetric = TRUE)
  list(eigenvalues = e$values, top_vector = stats::setNames(e$vectors[, 1L] * scale, names(values)),
       gradient = max(abs(gradient)))
}
invisible(lapply(names(points), function(name) {
  result <- curvature(points[[name]])
  cat(sprintf("\n%s: max |score| %.2e; largest Hessian eigenvalues: %s; smallest: %s\n", name,
              result$gradient, paste(sprintf("%.4g", head(result$eigenvalues, 3)), collapse = ", "),
              paste(sprintf("%.4g", tail(result$eigenvalues, 2)), collapse = ", ")))
  if (name == "low") print(round(result$top_vector, 4))
}))
em_from <- function(point, steps) {
  history <- numeric(steps)
  final <- Reduce(function(p, i) {
    q <- ref_m_step(ref_posteriors(p, data, vars), p, data, vars, "VVI")
    history[[i]] <<- ref_log_likelihood(q, data, vars)
    q
  }, seq_len(steps), point)
  list(point = final, history = history)
}
plain <- em_from(points$low, 3000)
cat("\nplain EM from low, LL at steps 1, 10, 100, 1000, 2000, 3000:\n")
print(plain$history[c(1, 10, 100, 1000, 2000, 3000)], digits = 12)
set.seed(7)
nudged <- vapply(1:6, function(r) {
  values <- ref_pack(points$low, "VVI")
  values <- values * (1 + stats::rnorm(length(values), 0, 0.01))
  run <- em_from(ref_unpack(values, points$low, "VVI"), 1500)
  run$history[[1500]]
}, numeric(1))
cat("\nEM after 1% random nudges of low (1500 steps):\n"); print(nudged, digits = 10)
