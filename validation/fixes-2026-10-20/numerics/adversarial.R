pkgload::load_all(quiet = TRUE)
record <- function(label, actual, expected) {
  actual <- tryCatch(force(actual), error = identity)
  if (inherits(actual, "condition")) {
    return(data.frame(case = label, actual = NA_real_, expected = expected,
                      error = conditionMessage(actual), absolute_error = NA_real_))
  }
  data.frame(case = label, actual = actual, expected = expected,
             error = NA_character_, absolute_error = abs(actual - expected))
}
# F identity supplies an independent exact reference over repeated positive
# and negative weights, avoiding comparison against another inversion routine.
cases <- expand.grid(ratio = c(1e-12, 1e-11, 1e-10, 1e-9, 1, 1e9, 1e10, 1e11,
                               1e12, 1e16, 1e24), m = c(1, 3), k = c(1, 4))
zero <- lapply(seq_len(nrow(cases)), function(index) {
  case <- cases[index, ]
  expected <- with(case, stats::pf(ratio * k / m, m, k, lower.tail = FALSE))
  record(sprintf("zero-ratio-%g-m%d-k%d", case$ratio, case$m, case$k),
         with(case, .latents_imhof_upper(0, c(rep(1, m), rep(-ratio, k)))), expected)
})
negative <- lapply(c(-100, -1, -1e-12, 0, 1e-12, 1, 100), function(q) {
  c(list(record(sprintf("single-positive-q%g", q), .latents_imhof_upper(q, 2),
                stats::pchisq(q / 2, 1, lower.tail = FALSE))),
    list(record(sprintf("single-negative-q%g", q), .latents_imhof_upper(q, -2),
                if (q < 0) stats::pchisq(-q / 2, 1) else 0)))
})
convolution <- function(q, a, b) {
  f <- function(t) stats::pchisq(pmax((q - b * t^2) / a, 0), 1, lower.tail = FALSE) *
    2 * stats::dnorm(t)
  kink <- if (q < 0) sqrt(q / b) else 0
  2 * stats::pnorm(kink) - 1 +
    stats::integrate(f, kink, Inf, rel.tol = 1e-11, abs.tol = 1e-13)$value
}
grid <- expand.grid(q = c(-5, -1e-12, 1e-12, 5), b = c(-0.1, -1, -1e4))
mixed <- lapply(seq_len(nrow(grid)), function(index) {
  q <- grid$q[index]
  b <- grid$b[index]
  record(sprintf("mixed-q%g-b%g", q, b), .latents_imhof_upper(q, c(1, b)),
         convolution(q, 1, b))
})
result <- do.call(rbind, c(zero, unlist(negative, recursive = FALSE), mixed))
print(result, row.names = FALSE, digits = 12)
utils::write.csv(result, "validation/fixes-2026-10-20/numerics/adversarial.csv", row.names = FALSE)
cat("Errors:", sum(!is.na(result$error)), " Maximum absolute error:",
    max(result$absolute_error, na.rm = TRUE), "\n")
stopifnot(all(is.na(result$error)), max(result$absolute_error) < 1e-10)
