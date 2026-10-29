# Why test-release-engine.R's grand-scaled EVV/VVV comparison fails on Linux
# and Windows (CI on 344dc7c: posteriors 3e-5 apart against 1e-6) and passes
# on macOS. Hypothesis: at tol = 1e-10 these flat fits stop up to ~1e-4 short
# of the fixed point, so the two scalings (and two BLAS) stop at different
# points. Measured here: the posterior gap between tol 1e-10 and 1e-15 fits of
# the same data, and between raw and grand-scaled fits at each tolerance.
#   Rscript validation/matrix-refresh-2026-10-28/release-engine-tolerance.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))

gaps_for <- function(model) {
  data <- srl[seq_len(150L), c("cognitive_strategies", "intrinsic_value", "test_anxiety")]
  data$test_anxiety <- data$test_anxiety * 1000
  fit <- function(tol, ...) {
    withCallingHandlers(lpa(data, names(data), n_profiles = 2, n_starts = 3, seed = 3, model = model,
                            search = "kmeans", tol = tol, max_iter = 1e6, ...),
                        warning = function(w) invokeRestart("muffleWarning"))
  }
  raw10 <- fit(1e-10)
  raw15 <- fit(1e-15)
  scaled10 <- fit(1e-10, centering = "grand", scaling = "grand")
  scaled15 <- fit(1e-15, centering = "grand", scaling = "grand")
  # Two profiles: compare under the better of the two label orders.
  gap <- function(a, b) min(max(abs(a$subject_posteriors - b$subject_posteriors)),
                            max(abs(a$subject_posteriors - b$subject_posteriors[, 2:1])))
  data.frame(model = model,
             raw_1e10_vs_1e15 = gap(raw10, raw15),
             raw_vs_scaled_1e10 = gap(raw10, scaled10),
             raw_vs_scaled_1e15 = gap(raw15, scaled15),
             iterations_1e10 = sum(raw10$starts$iterations %||% NA),
             iterations_1e15 = sum(raw15$starts$iterations %||% NA))
}

print(do.call(rbind, lapply(c("EVV", "VVV", "VVI", "EEE"), gaps_for)), digits = 3)
