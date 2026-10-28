Sys.setenv(VECLIB_MAXIMUM_THREADS = 1, OMP_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)
pkgload::load_all(quiet = TRUE)
library(mclust)
out <- "/private/tmp/latents-independent-spike-review-20261022"
original <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/c17a2edf-450c-4631-b62e-18ced1b5c86e/scratchpad/mclust/spikes.R"
expressions <- parse(original)
condition_assignment <- which(vapply(expressions, function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    identical(expr[[2L]], as.name("conditions"))
}, logical(1)))
eval(expressions[[condition_assignment]])
archived <- readRDS("validation/review-fit-search-2026-10-22/mclust-benchmark-inputs.rds")
source("equivalence/search/spike-measure.R")
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1L], conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
score <- function(data, fit) {
  x <- as.matrix(data)
  p <- fit$parameters
  density <- vapply(seq_along(p$pro), function(k) {
    centered <- sweep(x, 2L, p$mean[, k], "-")
    ch <- chol(p$variance$sigma[, , k])
    standardized <- forwardsolve(t(ch), t(centered))
    log(p$pro[k]) - ncol(x) * log(2 * pi)/2 - sum(log(diag(ch))) -
      colSums(standardized^2)/2
  }, numeric(nrow(x)))
  origin <- apply(density, 1L, max)
  sum(origin + log(rowSums(exp(density - origin))))
}
cells <- data.frame(condition = c("lpa_vvi_k3", "lpa_vvv_k2", "lpa_vvv_k2"),
                    replicate = c(1L, 1L, 12L), k = c(4L, 3L, 3L))
results <- lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  spec <- conditions[[cell$condition]]
  seed <- 20261021L + 1000L * spec$index + cell$replicate
  data <- simulab::simulate_lpa(n = spec$n, means = spec$means,
                                sds = spec$sds, seed = seed)
  data <- data[grep("^indicator_", names(data))]
  fits <- lapply(c(1e-5, 1e-8, 1e-10), function(tolerance) {
    captured <- capture(Mclust(data, G = cell$k, modelNames = spec$model,
                              control = emControl(tol = tolerance), verbose = FALSE))
    fit <- captured$value
    overall <- vapply(data, function(x) mean((x - mean(x))^2), numeric(1))
    independent_ratio <- min(vapply(seq_len(cell$k), function(k) {
      standardized <- fit$parameters$variance$sigma[, , k] /
        sqrt(outer(overall, overall))
      min(eigen(standardized, symmetric = TRUE, only.values = TRUE)$values)
    }, numeric(1)))
    shared_ratio <- spike_measure(list(covariances = fit$parameters$variance$sigma), overall)
    stopifnot(abs(independent_ratio - shared_ratio) < 1e-12)
    data.frame(cell, tolerance = tolerance, reported_ll = fit$loglik,
      rescored_ll = score(data, fit), ratio = independent_ratio,
      min_mass = min(colSums(fit$z)), warnings = paste(captured$warnings, collapse = ";"))
  })
  summary <- do.call(rbind, fits)
  reference <- subset(archived, condition == cell$condition &
                       replicate == cell$replicate & k == cell$k & !prior)
  stopifnot(nrow(reference) == 1L,
            abs(summary$reported_ll[1L] - reference$mclust_ll) < 1e-8,
            abs(summary$ratio[1L] - reference$ratio) < 1e-10)
  summary
})
results <- do.call(rbind, results)
print(results, row.names = FALSE, digits = 12)
write.csv(results, file.path(out, "mclust-small-check.csv"), row.names = FALSE)

# Algebraic counterexample to documentation claiming mass < parameter count
# forces singularity. Soft responsibilities preserve affine support even with
# arbitrarily small total mass.
x <- rbind(c(0, 0), c(1, 0), c(0, 1))
responsibility <- rep(.01, 3)
mean <- colSums(x * responsibility) / sum(responsibility)
centered <- sweep(x, 2L, mean, "-")
covariance <- crossprod(centered * sqrt(responsibility)) / sum(responsibility)
print(list(mass = sum(responsibility), gaussian_parameters = 5,
           covariance = covariance, eigenvalues = eigen(covariance)$values))
stopifnot(sum(responsibility) < 5, min(eigen(covariance)$values) > 0)

# Weight normalization alone changes mass when zero-weight units are appended.
w <- rep(1, 10)
normalize <- function(w) w * length(w) / sum(w)
print(list(original_mass = sum(normalize(w)),
           appended_zero_mass = sum(normalize(c(w, rep(0, 90)))),
           original_kish = sum(w)^2/sum(w^2),
           appended_kish = sum(c(w, rep(0, 90)))^2/sum(c(w, rep(0, 90))^2)))
stopifnot(sum(normalize(c(w, rep(0, 90)))) == 100)
