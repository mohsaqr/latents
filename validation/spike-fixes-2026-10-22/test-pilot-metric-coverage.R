# Run from the repository root. Pure saved-record verification; no fitting.
stopifnot(file.exists("equivalence/search/analyse-pilot.R"))
out <- tempfile("latents-pilot-metrics-")
dir.create(out)
legacy <- data.frame(condition = c("lpa_vvi_k3", "lpa_vvi_k3", "lca_k3", "lpa_vvv_k2"),
  replicate = c(1L, 2L, 1L, 1L), k = 3L, fit = "kmeans_10",
  log_likelihood = c(-100, -10, -90, -50),
  min_relative_variance = c(.1, NA, NA, .001), seconds = 1)
finals <- legacy[1:3, ]
finals$fit <- "emem_10"
finals$log_likelihood <- c(-1, 0, -80)
finals$min_relative_variance <- NA_real_
saveRDS(list(rows = legacy, finals = finals), file.path(out, "legacy_001.rds"))
new <- data.frame(condition = c("custom_categorical", "custom_gaussian", "unknown", "failed"),
  replicate = 1L, k = 2L, fit = "emem_10", log_likelihood = c(-20, -30, -40, -Inf),
  min_relative_variance = c(NA, Inf, .1, .2), seconds = 1,
  has_gaussian = c(FALSE, TRUE, NA, TRUE))
saveRDS(list(rows = new, finals = NULL), file.path(out, "new_002.rds"))
log <- system2(file.path(R.home("bin"), "Rscript"),
  c("equivalence/search/analyse-pilot.R", shQuote(out)), stdout = TRUE, stderr = TRUE)
stopifnot(is.null(attr(log, "status")))
result <- readRDS(file.path(out, "summary.rds"))
rows <- result$rows
stopifnot(rows$best[1L] == -100, rows$status[1L] == "best",
  is.na(rows$best[2L]), rows$status[2L] == "unclassified",
  rows$best[3L] == -80, rows$status[3L] == "miss",
  rows$status[4L] == "spike", rows$status[5L] == "best",
  rows$status[6L] == "unclassified", rows$status[7L] == "unclassified",
  rows$status[8L] == "error", !anyNA(rows$status),
  sum(result$metric_counts$Freq[result$metric_counts$source == "finalist" &
    result$metric_counts$metric_status == "unclassified"]) == 2L)
# All-unknown legacy records and no finalists must also yield an explicit
# unclassified status when the reference is completely unavailable.
unknown <- tempfile("latents-pilot-allunknown-")
dir.create(unknown)
saveRDS(list(rows = legacy[2L, ], finals = NULL), file.path(unknown, "legacy_001.rds"))
log_unknown <- system2(file.path(R.home("bin"), "Rscript"),
  c("equivalence/search/analyse-pilot.R", shQuote(unknown)), stdout = TRUE, stderr = TRUE)
stopifnot(is.null(attr(log_unknown, "status")))
unknown_result <- readRDS(file.path(unknown, "summary.rds"))
stopifnot(unknown_result$rows$status == "unclassified", is.na(unknown_result$rows$best))

# Verify future pilot records mark failed and successful Gaussian/categorical
# fits explicitly without running the simulation driver or loading packages.
expressions <- parse("equivalence/search/pilot.R")
is_describe <- vapply(expressions, function(expression) {
  is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
    identical(expression[[2L]], as.name("describe"))
}, logical(1))
stopifnot(sum(is_describe) == 1L)
sandbox <- new.env(parent = baseenv())
sandbox$spec <- list(family = "lpa")
sandbox$spike_measure <- function(fit, overall) .2
sandbox$get_results <- function(fit, what) {
  data.frame(start = 1L, kind = "ward", stage = "final", log_likelihood = -10,
             converged = TRUE)
}
sandbox$`%||%` <- function(left, right) if (is.null(left)) right else left
eval(expressions[[which(is_describe)]], envir = sandbox)
fit <- list(subject_posteriors = matrix(1, 4L, 1L), log_likelihood = -10,
             converged = TRUE, boundary = FALSE, n_best_replicated = 1L)
data <- data.frame(x = 1:4)
gaussian <- sandbox$describe(fit, data, "x", 1, "emem_10")
failed <- sandbox$describe(simpleError("test failure"), data, "x", 1, "emem_10")
stopifnot(gaussian$row$has_gaussian, gaussian$finals$has_gaussian,
  failed$row$has_gaussian, is.na(gaussian$finals$min_relative_variance))
sandbox$spec$family <- "lca"
categorical <- sandbox$describe(fit, data, "x", 1, "emem_10")
stopifnot(!categorical$row$has_gaussian, !categorical$finals$has_gaussian)
cat("Synthetic pilot metric coverage assertions passed.\n")
