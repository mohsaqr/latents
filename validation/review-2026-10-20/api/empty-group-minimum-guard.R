library(latents, lib.loc = "validation/review-2026-10-20/rlib")
set.seed(792)
fit <- function(data) suppressMessages(multilpa(data, "y", "g", n_profiles = 1,
  n_group_classes = 1, missing = "fiml", n_starts = 1, seed = 3))
evaluate <- function(fitted, type) {
  tryCatch(parameter_inference(fitted, vcov_type = type), error = identity)
}
invisible(lapply(1:2, function(ng) {
  d <- data.frame(g = rep(seq_len(ng), each = 50), y = rnorm(50 * ng))
  padded <- rbind(d, data.frame(g = ng + 1:95, y = NA_real_))
  a <- fit(d)
  b <- fit(padded)
  cat("\n=== Informative groups:", ng, "===\n")
  cat("Original groups/parameters:", a$n_groups, a$n_parameters, "\n")
  cat("Padded groups/informative/parameters:", b$n_groups, b$n_informative_groups,
      b$n_parameters, "\n")
  cat("Log likelihood:", a$log_likelihood, b$log_likelihood, "\n")
  stopifnot(abs(a$log_likelihood - b$log_likelihood) < 1e-9)
  invisible(lapply(c("robust", "opg"), function(type) {
    baseline <- evaluate(a, type)
    extra <- evaluate(b, type)
    cat("Type:", type, "\n")
    cat("Baseline:", paste(class(baseline), collapse = ", "), "\n")
    if (inherits(baseline, "condition")) cat(conditionMessage(baseline), "\n")
    cat("Padded:", paste(class(extra), collapse = ", "), "\n")
    if (inherits(extra, "condition")) cat(conditionMessage(extra), "\n") else print(extra)
    stopifnot(inherits(baseline, "latents_too_few_groups"))
    if (type == "robust") stopifnot(is.data.frame(extra))
  }))
}))
