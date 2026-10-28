pkgload::load_all('.', quiet = TRUE)
fit <- readRDS('validation/release-audit-2026-10-28/inference/units-probe.rds')[[1L]]$fit
result <- lapply(c('observed', 'robust', 'opg'), function(type) {
  info <- latents:::.lta_inference(fit, type)
  keep <- info$identified
  information <- info$information[keep, keep]
  scores <- info$group_scores[, keep]
  inverse <- solve(information)
  expected <- switch(type, observed = inverse,
    robust = inverse %*% crossprod(scores) %*% inverse * fit$n_groups/(fit$n_groups - 1),
    opg = solve(crossprod(scores)))
  obtained <- info$vcov[keep, keep]
  difference <- max(abs(obtained - expected)/sqrt(outer(diag(expected), diag(expected))))
  stopifnot(difference < 1e-10)
  data.frame(type = type, normalized_covariance_error = difference)
})
print(do.call(rbind, result))
