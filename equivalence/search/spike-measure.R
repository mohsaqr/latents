# The benchmark's spike measure (not part of the package): a fit's smallest
# class variance relative to its indicator's variance (divisor n), with full
# covariances read through their eigenvalues on the standardized scale. A
# reported solution below 0.01 is counted as a spike, as in the 2026-10-21
# benchmark, so the rates before and after the admissibility rule was removed
# are measured alike.
spike_measure <- function(state, overall) {
  if (!is.null(state$parameters) && is.list(state$parameters)) state <- state$parameters
  holders <- c(list(state), if (is.list(state$measurement)) state$measurement)
  d <- length(overall)
  scale <- 1 / sqrt(overall)
  ratios <- unlist(lapply(holders, function(holder) {
    if (!is.list(holder)) return(NULL)
    variances <- holder$variances
    from_variances <- if (is.matrix(variances) && ncol(variances) == d) {
      min(sweep(variances, 2L, overall, "/"))
    }
    covariances <- holder$covariances
    from_covariances <- if (is.array(covariances) && length(dim(covariances)) == 3L &&
                            dim(covariances)[1L] == d) {
      min(vapply(seq_len(dim(covariances)[3L]), function(class) {
        standardized <- covariances[, , class] * tcrossprod(scale)
        min(eigen(standardized, symmetric = TRUE, only.values = TRUE)$values)
      }, numeric(1)))
    }
    c(from_variances, from_covariances)
  }))
  if (length(ratios) == 0L) NA_real_ else min(ratios)
}
