# Shared by the growth mixture comparisons against lcmm::hlme
# (test-growth-mixture-lcmm.R, test-growth-predict-lcmm.R).

# latents' parameter list from an hlme fit, classes in hlme's order.
params_from_hlme <- function(fit, model) {
  spec <- fit$spec
  n_classes <- spec$n_classes
  best <- model$best
  counts <- model$N
  at <- 0L
  take <- function(count) {
    values <- best[at + seq_len(count)]
    at <<- at + count
    values
  }
  membership <- take(counts[1L])
  fixed <- take(counts[2L])
  varcov <- take(counts[3L])
  proportional <- take(counts[4L])
  # hlme(cor = AR(time)): cor1 is the rate, cor2 the standard deviation of
  # the process; cor = BM(time): cor1 is its standard deviation.
  correlation <- if (length(counts) >= 5L) take(counts[5L]) else numeric()
  params <- fit$params
  # Fixed effects: hlme lists each class-specific term once per class
  # ("time class1"), or once without a suffix when there is one class.
  label_of <- function(term) if (identical(term, "(Intercept)")) "intercept" else term
  params$beta[] <- vapply(seq_len(n_classes), function(k) {
    vapply(rownames(params$beta), function(term) {
      name <- if (n_classes == 1L) label_of(term) else
        sprintf("%s class%d", label_of(term), k)
      unname(fixed[[name]])
    }, numeric(1))
  }, numeric(nrow(params$beta)))
  if (length(params$common) > 0L) {
    params$common[] <- unname(fixed[names(params$common)])
  }
  n_random <- ncol(spec$random_design)
  covariance <- if (isTRUE(spec$random_diagonal)) diag(varcov, n_random) else {
    m <- matrix(0, n_random, n_random)
    m[upper.tri(m, diag = TRUE)] <- varcov
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  params$random_covariance <- list(covariance)
  params$random_scale <- if (length(proportional) > 0L)
    c(unname(proportional), 1) else rep(1, n_classes)
  params$sigma2 <- rep(unname(best[["stderr"]])^2, n_classes)
  if (length(correlation) == 2L) {
    params$correlation <- c(rate = abs(unname(correlation[1L])),
                            sd = abs(unname(correlation[2L])))
  } else if (length(correlation) == 1L) {
    params$correlation <- c(sd = abs(unname(correlation[1L])))
  }
  # hlme's membership logits are against its last class; latents' against
  # the first.
  if (n_classes > 1L) {
    terms <- colnames(spec$w)
    logits <- vapply(seq_len(n_classes), function(k) {
      if (k == n_classes) return(numeric(length(terms)))
      unname(membership[sprintf("%s class%d", ifelse(terms == "(Intercept)",
                                                     "intercept", terms), k)])
    }, numeric(length(terms))) |> matrix(length(terms))
    params$gamma[] <- logits - logits[, 1L]
  }
  params
}
