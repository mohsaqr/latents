# Simulated data and block layouts for the additive-engine tests. The dense
# reference densities are in tests/equivalence/helper-additive-dense.R.

#' Simulate raw ratings from the additive model
#' @return A data frame with `group`, `y1`, ..., one row per observation.
additive_draw <- function(parameters, sizes, seed) {
  # Scope the seed: restore the caller's RNG state on exit.
  old_seed <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  on.exit(if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(seed)
  d <- ncol(parameters$means)
  classes <- sample.int(nrow(parameters$means), length(sizes), replace = TRUE,
                        prob = parameters$weights)
  blocks <- lapply(seq_along(sizes), function(j) {
    intercept <- stats::rnorm(d, parameters$means[classes[j], ],
                              sqrt(parameters$between[classes[j], ]))
    noise <- matrix(stats::rnorm(sizes[j] * d), sizes[j], d) %*%
      diag(sqrt(parameters$within), d)
    ratings <- sweep(noise, 2L, intercept, "+")
    colnames(ratings) <- paste0("y", seq_len(d))
    data.frame(group = sprintf("g%02d", j), ratings)
  })
  do.call(rbind, blocks)
}

#' Split a rating frame into one matrix per group, in first-occurrence order.
additive_blocks <- function(data, vars, id) {
  groups <- unique(data[[id]])
  lapply(groups, function(g) as.matrix(subset(data, data[[id]] == g)[, vars]))
}
