# three_step(occasion =) on an lta() fit: a person-level distal outcome
# related to the profile at one occasion of the chain, or to the pair of
# profiles at two occasions. Independent references written out here:
# - each occasion's posteriors from the fitted normal densities, the pooled
#   classification-error matrix and its inverse;
# - the step-three chain's posterior of the profile(s) at the chosen
#   occasion(s) by a sum over every path of the chain (no forward-backward);
# - the ML step-three likelihood maximized by optim() over that path sum,
#   and its sandwich from numDeriv's per-person gradients and Hessian;
# - the BCH means and their sandwich from the inverse error rows.
# The Latent GOLD kit for the same model is prepared in
# equivalence/latentgold-step3-distal/ and awaits a Latent GOLD run.

distal_fixture <- function(n = 300L, occasions = 4L, seed = 41L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  state <- matrix(NA_integer_, n, occasions)
  state[, 1L] <- ifelse(stats::runif(n) < 0.55, 1L, 2L)
  for (t in seq_len(occasions)[-1L]) {  # simulation only
    moves <- stats::runif(n) < ifelse(state[, t - 1L] == 1L, 0.2, 0.15)
    state[, t] <- ifelse(moves, 3L - state[, t - 1L], state[, t - 1L])
  }
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n))
  s <- as.vector(t(state))
  data$y1 <- stats::rnorm(length(s), c(-1, 1)[s], 1)
  data$y2 <- stats::rnorm(length(s), c(-0.8, 0.8)[s], 1)
  # The distal outcome depends on the profile at occasion 3 (z3) and on the
  # pair at occasions 2 and 3 (z23); one value per person.
  z3 <- stats::rnorm(n, c(0, 1.5)[state[, 3L]], c(1, 0.8)[state[, 3L]])
  z23 <- stats::rnorm(n, c(0, 1, 2, 3)[(state[, 2L] - 1L) * 2L + state[, 3L]], 1)
  data$z3 <- rep(z3, each = occasions)
  data$z23 <- rep(z23, each = occasions)
  data$passed <- rep(ifelse(z3 + stats::rnorm(n) > 0.8, "yes", "no"), each = occasions)
  # A skipped occasion 3 for persons 1-15, sequences ending at occasion 3 for
  # persons 16-30: the chain carries them, BCH needs an assignment.
  data <- data[!(data$id %in% 1:15 & data$time == 3L), ]
  data <- data[!(data$id %in% 16:30 & data$time == 4L), ]
  row.names(data) <- NULL
  data
}

distal_fit <- function(data) {
  lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", occasions = "grid",
      n_starts = 2, seed = 1)
}

# The reference: everything below reads the fit's tables only.
distal_reference <- function(fit, data) {
  shares <- colMeans(fit$subject_posteriors[fit$informative, ])
  density <- vapply(1:2, function(k) {
    y1 <- stats::dnorm(data$y1, fit$means[k, "y1"], sqrt(fit$variances[k, "y1"]))
    y2 <- stats::dnorm(data$y2, fit$means[k, "y2"], sqrt(fit$variances[k, "y2"]))
    ifelse(is.na(y1), 1, y1) * ifelse(is.na(y2), 1, y2)
  }, numeric(nrow(data)))
  joint <- sweep(density, 2L, shares, "*")
  posteriors <- joint / rowSums(joint)
  measured <- !is.na(data$y1)
  modal <- max.col(posteriors, ties.method = "first")
  errors <- t(vapply(1:2, function(k) {
    vapply(1:2, function(j) sum(posteriors[measured & modal == j, k]), numeric(1)) /
      sum(posteriors[measured, k])
  }, numeric(2)))
  persons <- unique(data$id)
  assigned <- matrix(NA_integer_, length(persons), 4L)
  rows <- matrix(NA_integer_, length(persons), 4L)
  index <- cbind(match(data$id, persons), data$time)
  assigned[index[measured, , drop = FALSE]] <- modal[measured]
  rows[index[measured, , drop = FALSE]] <- which(measured)
  last <- tapply(data$time, match(data$id, persons), max)
  initial <- fit$initial_probabilities[1L, ]
  transition <- fit$transition_probabilities[, , 1L]
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  # Every path's prior and its assignments' probability, person by person.
  path_weight <- t(apply(assigned, 1L, function(m) {
    apply(paths, 1L, function(s) {
      prior <- initial[s[1L]] * prod(transition[cbind(s[-4L], s[-1L])])
      emitted <- prod(ifelse(is.na(m), 1, errors[cbind(s, ifelse(is.na(m), 1L, m))]))
      prior * emitted
    })
  }))
  chain <- function(occasion) {
    cells <- if (length(occasion) == 1L) paths[, occasion] else
      (paths[, occasion[1L]] - 1L) * 2L + paths[, occasion[2L]]
    weight <- vapply(sort(unique(cells)), function(cell) {
      rowSums(path_weight[, cells == cell, drop = FALSE])
    }, numeric(nrow(path_weight)))
    weight / rowSums(weight)
  }
  list(posteriors = posteriors, modal = modal, errors = errors, assigned = assigned,
       rows = rows, last = last, chain = chain, persons = persons)
}

person_outcome <- function(data, outcome, persons) {
  data[[outcome]][match(persons, data$id)]
}

# ML step three by optim() on the path sum; sandwich by numDeriv.
ml_reference <- function(weights, z, start) {
  k <- ncol(weights)
  per_person <- function(theta) {
    means <- theta[seq_len(k)]
    sds <- exp(theta[k + seq_len(k)])
    log(rowSums(weights * vapply(seq_len(k), function(j) {
      stats::dnorm(z, means[j], sds[j])
    }, numeric(length(z)))))
  }
  found <- stats::optim(start, function(theta) -sum(per_person(theta)), method = "BFGS",
                        control = list(reltol = 1e-15, maxit = 5000L))
  scores <- numDeriv::jacobian(per_person, found$par)
  information <- -numDeriv::hessian(function(theta) sum(per_person(theta)), found$par)
  inverse <- solve(information)
  g <- length(z)
  covariance <- inverse %*% crossprod(scores) %*% inverse * g / (g - 1)
  list(means = found$par[seq_len(k)], variances = exp(2 * found$par[k + seq_len(k)]),
       se = sqrt(diag(covariance))[seq_len(k)])
}

bch_reference <- function(weights, z) {
  totals <- colSums(weights)
  means <- colSums(weights * z) / totals
  influence <- sweep(weights * outer(z, means, "-"), 2L, totals, "/")
  g <- length(z)
  list(means = means, se = sqrt(colSums(influence^2) * g / (g - 1)))
}

test_that("ML at one occasion equals optim() over the sum over every path", {
  skip_if_not_installed("numDeriv")
  data <- distal_fixture()
  fit <- distal_fit(data)
  reference <- distal_reference(fit, data)
  ours <- three_step(fit, data, "z3", occasion = 3, method = "ml")
  expect_identical(attr(ours, "occasion"), 3L)
  # Every person whose sequence reaches occasion 3, including those with
  # nothing measured there (the chain infers their profile).
  kept <- reference$last >= 3L
  weights <- reference$chain(3L)[kept, ]
  z <- person_outcome(data, "z3", reference$persons)[kept]
  theirs <- ml_reference(weights, z, c(0, 1.5, 0, 0))
  expect_equal(ours$estimate, theirs$means, tolerance = 1e-6)
  expect_equal(ours$variance, theirs$variances, tolerance = 1e-5)
  expect_equal(ours$standard_error, theirs$se, tolerance = 1e-4)
  expect_equal(sum(ours$effective_n), sum(kept), tolerance = 1e-9)
})

test_that("ML at a pair of occasions equals optim() over the sum over every path", {
  skip_if_not_installed("numDeriv")
  data <- distal_fixture()
  fit <- distal_fit(data)
  reference <- distal_reference(fit, data)
  ours <- three_step(fit, data, "z23", occasion = c(2, 3), method = "ml")
  expect_identical(ours$from, c(1L, 1L, 2L, 2L))
  expect_identical(ours$to, c(1L, 2L, 1L, 2L))
  kept <- reference$last >= 3L
  weights <- reference$chain(c(2L, 3L))[kept, ]
  z <- person_outcome(data, "z23", reference$persons)[kept]
  theirs <- ml_reference(weights, z, c(0, 1, 2, 3, 0, 0, 0, 0))
  expect_equal(ours$estimate, theirs$means, tolerance = 1e-5)
  expect_equal(ours$standard_error, theirs$se, tolerance = 1e-3)
})

test_that("BCH at one occasion and at a pair equals the inverse-error weighting", {
  data <- distal_fixture()
  fit <- distal_fit(data)
  reference <- distal_reference(fit, data)
  inverse <- solve(reference$errors)
  m <- reference$assigned
  # One occasion: persons with an assignment there.
  kept <- !is.na(m[, 3L])
  z <- person_outcome(data, "z3", reference$persons)[kept]
  theirs <- bch_reference(inverse[m[kept, 3L], ], z)
  ours <- three_step(fit, data, "z3", occasion = 3)
  expect_equal(ours$estimate, theirs$means, tolerance = 1e-10)
  expect_equal(ours$standard_error, theirs$se, tolerance = 1e-10)
  # A pair: the product of the two occasions' inverse rows, later fastest.
  kept <- !is.na(m[, 2L]) & !is.na(m[, 3L])
  first <- inverse[m[kept, 2L], ]
  second <- inverse[m[kept, 3L], ]
  weights <- cbind(first[, 1] * second, first[, 2] * second)
  z <- person_outcome(data, "z23", reference$persons)[kept]
  theirs <- bch_reference(weights, z)
  ours <- three_step(fit, data, "z23", occasion = c(2, 3))
  expect_equal(ours$estimate, theirs$means, tolerance = 1e-10)
  expect_equal(ours$standard_error, theirs$se, tolerance = 1e-10)
  # Pairwise differences use the same contributions.
  pairs <- three_step(fit, data, "z23", occasion = c(2, 3), contrast = "pairs")
  expect_equal(pairs$estimate[pairs$class == 4L & pairs$reference_class == 1L],
               theirs$means[4L] - theirs$means[1L], tolerance = 1e-10)
  expect_identical(pairs$reference_from[pairs$class == 4L & pairs$reference_class == 1L], 1L)
})

test_that("the uncorrected methods use each occasion's own classification", {
  data <- distal_fixture()
  fit <- distal_fit(data)
  reference <- distal_reference(fit, data)
  m <- reference$assigned
  kept <- !is.na(m[, 3L])
  z <- person_outcome(data, "z3", reference$persons)[kept]
  modal <- three_step(fit, data, "z3", occasion = 3, method = "modal")
  expect_equal(modal$estimate, as.vector(tapply(z, m[kept, 3L], mean)), tolerance = 1e-12)
  shares <- reference$posteriors[reference$rows[kept, 3L], ]
  proportional <- three_step(fit, data, "z3", occasion = 3, method = "proportional")
  expect_equal(proportional$estimate, as.vector(colSums(shares * z) / colSums(shares)),
               tolerance = 1e-12)
})

test_that("a categorical distal outcome by ML equals EM on the path sum", {
  data <- distal_fixture()
  fit <- distal_fit(data)
  reference <- distal_reference(fit, data)
  kept <- reference$last >= 3L
  weights <- reference$chain(3L)[kept, ]
  yes <- person_outcome(data, "passed", reference$persons)[kept] == "yes"
  # Direct maximization of sum log sum_k a_ik p_k^y (1 - p_k)^(1 - y).
  objective <- function(logit) {
    p <- stats::plogis(logit)
    -sum(log(rowSums(weights * outer(yes, p, function(y, q) ifelse(y, q, 1 - q)))))
  }
  found <- stats::optim(c(0, 0), objective, method = "BFGS",
                        control = list(reltol = 1e-15, maxit = 5000L))
  ours <- three_step(fit, data, "passed", occasion = 3, method = "ml")
  expect_equal(ours$estimate[ours$category == "yes"], stats::plogis(found$par),
               tolerance = 1e-6)
})
