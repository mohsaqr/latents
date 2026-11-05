# The categorical boundary flag against poLCA (Linzer & Lewis 2011). Five
# binary items are fitted by both programs, twice: once with two items that
# copy the class, so every candidate has response probabilities at 0 or 1,
# and once with none, so the two-class model is interior and the three-class
# model is not. poLCA has no probability floor: a probability at the bound is
# exactly 0 or about 1e-40 there. It is called a bound at the package's own
# reporting floor, 1e-6 (.multilpa_probability_floor()).
#
# The package keeps probabilities at least min_probability = 1e-10 from 0 and
# 1, which moves a bounded fit's log likelihood by at most
# n x items x 1e-10 = 250 x 5 x 1e-10 = 1.25e-7 (observed: 5e-8). The
# likelihoods are therefore compared to 1e-6.

boundary_items <- function(seed, exact) {
  set.seed(seed)
  n <- 250L
  class <- sample(1:2, n, TRUE)
  draw <- function(p) ifelse(stats::runif(n) < p[class], "yes", "no")
  data.frame(a = if (exact) ifelse(class == 1L, "yes", "no") else draw(c(.8, .25)),
             b = if (exact) ifelse(class == 1L, "up", "down") else draw(c(.75, .2)),
             c = draw(c(.85, .2)), d = draw(c(.9, .3)), e = draw(c(.75, .2)))
}

# One candidate fitted by both programs, from many starts, to a tight tolerance.
boundary_pair <- function(data, k) {
  fit <- lca(data, names(data), k, n_starts = 20, seed = 1, tol = 1e-12,
             max_iter = 20000)
  coded <- as.data.frame(lapply(data, \(v) as.integer(factor(v))))
  set.seed(7)
  reference <- poLCA::poLCA(cbind(a, b, c, d, e) ~ 1, coded, nclass = k, nrep = 30,
                            verbose = FALSE, maxiter = 50000, tol = 1e-14,
                            calc.se = FALSE)
  list(fit = fit, reference = reference, data = data)
}

test_that("the categorical boundary flag agrees with poLCA at the same maximum", {
  skip_if_not_installed("poLCA")
  cases <- data.frame(seed = c(1106, 1106, 2211, 2211), exact = c(TRUE, TRUE, FALSE, FALSE),
                      k = c(2L, 3L, 2L, 3L))
  pairs <- Map(\(seed, exact, k) boundary_pair(boundary_items(seed, exact), k),
               cases$seed, cases$exact, cases$k)
  floor <- 1e-6
  latents_flag <- vapply(pairs, \(pair) .latents_boundary_active(pair$fit), logical(1))
  polca_flag <- vapply(pairs, \(pair) min(unlist(pair$reference$probs)) <= floor,
                       logical(1))
  expect_identical(latents_flag, polca_flag)
  # The cases include both outcomes, so agreement is not two constant flags.
  expect_identical(latents_flag, c(TRUE, TRUE, FALSE, TRUE))
  lapply(pairs, \(pair) {
    expect_lt(abs(pair$fit$log_likelihood - pair$reference$llik), 1e-6)
    expect_identical(as.integer(pair$fit$n_parameters), as.integer(pair$reference$npar))
  })
  # The grid reports the flag the fits carry.
  interior <- boundary_items(2211, FALSE)
  grid <- enumerate_lca(interior, names(interior), n_classes = 2:3, n_starts = 20,
                        seed = 1, tol = 1e-12, max_iter = 20000)
  expect_identical(get_candidates(grid)$boundary, polca_flag[3:4])
})

test_that("the interior two-class solution equals poLCA's estimates", {
  skip_if_not_installed("poLCA")
  pair <- boundary_pair(boundary_items(2211, FALSE), 2L)
  # Classes matched by size; categories in factor() order, as coded for poLCA.
  ours <- order(pair$fit$profile_probabilities)
  theirs <- order(pair$reference$P)
  expect_equal(unname(pair$fit$profile_probabilities[ours]),
               unname(pair$reference$P[theirs]), tolerance = 1e-6)
  lapply(names(pair$data), \(item) {
    categories <- levels(factor(pair$data[[item]]))
    expect_equal(unname(pair$fit$response_probabilities[[item]][ours, categories]),
                 unname(pair$reference$probs[[item]][theirs, ]), tolerance = 1e-6)
  })
})
