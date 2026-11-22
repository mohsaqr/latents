# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1026-categorical-r2.R,
# so that they run in CI's equivalence job rather than under R CMD check.

r2_shared_counts <- function() {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(3)
  # One overdispersed profile and one under-dispersed one: their summed score
  # at alpha = 0 is negative, yet the shared alpha's maximum is interior.
  y <- c(stats::rnbinom(200, size = 0.3, mu = 2), rep(30 + c(-1, 0, 1), length.out = 50))
  list(y = y, x = c(stats::rnorm(200, 0), stats::rnorm(50, 12)))
}

r2_shared_alpha <- function(y, classes) {
  means <- tapply(y, classes, mean)
  exp(stats::optimize(function(log_alpha) {
    sum(stats::dnbinom(y, size = exp(-log_alpha), mu = means[classes], log = TRUE))
  }, c(-8, 4), maximum = TRUE, tol = 1e-10)$maximum)
}

r2_partial_data <- function(seed = 31, n = 600) {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  u1 <- vapply(class, function(k) sample(1:3, 1L, prob = if (k == 1) c(.6, .3, .1) else
    c(.1, .3, .6)), integer(1))
  u2 <- ifelse(stats::runif(n) < 0.5, as.integer(u1 >= 2), stats::rbinom(n, 1, c(.3, .7)[class]))
  data <- data.frame(u1 = u1, u2 = u2, u3 = stats::rbinom(n, 1, c(.2, .8)[class]),
                     u4 = stats::rbinom(n, 1, c(.3, .75)[class]))
  data$u2[stats::runif(n) < ifelse(data$u1 == 3, 0.5, 0.1)] <- NA
  data$u1[sample(which(!is.na(data$u2)), 40)] <- NA
  data$u4[sample(n, 30)] <- NA
  data
}

r2_partial_fit <- function(data = r2_partial_data()) {
  lca(data, vars = names(data), n_classes = 2, local_dependence = "u1:u2",
      missing = "fiml", n_starts = 10, seed = 1, tol = 1e-13, max_iter = 20000)
}

# The observed-data log likelihood per row, by hand: a half-observed pair row
# contributes the joint probabilities that agree with its observed item.
r2_partial_rows <- function(theta, data, grid) {
  softmax <- function(v) { e <- exp(c(v, 0)); e / sum(e) }
  first <- sub(":.*$", "", grid)
  second <- sub("^.*:", "", grid)
  share <- stats::plogis(theta[1])
  at <- 1
  take <- function(m) { v <- theta[at + seq_len(m)]; at <<- at + m; v }
  blocks <- lapply(c(length(grid), 2, 2), function(K) lapply(1:2, function(k) softmax(take(K - 1))))
  logs <- vapply(1:2, function(k) {
    joint <- blocks[[1]][[k]]
    pair <- ifelse(!is.na(data$u1) & !is.na(data$u2),
      joint[match(paste(data$u1, data$u2, sep = ":"), grid)],
      ifelse(!is.na(data$u1), vapply(data$u1, function(v) sum(joint[first == as.character(v)]), 0),
        ifelse(!is.na(data$u2), vapply(data$u2, function(v) sum(joint[second == as.character(v)]), 0), 1)))
    third <- ifelse(is.na(data$u3), 1, blocks[[2]][[k]][data$u3 + 1])
    fourth <- ifelse(is.na(data$u4), 1, blocks[[3]][[k]][data$u4 + 1])
    log(c(share, 1 - share)[k]) + log(pair) + log(third) + log(fourth)
  }, numeric(nrow(data)))
  top <- apply(logs, 1, max)
  top + log(rowSums(exp(logs - top)))
}

r2_partial_theta <- function(fit) {
  logits <- function(p) log(p[-length(p)] / p[length(p)])
  c(stats::qlogis(mean(fit$subject_posteriors[, 1])),
    unlist(lapply(c("u1:u2", "u3", "u4"), function(name) {
      unlist(lapply(1:2, function(k) logits(fit$response_probabilities[[name]][k, ])))
    })))
}

shared_reference <- "independent: profile likelihood of a shared NB dispersion, stats::dnbinom() maximised by stats::optimize()"

test_that("the NB M-step keeps a shared interior dispersion that beats the floor (C2-1)", {
  counts <- r2_shared_counts()
  weight <- cbind(rep(c(1, 0), c(200, 50)), rep(c(0, 1), c(200, 50)))
  means <- colSums(weight * counts$y) / colSums(weight)
  objective <- function(point) sum(vapply(1:2, function(k) {
    sum(weight[, k] * stats::dnbinom(counts$y, size = 1 / point$dispersion[k],
                                     mu = point$means[k], log = TRUE))
  }, numeric(1)))
  # The summed overdispersion score at zero is negative here.
  expect_lt(sum(colSums(weight * (outer(counts$y, means, "-")^2 - counts$y))), 0)
  interior <- r2_shared_alpha(counts$y, rep(1:2, c(200, 50)))
  start <- list(means = means, dispersion = rep(interior, 2))
  solved <- latents:::.latents_negative_binomial_maximize(counts$y, weight, means, rep(interior, 2),
                                                          shared = TRUE)
  expect_gte(objective(solved), objective(start) - 1e-8)
  # optimize() stops at its own tolerance (1e-10 on the log scale) and the
  # M-step at its own; 1e-5 relative, the original test's tolerance.
  expect_reference(solved$dispersion, rep(interior, 2), reference = shared_reference,
                   covers = c("lpa()", "multilpa()"),
                   quantity = "shared NB dispersion of the M-step, known classes",
                   tolerance = 1e-5, measure = "relative", key = c("class 1", "class 2"),
                   note = "Package side: the internal negative-binomial M-step (.latents_negative_binomial_maximize, shared = TRUE) of lpa()/multilpa(count_dispersion = \"equal\").")
})

test_that("lpa() fits a shared NB dispersion at its interior maximum (C2-1)", {
  counts <- r2_shared_counts()
  data <- data.frame(y = counts$y, x = counts$x)
  fit <- lpa(data, vars = c("x", "y"), count = "y", count_model = "negative_binomial",
             count_dispersion = "equal", n_profiles = 2, n_starts = 5, seed = 1)
  classes <- max.col(fit$subject_posteriors)
  expect_true(fit$converged)
  # The classes are separated, so the profile likelihood on the modal classes
  # is the maximum; 1e-5 relative, the original test's tolerance.
  expect_reference(unname(fit$count_dispersion[, "y"]),
                   rep(r2_shared_alpha(counts$y, classes), 2), reference = shared_reference,
                   covers = "lpa()", quantity = "shared NB dispersion of lpa(), modal classes",
                   tolerance = 1e-5, measure = "relative", key = c("profile 1", "profile 2"))
})

test_that("lta() fits a shared NB dispersion at its interior maximum (C2-1)", {
  skip_on_cran()
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(3)
  n <- 150
  state <- rep(sample(1:2, n, TRUE, prob = c(0.75, 0.25)), each = 3)
  count <- ifelse(state == 1, stats::rnbinom(length(state), size = 0.3, mu = 2),
                  30 + sample(-1:1, length(state), TRUE))
  data <- data.frame(id = rep(seq_len(n), each = 3), time = rep(1:3, n),
                     y = stats::rnorm(length(state), c(0, 12)[state]), k = count)
  fit <- quietly(lta(data, c("y", "k"), "id", 2, time = "time", count = "k",
                     count_model = "negative_binomial", count_dispersion = "equal",
                     n_starts = 3, seed = 1), "latents_boundary")
  # The states are separated, so the hand profile likelihood on the true
  # states is the maximum. The old rule returned the 1e-8 floor and a log
  # likelihood 23.5 lower. 1e-4 relative, the original test's tolerance: the
  # hand reference uses the true states, the fit its posteriors.
  expect_reference(unname(fit$measurement[[1L]]$count_dispersion[, "k"]),
                   rep(r2_shared_alpha(data$k, state), 2), reference = shared_reference,
                   covers = "lta()", quantity = "shared NB dispersion of lta(), true states",
                   tolerance = 1e-4, measure = "relative", key = c("state 1", "state 2"))
})

test_that("a half-observed pair contributes its observed item to the likelihood (C2-2)", {
  data <- r2_partial_data()
  fit <- r2_partial_fit(data)
  expect_false(is.null(attr(fit$categorical_data, "latents_partial")))
  grid <- fit$categorical_levels[["u1:u2"]]
  theta <- r2_partial_theta(fit)
  expect_reference(as.numeric(logLik(fit)), sum(r2_partial_rows(theta, data, grid)),
                   reference = "independent: observed-data likelihood with half-observed pairs, base R (r2_partial_rows)",
                   covers = c("lca()", "logLik.multilpa()"),
                   quantity = "log likelihood with a half-observed local-dependence pair",
                   tolerance = 1e-10, measure = "relative")
  # The fit is a stationary point of the hand likelihood: the package's claim
  # is a zero gradient; numDeriv's gradient there, to the original absolute
  # bound 1e-4 (the fit stops at tol 1e-13 on the likelihood, not on the
  # gradient).
  skip_without_reference("numDeriv")
  gradient <- getExportedValue("numDeriv", "grad")(function(t) sum(r2_partial_rows(t, data, grid)), theta)
  expect_reference(rep(0, length(theta)), gradient, reference = "numDeriv", covers = "lca()",
                   quantity = "gradient of the hand likelihood at the fit",
                   tolerance = 1e-4, measure = "absolute",
                   note = "Package side: the zero gradient that the fitted maximum claims; reference: numDeriv::grad() of the independent hand likelihood r2_partial_rows().")
})

test_that("observed and robust SEs with half-observed pairs match numDeriv (C2-2)", {
  skip_on_cran()
  skip_without_reference("numDeriv")
  hessian <- getExportedValue("numDeriv", "hessian")
  jacobian <- getExportedValue("numDeriv", "jacobian")
  data <- r2_partial_data()
  fit <- r2_partial_fit(data)
  grid <- fit$categorical_levels[["u1:u2"]]
  theta <- r2_partial_theta(fit)
  observed <- solve(-hessian(function(t) sum(r2_partial_rows(t, data, grid)), theta))
  rows <- jacobian(function(t) r2_partial_rows(t, data, grid), theta)
  robust <- observed %*% crossprod(rows) %*% observed * nrow(data) / (nrow(data) - 1)
  softmax <- function(v) { e <- exp(c(v, 0)); e / sum(e) }
  natural <- function(t) {
    at <- 1
    take <- function(m) { v <- t[at + seq_len(m)]; at <<- at + m; v }
    unlist(lapply(c(length(grid), 2, 2), function(K) unlist(lapply(1:2, function(k) softmax(take(K - 1))))))
  }
  map <- jacobian(natural, theta)
  keys <- unlist(lapply(list(list("u1:u2", grid), list("u3", c("0", "1")), list("u4", c("0", "1"))),
    function(item) unlist(lapply(1:2, function(k) paste(item[[1]], k, item[[2]], sep = "|")))))
  # The delta-method standard errors of the response probabilities, from the
  # numerical information of the hand likelihood.
  reference_table <- function(covariance) {
    data.frame(key = keys, standard_error = sqrt(diag(map %*% covariance %*% t(map))))
  }
  # The package's table, keyed item|class|category like the reference.
  package_table <- function(table) {
    rows <- subset(table, parameter == "response" & level == "measurement")
    parts <- strsplit(rows$term, ":", fixed = TRUE)
    plain <- function(x) gsub("%3A", ":", x, fixed = TRUE)
    data.frame(key = paste(plain(vapply(parts, `[`, "", 1)), sub("profile_", "", rows$outcome),
                           plain(vapply(parts, `[`, "", 2)), sep = "|"),
               standard_error = rows$standard_error)
  }
  covers <- c("lca()", "parameter_inference.multilpa()")
  note <- "Reference: numDeriv::hessian()/jacobian() of the independent hand likelihood r2_partial_rows(), mapped to probabilities by the delta method."
  # 1e-6 relative per standard error, the original test's bound on |ratio - 1|.
  expect_reference(package_table(as.data.frame(parameter_inference(fit, boundary = "fix"))),
                   reference_table(observed), reference = "numDeriv", covers = covers,
                   quantity = "observed-information standard errors of the response probabilities",
                   tolerance = 1e-6, measure = "relative", key = "key", note = note)
  expect_reference(package_table(as.data.frame(parameter_inference(fit, vcov_type = "robust",
                                                                   boundary = "fix"))),
                   reference_table(robust), reference = "numDeriv", covers = covers,
                   quantity = "robust standard errors of the response probabilities",
                   tolerance = 1e-6, measure = "relative", key = "key", note = note)
})
