# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Dispersion and additive-dispersion families: the same engine with
# class-specific within-group variances. References: dense covariance
# densities, numerical derivatives, closed-form nesting relations.

families_truth <- list(
  means = rbind(c(-0.8, 0.3), c(0.9, -0.4)),
  between = rbind(c(0.3, 0.3), c(0.3, 0.3)),
  within = rbind(c(0.4, 0.5), c(1.6, 2.0)),
  weights = c(0.5, 0.5))
families_data <- local({
  # Class-specific within variances: draw each group's ratings with its
  # class's residual SDs.
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(31)
  sizes <- rep(c(4L, 6L, 8L, 10L), 40L)
  classes <- sample.int(2L, length(sizes), replace = TRUE)
  do.call(rbind, lapply(seq_along(sizes), function(j) {
    h <- classes[j]
    intercept <- stats::rnorm(2L, families_truth$means[h, ],
                              sqrt(families_truth$between[h, ]))
    ratings <- sweep(matrix(stats::rnorm(sizes[j] * 2L), sizes[j], 2L) %*%
                       diag(sqrt(families_truth$within[h, ])), 2L, intercept, "+")
    data.frame(group = sprintf("g%03d", j), class = h,
               y1 = ratings[, 1L], y2 = ratings[, 2L])
  }))
})
families_vars <- c("y1", "y2")

test_that("scores equal numerical derivatives for every family structure", {
  skip_without_reference("numDeriv")
  stats <- latents:::.additive_prepare(families_data, families_vars, "group")
  classes <- c("group_class_1", "group_class_2")
  point <- list(means = rbind(c(-0.5, 0.2), c(0.7, -0.1)),
                between = rbind(c(0.2, 0.4), c(0.35, 0.25)),
                within = rbind(c(0.6, 0.8), c(1.3, 1.7)), weights = c(0.4, 0.6))
  structures <- list(dispersion = latents:::.additive_structure("dispersion"),
                     additive_dispersion_varying = latents:::.additive_structure("additive_dispersion", "varying"),
                     additive_dispersion_equal = latents:::.additive_structure("additive_dispersion", "equal"))
  sides <- lapply(names(structures), \(name) {
    structure <- structures[[name]]
    # Project the point onto the structure: shared blocks take one row.
    shared <- function(m, mode) if (mode == "equal") m[c(1, 1), ] else m
    projected <- list(means = shared(point$means, structure$means),
                      between = shared(point$between, structure$between),
                      within = shared(point$within, structure$within),
                      weights = point$weights)
    theta <- latents:::.additive_pack(projected, structure, classes, families_vars)
    analytic <- colSums(latents:::.additive_group_scores(stats, projected, structure,
                                                         classes, families_vars))
    numeric_gradient <- numDeriv::grad(\(v) latents:::.additive_expectation(
      stats, latents:::.additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
    list(actual = data.frame(structure = name, parameter = names(theta), score = unname(analytic)),
         expected = data.frame(structure = name, parameter = names(theta), score = numeric_gradient))
  })
  # numDeriv's error is absolute near zero; a relative measure would judge
  # the rounding noise of a near-zero score. The original tolerance, 1e-7.
  expect_reference(do.call(rbind, lapply(sides, `[[`, "actual")),
                   do.call(rbind, lapply(sides, `[[`, "expected")), reference = "numDeriv",
                   covers = "multilpa()", quantity = "score off the optimum, dispersion families",
                   tolerance = 1e-7, measure = "scaled", key = c("structure", "parameter"),
                   note = "Package side: the internal analytic group scores (.additive_group_scores) behind the standard errors of multilpa(family = 'dispersion' / 'additive_dispersion'); reference: numDeriv::grad() of the package's likelihood.")
})

test_that("fitted points are stationary and EM is monotone in the new families", {
  skip_without_reference("numDeriv")
  stats <- latents:::.additive_prepare(families_data, families_vars, "group")
  cases <- list(dispersion = list("dispersion", NULL),
                additive_dispersion_varying = list("additive_dispersion", "varying"),
                additive_dispersion_equal = list("additive_dispersion", "equal"))
  gradients <- lapply(names(cases), \(name) {
    case <- cases[[name]]
    fit <- multilpa(families_data, families_vars, "group", n_group_classes = 2,
                    family = case[[1]],
                    between_variance = case[[2]] %||% "equal",
                    tol = 1e-10, seed = 3)
    expect_true(fit$converged)
    expect_false(any(fit$between_zero))
    expect_true(all(diff(fit$log_likelihood_history) >=
                      -1e-10 * (1 + abs(fit$log_likelihood))))
    theta <- coef(fit, scale = "unconstrained")
    gradient <- numDeriv::grad(\(v) latents:::.additive_expectation(
      stats, latents:::.additive_unpack(v, 2L, 2L, fit$structure))$log_likelihood, theta)
    list(actual = data.frame(family = name, parameter = names(theta), gradient = 0),
         expected = data.frame(family = name, parameter = names(theta), gradient = gradient))
  })
  # A maximum has a zero gradient: the package side is that zero, the
  # reference the numerical gradient at coef(fit). The original bound,
  # max |gradient| < 1e-3, absolute: the fit stops at tol = 1e-10 on the
  # likelihood, not on the gradient.
  expect_reference(do.call(rbind, lapply(gradients, `[[`, "actual")),
                   do.call(rbind, lapply(gradients, `[[`, "expected")), reference = "numDeriv",
                   covers = c("multilpa()", "coef.multilpa_additive()"),
                   quantity = "gradient at the fitted maximum, dispersion families",
                   tolerance = 1e-3, measure = "absolute", key = c("family", "parameter"),
                   note = "Reference: numDeriv::grad() of the package's likelihood at coef(fit, scale = 'unconstrained').")
})

test_that("shared means are pooled by precision when classes are unequal", {
  skip_without_reference("numDeriv")
  # 75/25 classes: an unweighted average of class means would not be the
  # maximum, so stationarity pins down the pooling weights.
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(77)
  sizes <- rep(c(5L, 9L), 80L)
  classes <- ifelse(stats::runif(length(sizes)) < 0.75, 1L, 2L)
  sd_within <- c(0.6, 1.8)
  data <- do.call(rbind, lapply(seq_along(sizes), function(j) {
    intercept <- stats::rnorm(1L, 0.5, 0.5)
    data.frame(group = j, y = intercept +
                 stats::rnorm(sizes[j], 0, sd_within[classes[j]]))
  }))
  # EM alone, without the Newton finish that would mask a wrong M-step: its
  # fixed point must be stationary and every step must raise the likelihood.
  stats <- latents:::.additive_prepare(data, "y", "group")
  structure <- latents:::.additive_structure("dispersion")
  start <- latents:::.additive_initialize(stats, 2L, structure, 1e-6, 1L)
  em <- latents:::.additive_em(stats, start, structure, 1e-6, 20000L, 1e-14,
                               matrix(FALSE, 2L, 1L))
  expect_true(em$converged)
  expect_gt(max(em$parameters$weights), 0.65)
  expect_true(all(diff(em$history) >= -1e-10 * (1 + abs(em$history))[-1L]))
  theta <- latents:::.additive_pack(em$parameters, structure,
                                    c("group_class_1", "group_class_2"), "y")
  gradient <- numDeriv::grad(\(v) latents:::.additive_expectation(
    stats, latents:::.additive_unpack(v, 2L, 1L, structure))$log_likelihood, theta)
  # A fixed point of EM has a zero gradient. The original bound,
  # max |gradient| < 1e-4, absolute.
  expect_reference(setNames(rep(0, length(theta)), names(theta)), setNames(gradient, names(theta)),
                   reference = "numDeriv", covers = "multilpa()",
                   quantity = "gradient at the EM fixed point, dispersion family, unequal classes",
                   tolerance = 1e-4, measure = "absolute",
                   note = "Package side: the zero gradient of the internal EM's (.additive_em) fixed point, the M-step of multilpa(family = 'dispersion'); reference: numDeriv::grad() of the package's likelihood there.")
})


test_that("class-specific within variances: compact density equals dense", {
  stats <- latents:::.additive_prepare(families_data, families_vars, "group")
  compact <- latents:::.additive_log_density(stats, families_truth)
  blocks <- additive_blocks(families_data, families_vars, "group")
  dense <- t(vapply(blocks, \(x) vapply(1:2, \(h) additive_dense_density(
    x, families_truth$means[h, ], families_truth$between[h, ],
    families_truth$within[h, ]), numeric(1)), numeric(2)))
  expect_reference(unname(compact), dense,
                   reference = "independent: dense Cholesky density (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "group log density, class-specific within variances",
                   tolerance = 1e-10, measure = "relative",
                   note = "Package side: the internal sufficient-statistic density (.additive_log_density) of multilpa(family = 'dispersion'); rows are groups, columns classes.")
})
