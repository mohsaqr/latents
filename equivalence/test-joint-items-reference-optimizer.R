# Independent optimizer contract: retain a complete best point, not a lower last step.

ji_reference_environment <- function() {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(file.path("features", "joint-items", "reference.R"), envir = environment)
  environment
}

test_that("independent EM retains coherent best states and passes negative numerical steps", {
  environment <- ji_reference_environment()
  case <- readRDS(file.path("features", "joint-items", "kit", "cases.rds"))$huge_01
  coded <- environment$ref_code(case$data, case$spec$vars, case$spec$pairs)
  # Own selected strict fit from seed1:2/15 starts; stored only as an optimizer
  # regression point, never supplied to production fitting or comparison.
  start <- list(priors = c(profile_1 = 0.348694520745934, profile_2 = 0.391087642064523,
profile_3 = 0.260217837189543), probabilities = list(structure(c(0.0831329364715534,
0.48124440098418, 0.209260299406478, 0.0270543958349065, 0.073282568690033,
0.0534412766129096, 0.0458646270039591, 0.0424126177774554, 0.0454242695848892,
0.0394691142853711, 0.0602975088493613, 0.0711162837705654, 0.200190387172842,
0.215558139643673, 0.275125583622522, 0.101949505648072, 0.0203231536243177,
0.0632263796179662, 0.0597150865451779, 0.0238168028348363, 0.0398251686570101,
0.0970964291958771, 0.0111666982261091, 0.0489046009065038, 0.345527517842241,
0.0718981093700339, 0.193676137821155), dim = c(3L, 9L)), structure(c(0.20883390748885,
0.853036413751767, 0.792292962273534, 0.79116609251115, 0.146963586248233,
0.207707037726466), dim = 3:2), structure(c(0.244402586112791,
0.789300872396838, 0.75856575564956, 0.755597413887209, 0.210699127603162,
0.24143424435044), dim = 3:2), structure(c(0.293830622826632,
0.755498222740897, 0.219726879197973, 0.706169377173368, 0.244501777259103,
0.780273120802027), dim = 3:2), structure(c(0.150920795270859,
0.804824346991215, 0.252165225067617, 0.849079204729141, 0.195175653008785,
0.747834774932382), dim = 3:2), structure(c(0.19885118460236,
0.849387855752601, 0.794628800293088, 0.80114881539764, 0.150612144247399,
0.205371199706913), dim = 3:2), structure(c(0.254054594447607,
0.757746789080023, 0.217575976089875, 0.745945405552393, 0.242253210919977,
0.782424023910125), dim = 3:2)))
  initial <- environment$ref_loglik(coded$codes, start$priors, start$probabilities)
  first_step <- environment$ref_em(coded$codes, coded$levels, 3L, start = start, max_iter = 2L)
  expect_lt(first_step$termination_gain, -1e-9)
  expect_lt(first_step$termination_log_likelihood, initial)
  expect_identical(first_step$log_likelihood, initial)
  expect_identical(first_step$priors, start$priors)
  expect_identical(first_step$probabilities, start$probabilities)
  expect_identical(first_step$best_iteration, 1L)
  expect_identical(first_step$iterations, 2L)
  expect_false(first_step$converged)
  fit <- environment$ref_em(coded$codes, coded$levels, 3L, start = start, max_iter = 1000L)
  expect_gte(fit$log_likelihood, initial)
  expect_gt(fit$iterations, 2L)
  expect_lte(fit$best_iteration, fit$iterations)
  expect_gte(fit$termination_gain, 0)
  expect_lt(fit$termination_gain, 1e-12)
  expect_true(fit$converged)
  expect_identical(environment$ref_loglik(coded$codes, fit$priors, fit$probabilities),
                   fit$log_likelihood)
  joint <- environment$ref_log_joint(coded$codes, fit$priors, fit$probabilities)
  expect_identical(fit$posteriors, exp(joint - environment$ref_log_sum_exp(joint)))
})

test_that("independent EM retains ordinary improving states and iteration limits", {
  environment <- ji_reference_environment()
  withr::local_seed(632)
  latent <- rep(c(1L, 2L), each = 100L)
  codes <- vapply(seq_len(4L), function(v) {
    1L + as.integer(runif(200L) < ifelse(latent == 1L, .15, .85))
  }, integer(200L))
  levels <- rep(list(1:2), 4L)
  start <- list(priors = c(.4, .6), probabilities = rep(list(
    matrix(c(.6, .4, .4, .6), 2L)), 4L))
  initial <- environment$ref_loglik(codes, start$priors, start$probabilities)
  fit <- environment$ref_em(codes, levels, 2L, start = start, tol = 1e-10)
  expect_gt(fit$log_likelihood, initial)
  expect_gt(fit$best_iteration, 1L)
  expect_lte(fit$best_iteration, fit$iterations)
  expect_true(fit$converged)
  expect_identical(environment$ref_loglik(codes, fit$priors, fit$probabilities),
                   fit$log_likelihood)
  joint <- environment$ref_log_joint(codes, fit$priors, fit$probabilities)
  expect_identical(fit$posteriors, exp(joint - environment$ref_log_sum_exp(joint)))
  limited <- environment$ref_em(codes, levels, 2L, start = start, max_iter = 1L)
  expect_identical(limited$log_likelihood, initial)
  expect_identical(limited$best_iteration, 1L)
  expect_identical(limited$iterations, 1L)
  expect_false(limited$converged)
  expect_error(environment$ref_em(matrix(1L, 2L, 1L), list(1L), 1L,
    start = list(priors = 1, probabilities = list(matrix(0, 1L, 1L)))),
    "non-finite likelihood or posteriors", fixed = TRUE)
})
