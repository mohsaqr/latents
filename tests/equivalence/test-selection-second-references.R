# Second references for model selection (the first are in test-accessors-selection.R, which holds as many
# references as one case may): Rmixmod for the Gaussian maxima behind compare_models() and bootstrap_lrt(), and
# randomLCA and glca for the latent class maxima behind enumerate_lca() (poLCA is the first). The fits are those of
# test-accessors-selection.R, on the Old Faithful eruptions (datasets::faithful) and poLCA's carcinoma data.

selection_vars <- c("eruptions", "waiting")
selection_fit <- function(n_profiles, model) {
  lpa(datasets::faithful, selection_vars, n_profiles = n_profiles, model = model,
      n_starts = if (n_profiles == 1L) 1L else 5L, seed = 1, tol = 1e-12, max_iter = 100000L)
}
selection_eei_1 <- selection_fit(1L, "EEI")
selection_eei_2 <- selection_fit(2L, "EEI")
selection_vvi_2 <- selection_fit(2L, "VVI")

# A second reference for the Gaussian maxima: Rmixmod (MIXMOD's EM, 20 random starts to a relative change of
# 1e-12). Its diagonal models with free proportions are Gaussian_pk_L_B (EEI: one volume and one shape for every
# class) and Gaussian_pk_Lk_Bk (VVI), the mapping of equivalence/lpa-packages/compare-functions.R.
rmixmod_log_likelihood <- function(model, k) {
  names <- c(EEI = "Gaussian_pk_L_B", VVI = "Gaussian_pk_Lk_Bk")
  fit <- Rmixmod::mixmodCluster(datasets::faithful, nbCluster = k,
                                models = Rmixmod::mixmodGaussianModel(listModels = names[[model]]),
                                strategy = Rmixmod::mixmodStrategy(nbTry = 20L, epsilonInAlgo = 1e-12,
                                                                   nbIterationInAlgo = 100000L),
                                seed = 1L)
  fit@bestResult@likelihood
}

test_that("compare_models() reports the log likelihoods Rmixmod reaches, and their likelihood ratio", {
  skip_without_reference("Rmixmod")
  comparison <- compare_models(shared = selection_eei_2, varying = selection_vvi_2, test = "lrt")
  likelihoods <- c(shared = rmixmod_log_likelihood("EEI", 2L), varying = rmixmod_log_likelihood("VVI", 2L))
  expected <- data.frame(model = names(likelihoods), log_likelihood = unname(likelihoods))
  # Two EM searches stopped at a relative change of 1e-12: 1e-8 relative on a maximum, as for mclust above.
  expect_reference(comparison, expected, reference = "Rmixmod", covers = "compare_models()",
                   quantity = "log likelihood, EEI against VVI", tolerance = 1e-8, measure = "relative",
                   key = "model")
  # A difference of two maxima loses three digits to cancellation: 1e-6 relative.
  expect_reference(subset(comparison, model == "varying"),
                   data.frame(model = "varying", lr_statistic = 2 * unname(diff(likelihoods))),
                   reference = "Rmixmod", covers = "compare_models()",
                   quantity = "likelihood-ratio statistic of EEI within VVI", tolerance = 1e-6, measure = "relative",
                   key = "model")
})

test_that("bootstrap_lrt()'s observed statistic is twice Rmixmod's log-likelihood difference", {
  skip_without_reference("Rmixmod")
  test <- quietly(bootstrap_lrt(selection_eei_1, selection_eei_2, data = datasets::faithful,
                                iter = 9, seed = 1, n_starts = 2))
  statistic <- 2 * (rmixmod_log_likelihood("EEI", 2L) - rmixmod_log_likelihood("EEI", 1L))
  # A difference of two maxima: 1e-6 relative, as against mclust.
  expect_reference(as.data.frame(test), data.frame(alternative_profiles = 2L, statistic = statistic),
                   reference = "Rmixmod", covers = c("bootstrap_lrt()", "as.data.frame.multilpa_bootstrap_lrt()"),
                   quantity = "observed LR statistic, one against two EEI profiles", tolerance = 1e-6,
                   measure = "relative", key = "alternative_profiles")
})

test_that("enumerate_lca() reaches randomLCA's and glca's maxima on the carcinoma data", {
  skip_without_reference("randomLCA")
  skip_without_reference("glca")
  carcinoma <- get(utils::data(list = "carcinoma", package = "poLCA", envir = environment()),
                   envir = environment())
  zero_one <- as.data.frame(lapply(carcinoma, \(v) as.integer(v) - min(as.integer(v))))
  grid <- quietly(enumerate_lca(zero_one, names(zero_one), n_classes = 1:3, n_starts = 20,
                                seed = 1, tol = 1e-12, max_iter = 100000L))
  covers <- c("enumerate_lca()", "as.data.frame.multilpa_enumeration()")
  random_lca <- do.call(rbind, lapply(1:3, \(k) {
    fit <- randomLCA::randomLCA(zero_one, nclass = k, penalty = 0, EMtol = 1e-12, notrials = 50L, seed = 1L,
                                cores = 1L, calcSE = FALSE)
    data.frame(candidate = sprintf("classes%d", k), log_likelihood = fit$logLik)
  }))
  # randomLCA stops at a change of 1e-12 in the log likelihood: 1e-6 absolute, as against poLCA.
  expect_reference(as.data.frame(grid), random_lca, reference = "randomLCA", covers = covers,
                   quantity = "log likelihood per candidate, carcinoma", tolerance = 1e-6, measure = "absolute",
                   key = "candidate")
  # glca reads the items through its own item() in the formula, with categories numbered from one (carcinoma's own
  # coding).
  formula <- stats::as.formula(sprintf("item(%s) ~ 1", paste(names(carcinoma), collapse = ", ")),
                               env = asNamespace("glca"))
  glca_fits <- do.call(rbind, lapply(2:3, \(k) {
    fit <- glca::glca(formula, data = carcinoma, nclass = k, n.init = 50L, maxiter = 100000L, eps = 1e-12,
                      seed = 1L, verbose = FALSE, std.err = FALSE)
    data.frame(candidate = sprintf("classes%d", k), log_likelihood = fit$gof$loglik)
  }))
  expect_reference(subset(as.data.frame(grid), candidate != "classes1"), glca_fits, reference = "glca",
                   covers = covers, quantity = "log likelihood per candidate of two and three classes, carcinoma",
                   tolerance = 1e-6, measure = "absolute", key = "candidate")
})
