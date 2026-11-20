# Model selection and its accessors: compare_models(), enumerate_lca(),
# enumerate_classes() (group-class families and transition grids),
# enumerate_regressions(), get_tree(), bootstrap_lrt() and the bootstrap table
# of parameter_inference(method = "bootstrap").
#
# References:
# - mclust on the Old Faithful eruptions (datasets::faithful): its EM run to a
#   relative change of 1e-12 (mclustBIC() for the maxima, meEEI() continued
#   from latents' posteriors for the estimates); the criteria, the likelihood
#   ratio and its chi-square p-value computed from mclust's likelihoods by
#   their definitions;
# - poLCA on its carcinoma data (the calls of
#   equivalence/enumeration-packages/compare-functions.R);
# - stats::lm() for a one-class regression; lme4::lmer() for the one-class
#   additive model (one random intercept per indicator);
# - written in base R from each candidate's coef(): the Gaussian
#   mixture-regression likelihood and Bayes' rule posteriors (and from them
#   the tree's flows), the dense additive likelihood (helper-additive-dense.R),
#   and the forward algorithm of a latent transition model.
# The bootstrap draws cannot be supplied by an argument (rule 5): only the
# deterministic parts are compared (the observed statistic, the point estimates).

selection_vars <- c("eruptions", "waiting")
selection_n <- nrow(datasets::faithful)
selection_fit <- function(n_profiles, model) {
  lpa(datasets::faithful, selection_vars, n_profiles = n_profiles, model = model,
      n_starts = if (n_profiles == 1L) 1L else 5L, seed = 1, tol = 1e-12, max_iter = 100000L)
}
selection_eei_1 <- selection_fit(1L, "EEI")
selection_eei_2 <- selection_fit(2L, "EEI")
selection_vvi_2 <- selection_fit(2L, "VVI")

# mclust's maximum and parameter count for one cell; its BIC is 2 log L - p log n.
mclust_cell <- function(model, k) {
  bic <- mclust::mclustBIC(datasets::faithful, G = k, modelNames = model,
                           control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(1e6, 1e6)),
                           verbose = FALSE)[1L, model]
  p <- mclust::nMclustParams(model, 2L, k)
  c(n_parameters = p, log_likelihood = (bic + p * log(selection_n)) / 2)
}

test_that("compare_models() reports mclust's likelihoods, criteria and likelihood-ratio test", {
  skip_without_reference("mclust")
  comparison <- compare_models(shared = selection_eei_2, varying = selection_vvi_2, test = "lrt")
  cells <- rbind(shared = mclust_cell("EEI", 2L), varying = mclust_cell("VVI", 2L))
  deviance <- -2 * cells[, "log_likelihood"]
  bic <- deviance + cells[, "n_parameters"] * log(selection_n)
  criteria <- data.frame(model = rownames(cells), n_parameters = cells[, "n_parameters"],
                         log_likelihood = cells[, "log_likelihood"],
                         aic = deviance + 2 * cells[, "n_parameters"], bic_individual = bic)
  # BIC differences and Akaike-type weights exp(-delta / 2), normalised.
  weights <- data.frame(model = rownames(cells), delta_bic = bic - min(bic),
                        bic_weight = exp(-(bic - min(bic)) / 2) / sum(exp(-(bic - min(bic)) / 2)))
  statistic <- 2 * (cells["varying", "log_likelihood"] - cells["shared", "log_likelihood"])
  df <- cells["varying", "n_parameters"] - cells["shared", "n_parameters"]
  ratio <- data.frame(model = "varying", lr_statistic = statistic, lr_df = df,
                      lr_p_value = stats::pchisq(statistic, df, lower.tail = FALSE))
  note <- "The criteria, BIC weights and the chi-square p-value are computed from mclust's log likelihoods (mclustBIC(), EM to 1e-12) by their definitions."
  check <- function(table, covers, label) {
    expect_reference(table, criteria, reference = "mclust", covers = covers,
                     quantity = paste0("parameters, log likelihood, AIC and BIC, EEI against VVI", label),
                     tolerance = 1e-8, measure = "relative", key = "model", note = note)
    expect_reference(table, weights, reference = "mclust", covers = covers,
                     quantity = paste0("BIC differences and weights, EEI against VVI", label),
                     tolerance = 1e-6, measure = "relative", key = "model", note = note)
    expect_reference(subset(table, model == "varying"), ratio, reference = "mclust", covers = covers,
                     quantity = paste0("likelihood-ratio test of EEI within VVI", label),
                     tolerance = 1e-6, measure = "relative", key = "model", note = note)
  }
  # Two EM runs stopped at a relative change of 1e-12: the maxima agree to
  # about 1e-10 relative. A difference of two likelihoods (the statistic, the
  # BIC differences and weights, the p-value) loses three digits to
  # cancellation: 1e-6 relative there.
  check(comparison, "compare_models()", "")
  check(get_results(comparison), c("compare_models()", "get_results.latents_comparison()"), ", get_results()")
})

test_that("enumerate_lca() reaches poLCA's maxima on the carcinoma data", {
  skip_without_reference("poLCA")
  carcinoma <- get(utils::data(list = "carcinoma", package = "poLCA", envir = environment()),
                   envir = environment())
  zero_one <- as.data.frame(lapply(carcinoma, \(v) as.integer(v) - min(as.integer(v))))
  grid <- quietly(enumerate_lca(zero_one, names(zero_one), n_classes = 1:3, n_starts = 20,
                                seed = 1, tol = 1e-12, max_iter = 100000L))
  formula <- stats::as.formula(sprintf("cbind(%s) ~ 1", paste(names(carcinoma), collapse = ", ")))
  # poLCA prints its progress; capture.output() keeps the log clean.
  polca <- lapply(1:3, \(k) {
    utils::capture.output(fit <- withr::with_seed(1L, poLCA::poLCA(
      formula, carcinoma, nclass = k, nrep = if (k == 1L) 1L else 30L, maxiter = 100000L,
      tol = 1e-12, verbose = FALSE, calc.se = FALSE)))
    data.frame(candidate = sprintf("classes%d", k), n_parameters = fit$npar,
               log_likelihood = fit$llik, bic_individual = fit$bic)
  })
  expected <- do.call(rbind, polca)
  covers <- c("enumerate_lca()", "as.data.frame.multilpa_enumeration()")
  # Both EM searches stop at 1e-12; the tolerance of
  # equivalence/enumeration-packages (1e-6 absolute on the likelihood, twice
  # that on the BIC).
  expect_reference(as.data.frame(grid), expected[c("candidate", "n_parameters", "log_likelihood")],
                   reference = "poLCA", covers = covers,
                   quantity = "parameters and log likelihood per candidate, carcinoma",
                   tolerance = 1e-6, measure = "absolute", key = "candidate")
  expect_reference(as.data.frame(grid), expected[c("candidate", "bic_individual")],
                   reference = "poLCA", covers = covers, quantity = "BIC per candidate, carcinoma",
                   tolerance = 2e-6, measure = "absolute", key = "candidate")
})

# One to three classes of mpg on wt (datasets::mtcars).
regression_grid <- quietly(enumerate_regressions(mpg ~ wt, datasets::mtcars, n_classes = 1:3,
                                                 n_starts = 5, seed = 1))

# Bayes' rule and the likelihood of one candidate, from its coef() alone.
regression_reference <- function(k) {
  estimates <- coef(get_results(regression_grid, "model", n_classes = k))
  at <- \(name, class) estimates[[sprintf(name, class)]]
  logits <- c(0, vapply(seq_len(k)[-1L], \(j) at("membership.class_%d.(Intercept)", j), numeric(1)))
  shares <- exp(logits) / sum(exp(logits))
  joint <- vapply(seq_len(k), \(j) {
    shares[j] * stats::dnorm(datasets::mtcars$mpg,
                             at("coefficient.class_%d.(Intercept)", j) + at("coefficient.class_%d.wt", j) * datasets::mtcars$wt,
                             exp(at("log_sigma.class_%d", j)))
  }, numeric(nrow(datasets::mtcars)))
  list(log_likelihood = sum(log(rowSums(joint))), posteriors = joint / rowSums(joint))
}

test_that("one class of enumerate_regressions() is stats::lm()", {
  one <- get_results(regression_grid, "model", n_classes = 1)
  reference <- stats::lm(mpg ~ wt, datasets::mtcars)
  covers <- c("enumerate_regressions()", "get_results.latents_regression_enumeration()")
  # The one-class M-step is least squares in closed form; 1e-8 relative.
  expect_reference(subset(as.data.frame(regression_grid), n_classes == 1, select = c(n_classes, log_likelihood)),
                   data.frame(n_classes = 1L, log_likelihood = as.numeric(stats::logLik(reference))),
                   reference = "stats", covers = c("enumerate_regressions()", "as.data.frame.latents_regression_enumeration()"),
                   quantity = "log likelihood, one class", tolerance = 1e-8, measure = "relative", key = "n_classes")
  expect_reference(coef(one)[c("coefficient.class_1.(Intercept)", "coefficient.class_1.wt")], stats::coef(reference),
                   reference = "stats", covers = covers, quantity = "coefficients, one class",
                   tolerance = 1e-8, measure = "relative", key = c("(Intercept)", "wt"))
})

test_that("the regression grid's likelihoods and criteria follow from each candidate's estimates", {
  n <- nrow(datasets::mtcars)
  expected <- do.call(rbind, lapply(1:3, \(k) {
    reference <- regression_reference(k)
    p <- 4L * k - 1L   # two coefficients and a log sigma per class, k - 1 membership logits
    z <- reference$posteriors
    entropy <- -sum(z[z > 0] * log(z[z > 0]))
    deviance <- -2 * reference$log_likelihood
    data.frame(n_classes = k, log_likelihood = reference$log_likelihood, n_parameters = p,
               aic = deviance + 2 * p, bic = deviance + p * log(n), bic_rows = deviance + p * log(n),
               sabic = deviance + p * log((n + 2) / 24), icl = deviance + p * log(n) + 2 * entropy,
               # relative entropy, 1 - E / (n log k); undefined for one class
               entropy = if (k == 1L) NA_real_ else 1 - entropy / (n * log(k)))
  }))
  reference <- "independent: Gaussian mixture-regression likelihood and Bayes' rule posteriors from coef(), base R"
  # The same sums of logs computed twice: rounding only.
  expect_reference(as.data.frame(regression_grid), expected, reference = reference,
                   covers = c("enumerate_regressions()", "as.data.frame.latents_regression_enumeration()"),
                   quantity = "likelihood, parameters, criteria and entropy per candidate",
                   tolerance = 1e-10, measure = "relative", key = "n_classes")
  expect_reference(summary(regression_grid), expected, reference = reference,
                   covers = c("enumerate_regressions()", "summary.latents_regression_enumeration()"),
                   quantity = "likelihood, parameters, criteria and entropy per candidate, summary()",
                   tolerance = 1e-10, measure = "relative", key = "n_classes")
  expect_reference(get_results(regression_grid), expected, reference = reference,
                   covers = c("enumerate_regressions()", "get_results.latents_regression_enumeration()"),
                   quantity = "likelihood, parameters, criteria and entropy per candidate, get_results()",
                   tolerance = 1e-10, measure = "relative", key = "n_classes")
})

test_that("get_tree() carries the posterior mass consecutive candidates share", {
  posteriors <- lapply(1:3, \(k) regression_reference(k)$posteriors)
  n <- nrow(datasets::mtcars)
  # The flow from class a of k classes to class b of k + 1: the mean over units
  # of p_k(a) p_(k+1)(b); the shares are the mean posteriors.
  expected <- do.call(rbind, lapply(1:2, \(k) {
    flow <- crossprod(posteriors[[k]], posteriors[[k + 1L]]) / n
    cells <- expand.grid(from_class = seq_len(k), to_class = seq_len(k + 1L))
    from_share <- colMeans(posteriors[[k]])[cells$from_class]
    to_share <- colMeans(posteriors[[k + 1L]])[cells$to_class]
    mass <- flow[cbind(cells$from_class, cells$to_class)]
    data.frame(from_candidate = sprintf("classes%d", k), from_class = cells$from_class,
               to_candidate = sprintf("classes%d", k + 1L), to_class = cells$to_class,
               from_share = from_share, to_share = to_share, mass = mass,
               share_of_from = mass / from_share, share_of_to = mass / to_share)
  }))
  # Shares of mass run down to 1e-10: absolute.
  expect_reference(get_tree(regression_grid), expected,
                   reference = "independent: Gaussian mixture-regression likelihood and Bayes' rule posteriors from coef(), base R",
                   covers = c("get_tree()", "enumerate_regressions()"),
                   quantity = "flows between consecutive candidates", tolerance = 1e-10,
                   key = c("from_candidate", "from_class", "to_candidate", "to_class"))
})

# Thirty groups of five, two group classes, a random intercept on y1.
family_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  groups <- 30L
  size <- 5L
  group <- rep(seq_len(groups), each = size)
  class <- rep(sample(1:2, groups, TRUE), each = size)
  intercept <- rep(stats::rnorm(groups, 0, 0.7), each = size)
  data.frame(g = group, y1 = c(0, 2)[class] + intercept + stats::rnorm(groups * size),
             y2 = c(0, 1.5)[class] + stats::rnorm(groups * size))
})
family_grid <- quietly(enumerate_classes(family_data, c("y1", "y2"), "g", n_group_classes = 1:2,
                                         family = "additive", between_variance = "varying",
                                         n_starts = 3, seed = 1))

test_that("one group class of the additive family is two random-intercept models (lme4)", {
  skip_without_reference("lme4")
  reference <- vapply(c("y1", "y2"), \(variable) {
    as.numeric(stats::logLik(lme4::lmer(stats::reformulate("(1 | g)", variable), family_data, REML = FALSE)))
  }, numeric(1))
  # lme4's optimiser stops on the profiled deviance to its own tolerance,
  # about 1e-7 relative: 1e-6.
  expect_reference(subset(as.data.frame(family_grid), n_group_classes == 1, select = c(candidate, log_likelihood)),
                   data.frame(candidate = "additive_groups1", log_likelihood = sum(reference)),
                   reference = "lme4", covers = c("enumerate_classes()", "as.data.frame.latents_family_enumeration()"),
                   quantity = "log likelihood, one group class: the sum of two lmer() fits", tolerance = 1e-6,
                   measure = "relative", key = "candidate")
})

test_that("the family grid's likelihoods and criteria are the dense additive model's", {
  vars <- c("y1", "y2")
  groups <- split(family_data[vars], family_data$g)
  expected <- do.call(rbind, lapply(1:2, \(k) {
    estimates <- coef(get_candidate(family_grid, n_group_classes = k))
    within <- estimates[sprintf("within_variance.%s", vars)]
    # Each group's stacked observations under each group class, by the dense
    # covariance of helper-additive-dense.R.
    joint <- matrix(vapply(seq_len(k), \(h) vapply(groups, \(x) {
      log(estimates[[sprintf("weight.group_class_%d", h)]]) +
        additive_dense_density(as.matrix(x), estimates[sprintf("mean.group_class_%d.%s", h, vars)],
                               estimates[sprintf("between_variance.group_class_%d.%s", h, vars)], within)
    }, numeric(1)), numeric(length(groups))), length(groups))
    top <- apply(joint, 1L, max)
    log_likelihood <- sum(top + log(rowSums(exp(joint - top))))
    # Per class: two means and two between variances; two within variances; k - 1 weights.
    p <- 4L * k + 2L + (k - 1L)
    data.frame(candidate = sprintf("additive_groups%d", k), log_likelihood = log_likelihood,
               n_parameters = p, aic = -2 * log_likelihood + 2 * p,
               bic = -2 * log_likelihood + p * log(length(groups)),
               bic_individual = -2 * log_likelihood + p * log(nrow(family_data)))
  }))
  reference <- "independent: dense additive likelihood from coef() (helper-additive-dense.R)"
  # The same likelihood by a dense Cholesky factorisation: rounding only.
  expect_reference(as.data.frame(family_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "as.data.frame.latents_family_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
  expect_reference(summary(family_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "summary.latents_family_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate, summary()", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
  expect_reference(get_results(family_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "get_results.latents_family_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate, get_results()", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
})

# Sixty sequences of three occasions, two states that persist with probability 0.8.
transition_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  n <- 60L
  first <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  states <- Reduce(\(previous, occasion) ifelse(stats::runif(n) < 0.8, previous, 3L - previous),
                   2:3, accumulate = TRUE, init = first)
  state <- as.vector(t(do.call(cbind, states)))
  data.frame(id = rep(seq_len(n), each = 3L), time = rep(1:3, n),
             y1 = stats::rnorm(3L * n, c(0, 2.5)[state]), y2 = stats::rnorm(3L * n, c(0, 2)[state]))
})
transition_grid <- quietly(enumerate_classes(transition_data, c("y1", "y2"), "id", time = "time",
                                             n_profiles = 2:3, n_group_classes = 1,
                                             n_starts = 3, seed = 1))

test_that("the transition grid's likelihoods and criteria are the forward algorithm's", {
  vars <- c("y1", "y2")
  sequences <- split(transition_data, transition_data$id)
  expected <- do.call(rbind, lapply(2:3, \(k) {
    estimates <- coef(get_candidate(transition_grid, n_profiles = k))
    at <- \(name, ...) estimates[[sprintf(name, ...)]]
    initial <- vapply(seq_len(k), \(j) at("profile.initial_probability.profile_%d.group_class_1", j), numeric(1))
    # transition[i, j]: from profile i to profile j.
    transition <- outer(seq_len(k), seq_len(k), Vectorize(\(i, j) {
      at("profile.transition_probability.profile_%d.group_class_1:profile_%d", j, i)
    }))
    density <- \(row) vapply(seq_len(k), \(j) prod(vapply(vars, \(v) {
      stats::dnorm(row[[v]], at("measurement.mean.profile_%d.%s", j, v),
                   sqrt(at("measurement.variance.profile_%d.%s", j, v)))
    }, numeric(1))), numeric(1))
    # Forward recursion: alpha_t = (alpha_(t-1) %*% transition) * density_t.
    log_likelihood <- sum(vapply(sequences, \(rows) {
      alpha <- Reduce(\(alpha, t) drop(alpha %*% transition) * density(rows[t, ]),
                      seq_len(nrow(rows))[-1L], init = initial * density(rows[1L, ]))
      log(sum(alpha))
    }, numeric(1)))
    # Per profile two means and two variances; k - 1 initial and k (k - 1) transition probabilities.
    p <- 4L * k + (k - 1L) + k * (k - 1L)
    data.frame(candidate = sprintf("profiles%d_groups1", k), log_likelihood = log_likelihood,
               n_parameters = p, aic = -2 * log_likelihood + 2 * p,
               bic = -2 * log_likelihood + p * log(length(sequences)),
               bic_individual = -2 * log_likelihood + p * log(nrow(transition_data)))
  }))
  reference <- "independent: forward algorithm of a hidden Markov model from coef(), base R"
  # The same likelihood by an unscaled forward recursion over three occasions: rounding only.
  expect_reference(as.data.frame(transition_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "as.data.frame.latents_transition_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
  expect_reference(summary(transition_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "summary.latents_transition_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate, summary()", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
  expect_reference(get_results(transition_grid), expected, reference = reference,
                   covers = c("enumerate_classes()", "get_results.latents_transition_enumeration()"),
                   quantity = "likelihood, parameters and criteria per candidate, get_results()", tolerance = 1e-10,
                   measure = "relative", key = "candidate")
})

test_that("bootstrap_lrt()'s observed statistic is twice mclust's log-likelihood difference", {
  skip_without_reference("mclust")
  test <- quietly(bootstrap_lrt(selection_eei_1, selection_eei_2, data = datasets::faithful,
                                iter = 9, seed = 1, n_starts = 2))
  statistic <- 2 * (mclust_cell("EEI", 2L)[["log_likelihood"]] - mclust_cell("EEI", 1L)[["log_likelihood"]])
  expected <- data.frame(alternative_profiles = 2L, statistic = statistic)
  summarised <- summary(test)
  check <- function(table, covers, label) {
    expect_reference(table, expected, reference = "mclust", covers = covers,
                     quantity = sprintf("observed LR statistic, one against two EEI profiles, %s", label),
                     tolerance = 1e-6, measure = "relative", key = "alternative_profiles",
                     note = "Twice the difference of mclust's maxima (mclustBIC(), EM to 1e-12).")
  }
  # A difference of two likelihoods, each converged to about 1e-10 relative:
  # 1e-6 relative on the statistic. The p-value rests on replicates drawn
  # inside the function (rule 5) and is not compared.
  check(as.data.frame(test), c("bootstrap_lrt()", "as.data.frame.multilpa_bootstrap_lrt()"), "as.data.frame()")
  check(get_results(test), c("bootstrap_lrt()", "get_results.multilpa_bootstrap_lrt()"), "get_results()")
  check(as.data.frame(summarised), c("summary.multilpa_bootstrap_lrt()", "as.data.frame.summary_multilpa_bootstrap_lrt()"),
        "as.data.frame(summary())")
  check(get_results(summarised), c("summary.multilpa_bootstrap_lrt()", "get_results.summary_multilpa_bootstrap_lrt()"),
        "get_results(summary())")
})

test_that("the bootstrap table's estimates are the maximum mclust's EM keeps", {
  skip_without_reference("mclust")
  table <- quietly(parameter_inference(selection_eei_2, method = "bootstrap", iter = 9, n_starts = 2, seed = 1))
  z <- matrix(get_results(selection_eei_2, "posteriors")$posterior, ncol = 2)
  continued <- mclust::meEEI(data = as.matrix(datasets::faithful), z = z,
                             control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(1e6, 1e6)))
  parameters <- continued$parameters
  expected <- rbind(
    data.frame(level = "measurement", outcome = rep(sprintf("profile_%d", 1:2), each = 2),
               term = rep(selection_vars, 2), parameter = "mean", estimate = as.vector(parameters$mean)),
    data.frame(level = "measurement", outcome = "shared", term = selection_vars, parameter = "variance",
               estimate = unname(diag(parameters$variance$Sigma))),
    data.frame(level = "profile", outcome = sprintf("profile_%d", 1:2), term = "group_class_1",
               parameter = "probability", estimate = parameters$pro))
  # Both stop at a relative change of 1e-12 or below in the log likelihood;
  # the estimates then agree to about its square root: 1e-6 relative. The
  # intervals rest on draws made inside the function (rule 5) and are not compared.
  expect_reference(subset(as.data.frame(table), level != "group"), expected, reference = "mclust",
                   covers = c("parameter_inference.multilpa()", "as.data.frame.latents_bootstrap()"),
                   quantity = "point estimates of the bootstrap table, two EEI profiles", tolerance = 1e-6,
                   measure = "relative", key = c("level", "outcome", "term", "parameter"),
                   note = "mclust::meEEI() continued from latents' posteriors.")
})
