# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# The accessors of the additive family (multilpa(family = "additive"), Houle
# et al., 2026) and of the full cross-level family. References: for the
# additive fit, the dense-covariance likelihood of helper-additive-dense.R,
# its numDeriv Hessian inverted (standard errors), Wald intervals from their
# definitions, and Bayes' rule on the dense densities (group posteriors); for
# the cross-level fit, the published manifest likelihood written in this file,
# its maximum continued by stats::optim() from latents' estimates, and Bayes'
# rule on it (group posteriors). nobs() is the number of groups in the data.

# ---- additive family ---------------------------------------------------------

# 160 groups of four ratings, two well separated classes of about 80 groups,
# so the weak-class warning stays silent.
accessor_additive_data <- additive_draw(
  list(means = rbind(c(-1.5, -1.5), c(1.5, 1.5)), between = rbind(c(0.25, 0.5), c(0.5, 0.25)),
       within = c(1, 1), weights = c(0.5, 0.5)),
  rep(4L, 160L), seed = 1)
accessor_additive_fit <- quietly(multilpa(accessor_additive_data, c("y1", "y2"), "group",
                                          n_group_classes = 2, family = "additive",
                                          tol = 1e-10, seed = 3))
accessor_additive_coef <- coef(accessor_additive_fit)
accessor_additive_blocks <- additive_blocks(accessor_additive_data, c("y1", "y2"), "group")

# The free natural parameters in coef() order: four class means (class 1 y1,
# y2, class 2 y1, y2), two within variances, four between variances, the
# weight of class 1 (class 2's is one minus it).
accessor_additive_free <- names(accessor_additive_coef)[1:11]

# Per-group joint log densities (group x class) from the dense covariance.
accessor_additive_joint <- function(theta) {
  means <- matrix(theta[1:4], 2L, byrow = TRUE)
  between <- matrix(theta[7:10], 2L, byrow = TRUE)
  weights <- c(theta[[11L]], 1 - theta[[11L]])
  t(vapply(accessor_additive_blocks, \(x) vapply(1:2, \(h) {
    log(weights[h]) + additive_dense_density(x, means[h, ], between[h, ], theta[5:6])
  }, numeric(1)), numeric(2)))
}
accessor_additive_loglik <- function(theta) {
  joint <- accessor_additive_joint(theta)
  top <- apply(joint, 1L, max)
  sum(top + log(rowSums(exp(joint - top))))
}

# What each coef() name is in the parameter tables.
accessor_additive_rows <- data.frame(
  name = names(accessor_additive_coef),
  level = c(rep("between", 4L), "within", "within", rep("between", 4L), "group_class", "group_class"),
  group_class = c("group_class_1", "group_class_1", "group_class_2", "group_class_2", "shared", "shared",
                  "group_class_1", "group_class_1", "group_class_2", "group_class_2",
                  "group_class_1", "group_class_2"),
  indicator = c("y1", "y2", "y1", "y2", "y1", "y2", "y1", "y2", "y1", "y2", NA, NA),
  parameter = c(rep("mean", 4L), rep("variance", 6L), "weight", "weight"))

# Wald inference from standard errors: means on their own scale with a z test
# of zero; variances on the log scale and weights on the logit scale (delta
# method), back-transformed, with no test.
accessor_additive_wald <- function(errors) {
  rows <- accessor_additive_rows
  estimate <- unname(accessor_additive_coef[rows$name])
  error <- unname(errors[rows$name])
  z <- stats::qnorm(0.975)
  is_mean <- rows$parameter == "mean"
  interval <- function(estimate, error, parameter, sign) {
    switch(parameter,
           mean = estimate + sign * z * error,
           variance = estimate * exp(sign * z * error / estimate),
           weight = stats::plogis(stats::qlogis(estimate) + sign * z * error / (estimate * (1 - estimate))))
  }
  statistic <- ifelse(is_mean, estimate / error, NA_real_)
  data.frame(rows[c("level", "group_class", "indicator", "parameter")],
             estimate = estimate, standard_error = error, statistic = statistic,
             p_value = ifelse(is_mean, 2 * stats::pnorm(-abs(statistic)), NA_real_),
             conf_low = unlist(Map(interval, estimate, error, rows$parameter, -1)),
             conf_high = unlist(Map(interval, estimate, error, rows$parameter, 1)))
}

test_that("the dense likelihood at coef() is the additive fit's log likelihood", {
  # The same Gaussian likelihood by Cholesky factors of the full covariance: rounding only.
  expect_reference(as.numeric(logLik(accessor_additive_fit)),
                   accessor_additive_loglik(unname(accessor_additive_coef[accessor_additive_free])),
                   reference = "independent: dense-covariance additive likelihood (helper-additive-dense.R)",
                   covers = c("multilpa()", "logLik.multilpa_additive()", "coef.multilpa_additive()"),
                   quantity = "additive log likelihood at coef()", tolerance = 1e-10, measure = "relative")
})

test_that("the additive parameter tables are Wald inference from numDeriv's Hessian of the dense likelihood", {
  skip_without_reference("numDeriv")
  covariance <- solve(-numDeriv::hessian(accessor_additive_loglik,
                                         unname(accessor_additive_coef[accessor_additive_free])))
  errors <- stats::setNames(sqrt(diag(covariance)), accessor_additive_free)
  errors[["weight.group_class_2"]] <- errors[["weight.group_class_1"]]
  wald <- accessor_additive_wald(errors)
  # numDeriv's Richardson extrapolation gives the Hessian to about 1e-8
  # relative and the package's information is a central difference of its
  # analytic scores (step 1e-4): 1e-6 relative on standard errors and on the
  # intervals and tests built from them, as in test-accessors-lpa.R.
  check <- function(table, covers, label) {
    expect_reference(subset(table, parameter != "weight"), subset(wald, parameter != "weight"),
                     reference = "numDeriv", covers = covers,
                     quantity = sprintf("Wald inference for means and variances, %s", label),
                     tolerance = 1e-6, measure = "relative",
                     key = c("level", "group_class", "indicator", "parameter"),
                     note = "Reference: standard errors from the inverse of numDeriv::hessian() of the dense-covariance likelihood (helper-additive-dense.R) in the natural parameters; intervals: means estimate +/- z SE, variances on the log scale, from their definitions.")
    expect_reference(subset(table, parameter == "weight", select = -indicator),
                     subset(wald, parameter == "weight", select = -indicator),
                     reference = "numDeriv", covers = covers,
                     quantity = sprintf("Wald inference for class weights, %s", label),
                     tolerance = 1e-6, measure = "relative", key = c("level", "group_class", "parameter"),
                     note = "Reference: standard errors from the inverse of numDeriv::hessian() of the dense-covariance likelihood; intervals on the logit scale by the delta method.")
  }
  check(parameter_inference(accessor_additive_fit),
        c("multilpa()", "parameter_inference()", "parameter_inference.multilpa_additive()"),
        "parameter_inference()")
  check(as.data.frame(accessor_additive_fit),
        c("multilpa()", "as.data.frame.multilpa_additive()"), "as.data.frame()")
  check(get_results(accessor_additive_fit, "parameters"),
        c("multilpa()", "get_results.multilpa_additive()"), "get_results()")
  check(as.data.frame(summary(accessor_additive_fit)),
        c("summary.multilpa_additive()", "as.data.frame.summary_multilpa_additive()"),
        "as.data.frame(summary())")
  check(get_results(summary(accessor_additive_fit)),
        c("summary.multilpa_additive()", "get_results.summary_multilpa_additive()"),
        "get_results(summary())")
  intervals <- cbind(`2.5 %` = wald$conf_low, `97.5 %` = wald$conf_high)
  rownames(intervals) <- accessor_additive_rows$name
  expect_reference(confint(accessor_additive_fit), intervals, reference = "numDeriv",
                   covers = c("multilpa()", "confint.multilpa_additive()"),
                   quantity = "95% intervals of the additive parameters", tolerance = 1e-6, measure = "relative",
                   note = "Reference: the Wald intervals above, from numDeriv's Hessian of the dense likelihood.")
})

test_that("the additive group posteriors are Bayes' rule on the dense densities", {
  joint <- accessor_additive_joint(unname(accessor_additive_coef[accessor_additive_free]))
  posterior <- exp(joint - apply(joint, 1L, max))
  posterior <- posterior / rowSums(posterior)
  expected <- data.frame(group = unique(accessor_additive_data$group),
                         probability_group_class_1 = posterior[, 1L],
                         probability_group_class_2 = posterior[, 2L])
  # Probabilities in [0, 1], many near 1e-6 and below: absolute.
  expect_reference(get_results(accessor_additive_fit, "groups"), expected,
                   reference = "independent: Bayes' rule on dense-covariance densities (helper-additive-dense.R)",
                   covers = c("multilpa()", "get_results.multilpa_additive()"),
                   quantity = "group class posteriors, additive", tolerance = 1e-10, measure = "absolute",
                   key = "group")
})

test_that("nobs() of the additive fit is the number of groups", {
  expect_reference(nobs(accessor_additive_fit), length(unique(accessor_additive_data$group)),
                   reference = "independent: number of groups in the data",
                   covers = c("multilpa()", "nobs.multilpa_additive()"),
                   quantity = "number of observations, additive", tolerance = 0)
})

test_that("nobs() of the additive fit is the number of groups lme4 counts", {
  skip_without_reference("lme4")
  # A two-level fit's independent units are its groups: lme4's ngrps() of a
  # random-intercept model of the same ratings.
  reference <- lme4::lmer(y1 ~ 1 + (1 | group), accessor_additive_data, REML = FALSE)
  expect_reference(nobs(accessor_additive_fit), unname(lme4::ngrps(reference)[["group"]]), reference = "lme4",
                   covers = c("multilpa()", "nobs.multilpa_additive()"),
                   quantity = "number of observations (groups), additive", tolerance = 0)
})

# ---- full cross-level family -------------------------------------------------

# 80 teams of six members, two team types with different profile shares (the
# design of test-cross-level.R).
accessor_cross_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  type <- rep(1:2, each = 40L)
  profile_means <- rbind(c(-1.5, -1), c(1.5, 1))
  do.call(rbind, lapply(seq_along(type), function(j) {
    shares <- if (type[j] == 1L) c(0.8, 0.2) else c(0.25, 0.75)
    k <- sample(1:2, 6L, replace = TRUE, prob = shares)
    data.frame(team = sprintf("t%02d", j),
               y1 = stats::rnorm(6L, profile_means[k, 1L], 0.6),
               y2 = stats::rnorm(6L, profile_means[k, 2L], 0.6))
  }))
})
accessor_cross_fit <- quietly(multilpa(accessor_cross_data, c("y1", "y2"), "team", n_profiles = 2,
                                       n_group_classes = 2, family = "full_cross_level",
                                       tol = 1e-10, seed = 1))
accessor_cross_teams <- split(accessor_cross_data[c("y1", "y2")],
                              factor(accessor_cross_data$team, levels = unique(accessor_cross_data$team)))

# The manifest cross-level likelihood (Houle et al., 2026),
#   L_j = sum_h omega_h f(m_j | h) prod_i sum_k pi_{k|h} f(x_ij | k),
# on an unconstrained vector: profile means (profile x indicator, by column),
# log profile variances, group means (class x indicator), log group-mean
# variances, logit of profile 2 within each class, logit of class 2.
accessor_cross_parameters <- function(theta) {
  list(profile_means = matrix(theta[1:4], 2L), profile_variances = matrix(exp(theta[5:8]), 2L),
       group_means = matrix(theta[9:12], 2L), group_variances = matrix(exp(theta[13:16]), 2L),
       composition = cbind(1 - stats::plogis(theta[17:18]), stats::plogis(theta[17:18])),
       weights = c(1 - stats::plogis(theta[[19L]]), stats::plogis(theta[[19L]])))
}
accessor_cross_joint <- function(p) {
  t(vapply(accessor_cross_teams, function(x) {
    x <- as.matrix(x)
    members <- matrix(vapply(1:2, function(k) {
      exp(rowSums(stats::dnorm(x, rep(p$profile_means[k, ], each = nrow(x)),
                               rep(sqrt(p$profile_variances[k, ]), each = nrow(x)), log = TRUE)))
    }, numeric(nrow(x))), nrow(x))
    vapply(1:2, function(h) {
      log(p$weights[h]) +
        sum(stats::dnorm(colMeans(x), p$group_means[h, ], sqrt(p$group_variances[h, ]), log = TRUE)) +
        sum(log(members %*% p$composition[h, ]))
    }, numeric(1))
  }, numeric(2)))
}
accessor_cross_loglik <- function(theta) {
  joint <- accessor_cross_joint(accessor_cross_parameters(theta))
  top <- apply(joint, 1L, max)
  sum(top + log(rowSums(exp(joint - top))))
}
# latents' estimates on that vector, read from as.data.frame() and the composition table.
accessor_cross_start <- local({
  table <- as.data.frame(accessor_cross_fit)
  individual <- subset(table, level == "individual")
  group <- subset(table, level == "group")
  composition <- get_results(accessor_cross_fit, "composition")
  share <- function(h) subset(composition, group_class == h & profile == "profile_2")$share
  c(individual$mean, log(individual$variance), group$mean, log(group$variance),
    stats::qlogis(c(share("group_class_1"), share("group_class_2"))),
    stats::qlogis(subset(group, profile == "group_class_2")$prevalence[[1L]]))
})

test_that("the cross-level estimates are the maximum of the published likelihood, continued by optim()", {
  expect_reference(as.numeric(logLik(accessor_cross_fit)), accessor_cross_loglik(accessor_cross_start),
                   reference = "independent: manifest cross-level likelihood from the published formula (Houle et al., 2026), stats::dnorm()",
                   covers = c("multilpa()", "logLik.multilpa_cross_level()"),
                   quantity = "full cross-level log likelihood at the reported estimates",
                   tolerance = 1e-10, measure = "relative")
  continued <- stats::optim(accessor_cross_start, accessor_cross_loglik, method = "BFGS",
                            control = list(fnscale = -1, reltol = 1e-15, maxit = 1000L))
  expect_identical(continued$convergence, 0L)
  p <- accessor_cross_parameters(continued$par)
  expected <- data.frame(
    level = rep(c("individual", "group"), each = 4L),
    profile = c(rep(c("profile_1", "profile_2"), 2L), rep(c("group_class_1", "group_class_2"), 2L)),
    indicator = rep(rep(c("y1", "y2"), each = 2L), 2L),
    mean = c(as.vector(p$profile_means), as.vector(p$group_means)),
    variance = c(as.vector(p$profile_variances), as.vector(p$group_variances)),
    prevalence = c(rep(as.vector(p$weights %*% p$composition), 2L), rep(p$weights, 2L)))
  # BFGS on numerical gradients stops when the likelihood changes by less
  # than reltol; the estimates then agree to about its square root on a
  # flat optimum: 1e-4 relative.
  note <- "Reference: stats::optim(method = 'BFGS') maximising the published manifest likelihood written in this file, started from latents' estimates."
  expect_reference(as.data.frame(accessor_cross_fit), expected, reference = "stats",
                   covers = c("multilpa()", "as.data.frame.multilpa_cross_level()"),
                   quantity = "cross-level means, variances and prevalences", tolerance = 1e-4,
                   measure = "relative", key = c("level", "profile", "indicator"), note = note)
  expect_reference(summary(accessor_cross_fit)[["profiles"]], expected, reference = "stats",
                   covers = c("multilpa()", "summary.multilpa_cross_level()"),
                   quantity = "cross-level means, variances and prevalences, summary()", tolerance = 1e-4,
                   measure = "relative", key = c("level", "profile", "indicator"), note = note)
  composition <- data.frame(group_class = rep(c("group_class_1", "group_class_2"), 2L),
                            profile = rep(c("profile_1", "profile_2"), each = 2L),
                            share = as.vector(p$composition))
  expect_reference(get_results(accessor_cross_fit, "composition"), composition, reference = "stats",
                   covers = c("multilpa()", "get_results.multilpa_cross_level()"),
                   quantity = "cross-level composition", tolerance = 1e-4, measure = "relative",
                   key = c("group_class", "profile"), note = note)
})

test_that("the cross-level group posteriors are Bayes' rule on the published likelihood", {
  joint <- accessor_cross_joint(accessor_cross_parameters(accessor_cross_start))
  posterior <- exp(joint - apply(joint, 1L, max))
  posterior <- posterior / rowSums(posterior)
  expected <- data.frame(team = names(accessor_cross_teams),
                         probability_group_class_1 = posterior[, 1L],
                         probability_group_class_2 = posterior[, 2L])
  # Probabilities in [0, 1], many below 1e-10: absolute.
  expect_reference(get_results(accessor_cross_fit, "groups"), expected,
                   reference = "independent: manifest cross-level likelihood from the published formula (Houle et al., 2026), stats::dnorm()",
                   covers = c("multilpa()", "get_results.multilpa_cross_level()"),
                   quantity = "cross-level group class posteriors", tolerance = 1e-10, measure = "absolute",
                   key = "team")
})

test_that("nobs() of the cross-level fit is the number of teams", {
  expect_reference(nobs(accessor_cross_fit), length(accessor_cross_teams),
                   reference = "independent: number of groups in the data",
                   covers = c("multilpa()", "nobs.multilpa_cross_level()"),
                   quantity = "number of observations, cross-level", tolerance = 0)
})

test_that("nobs() of the cross-level fit is the number of teams lme4 counts", {
  skip_without_reference("lme4")
  reference <- lme4::lmer(y1 ~ 1 + (1 | team), accessor_cross_data, REML = FALSE)
  expect_reference(nobs(accessor_cross_fit), unname(lme4::ngrps(reference)[["team"]]), reference = "lme4",
                   covers = c("multilpa()", "nobs.multilpa_cross_level()"),
                   quantity = "number of observations (teams), cross-level", tolerance = 0)
})

test_that("parameter_inference() refuses the cross-level families", {
  # No numbers to compare: the working likelihood counts the ratings twice.
  expect_error(parameter_inference(accessor_cross_fit), class = "latents_unsupported_inference")
})
