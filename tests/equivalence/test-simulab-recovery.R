# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Recovery of known parameters: simulab draws data from parameters set here,
# latents fits them, and each estimate must lie within four of its standard
# errors of the parameter that generated the data. Several seeds per design;
# large samples, so the estimators' finite-sample bias is far inside the bound.
# The standard errors are latents' own (parameter_inference(), observed
# information; vcov() for mixture_regression()), which the other files check
# against numDeriv, Mplus and Latent GOLD. Four standard errors is a
# two-sided test at about 6e-5 per parameter, fixed before the first run.
# The fitted classes are matched to the generating ones by the order that
# minimises the squared distance of the estimates from the truth. Of the two
# probabilities of a binary item only the first is compared (the second is one
# minus it).
#
# This is evidence about the truth, not about other software: authority
# `independent` (simulab is the owner's package; the truth is the parameters).

recovery_seeds <- c(101L, 102L, 103L)

all_orders <- function(k) {
  if (k == 1L) return(list(1L))
  do.call(c, lapply(seq_len(k), \(first) lapply(all_orders(k - 1L), \(rest) c(first, setdiff(seq_len(k), first)[rest]))))
}

# Standard errors of the natural parameters, named as coef() names them.
inference_errors <- function(fit) {
  table <- parameter_inference(fit)
  # A group-class probability has no term: coef() names it level.parameter.outcome.
  names <- ifelse(is.na(table$term), sprintf("%s.%s.%s", table$level, table$parameter, table$outcome),
                  sprintf("%s.%s.%s.%s", table$level, table$parameter, table$outcome, table$term))
  stats::setNames(table$standard_error, names)
}

# `truth_for(order)`: the generating parameters named by the fitted labels,
# where order[g] is the fitted label of generating class g (a list of such
# vectors for two kinds of class). Returns the truth of the best order.
best_truth <- function(estimates, truth_for, orders) {
  candidates <- lapply(orders, truth_for)
  distance <- vapply(candidates, \(truth) sum((estimates[names(truth)] - truth)^2), numeric(1))
  candidates[[which.min(distance)]]
}

expect_recovered <- function(estimates, errors, truth, covers, design, seed) {
  stopifnot("every compared parameter has a standard error" = all(names(truth) %in% names(errors)),
            "every compared parameter has an estimate" = all(names(truth) %in% names(estimates)))
  invisible(lapply(names(truth), \(name) {
    expect_reference(estimates[[name]], truth[[name]],
                     reference = sprintf("independent: recovery of the parameters of simulab %s", design$generator),
                     covers = covers, quantity = sprintf("%s, %s, seed %d", design$label, name, seed),
                     tolerance = 4 * errors[[name]], measure = "absolute",
                     note = "Four standard errors of the estimate (latents' observed information), fixed before the run.")
  }))
}

# ---- lpa(): three profiles, four indicators, unit variances (EEI) ----------

lpa_design <- list(generator = "simulate_lpa()", label = "lpa(), 3 profiles EEI, n = 2000", n = 2000L,
                   proportions = c(0.4, 0.35, 0.25),
                   means = rbind(c(0, 0, 0, 0), c(4, 4, 0, 0), c(0, 4, 4, 4)))

test_that("lpa() recovers the parameters of simulab::simulate_lpa()", {
  design <- lpa_design
  vars <- sprintf("indicator_%d", 1:4)
  invisible(lapply(recovery_seeds, \(seed) {
    data <- as.data.frame(simulab::simulate_lpa(design$n, means = design$means, sds = 1,
                                                proportions = design$proportions, seed = seed))
    fit <- lpa(data, vars, n_profiles = 3, model = "EEI", n_starts = 10, seed = 1)
    estimates <- coef(fit)
    truth_for <- \(order) c(
      unlist(lapply(1:3, \(g) stats::setNames(design$means[g, ],
                                              sprintf("measurement.mean.profile_%d.%s", order[g], vars)))),
      stats::setNames(design$proportions, sprintf("profile.probability.profile_%d.group_class_1", order)),
      stats::setNames(rep(1, 4), sprintf("measurement.variance.shared.%s", vars)))
    expect_recovered(estimates, inference_errors(fit), best_truth(estimates, truth_for, all_orders(3L)),
                     c("lpa()", "coef.multilpa()"), design, seed)
  }))
})

# ---- lca(): three classes, six binary items ---------------------------------

lca_design <- list(generator = "simulate_lca()", label = "lca(), 3 classes, 6 binary items, n = 3000", n = 3000L,
                   proportions = c(0.5, 0.3, 0.2),
                   yes = rbind(rep(0.9, 6L), c(0.9, 0.9, 0.9, 0.1, 0.1, 0.1), rep(0.1, 6L)))

test_that("lca() recovers the parameters of simulab::simulate_lca()", {
  design <- lca_design
  probabilities <- array(0, c(3L, 6L, 2L))
  probabilities[, , 1L] <- design$yes
  probabilities[, , 2L] <- 1 - design$yes
  items <- sprintf("item_%d", 1:6)
  invisible(lapply(recovery_seeds, \(seed) {
    data <- as.data.frame(simulab::simulate_lca(design$n, probabilities = probabilities,
                                                proportions = design$proportions, seed = seed))
    fit <- lca(data, items, n_classes = 3, n_starts = 10, seed = 1)
    estimates <- coef(fit)
    truth_for <- \(order) c(
      unlist(lapply(1:3, \(g) stats::setNames(design$yes[g, ],
                                              sprintf("measurement.response.profile_%d.%s:1", order[g], items)))),
      stats::setNames(design$proportions, sprintf("profile.probability.profile_%d.group_class_1", order)))
    expect_recovered(estimates, inference_errors(fit), best_truth(estimates, truth_for, all_orders(3L)),
                     c("lca()", "coef.multilpa()"), design, seed)
  }))
})

# ---- lta(): a two-state hidden Markov chain, a three-category observation --

lta_design <- list(generator = "simulate_hmm()", label = "lta(), 2 states, 5 occasions, n = 2000", n = 2000L,
                   occasions = 5L, initial = c(0.6, 0.4),
                   transition = rbind(c(0.8, 0.2), c(0.3, 0.7)),
                   emission = rbind(c(0.8, 0.1, 0.1), c(0.1, 0.2, 0.7)))

test_that("lta() recovers the parameters of simulab::simulate_hmm()", {
  design <- lta_design
  categories <- sprintf("Observation %d", 1:3)
  invisible(lapply(recovery_seeds, \(seed) {
    data <- as.data.frame(simulab::simulate_hmm(design$n, transition = design$transition,
                                                chain_length = design$occasions, emission = design$emission,
                                                initial = design$initial, seed = seed))
    data$observation <- factor(data$observation, levels = categories)
    fit <- lta(data, "observation", id = "id", n_profiles = 2, time = "occasion", categorical = "observation",
               n_starts = 10, seed = 1)
    estimates <- coef(fit)
    truth_for <- \(order) c(
      stats::setNames(design$initial, sprintf("profile.initial_probability.profile_%d.group_class_1", order)),
      unlist(lapply(1:2, \(from) stats::setNames(
        design$transition[from, ],
        sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d", order, order[from])))),
      unlist(lapply(1:2, \(g) stats::setNames(
        design$emission[g, ], sprintf("measurement.response.profile_%d.observation:%s", order[g], categories)))))
    expect_recovered(estimates, inference_errors(fit), best_truth(estimates, truth_for, all_orders(2L)),
                     c("lta()", "coef.multilpa_transitions()"), design, seed)
  }))
})

# ---- multilpa(): 100 clusters of 20, two cluster classes, three profiles ---

multilpa_design <- list(generator = "simulate_ml_lpa()",
                        label = "multilpa(), 2 cluster classes, 3 profiles, 100 clusters of 20",
                        clusters = 100L, cluster_size = 20L, class_proportions = c(0.5, 0.5),
                        means = rbind(c(0, 0, 0), c(3, 3, 0), c(0, 3, 3)),
                        profile_probabilities = rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6)))

test_that("multilpa() recovers the parameters of simulab::simulate_ml_lpa()", {
  design <- multilpa_design
  # Every pairing of a profile order with a cluster-class order.
  orders <- do.call(c, lapply(all_orders(3L), \(profiles) lapply(all_orders(2L), \(classes) {
    list(profiles = profiles, classes = classes)
  })))
  invisible(lapply(recovery_seeds, \(seed) {
    data <- as.data.frame(simulab::simulate_ml_lpa(
      clusters = design$clusters, cluster_size = design$cluster_size, means = design$means,
      profile_probabilities = design$profile_probabilities,
      cluster_class_proportions = design$class_proportions, seed = seed))
    vars <- sprintf("indicator_%d", 1:3)
    fit <- multilpa(data, vars, id = "cluster", n_profiles = 3, n_group_classes = 2, variance_model = "equal",
                    n_starts = 10, seed = 1)
    estimates <- coef(fit)
    truth_for <- \(order) {
      p <- order$profiles
      h <- order$classes
      c(unlist(lapply(1:3, \(g) stats::setNames(design$means[g, ], sprintf("measurement.mean.profile_%d.%s", p[g], vars)))),
        stats::setNames(rep(1, 3), sprintf("measurement.variance.shared.%s", vars)),
        stats::setNames(design$class_proportions, sprintf("group.probability.group_class_%d", h)),
        unlist(lapply(1:2, \(c) stats::setNames(design$profile_probabilities[c, ],
                                                sprintf("profile.probability.profile_%d.group_class_%d", p, h[c])))))
    }
    expect_recovered(estimates, inference_errors(fit), best_truth(estimates, truth_for, orders),
                     c("multilpa()", "coef.multilpa()"), design, seed)
  }))
})

# ---- mixture_regression(): two regression lines -----------------------------

regression_design <- list(generator = "simulate_regression()",
                          label = "mixture_regression(), 2 classes of 1200 and 800 rows",
                          sizes = c(1200L, 800L), intercepts = c(1, 4), slopes = c(2, -1), sds = c(0.5, 0.7))

test_that("mixture_regression() recovers the parameters of simulab::simulate_regression()", {
  design <- regression_design
  invisible(lapply(recovery_seeds, \(seed) {
    # simulab draws each class's rows from its own line; the classes are stacked.
    data <- do.call(rbind, lapply(1:2, \(g) {
      drawn <- as.data.frame(simulab::simulate_regression(
        design$sizes[g], coefficients = c("(Intercept)" = design$intercepts[g], x = design$slopes[g]),
        error_sd = design$sds[g], seed = seed * 10L + g))
      drawn[c("x", "outcome")]
    }))
    fit <- mixture_regression(outcome ~ x, data, n_classes = 2, n_starts = 5, seed = 1)
    estimates <- coef(fit)
    errors <- sqrt(diag(vcov(fit)))
    share <- design$sizes / sum(design$sizes)
    truth_for <- \(order) c(
      stats::setNames(design$intercepts, sprintf("coefficient.class_%d.(Intercept)", order)),
      stats::setNames(design$slopes, sprintf("coefficient.class_%d.x", order)),
      stats::setNames(log(design$sds), sprintf("log_sigma.class_%d", order)),
      "membership.class_2.(Intercept)" = log(share[match(2L, order)] / share[match(1L, order)]))
    expect_recovered(estimates, errors, best_truth(estimates, truth_for, all_orders(2L)),
                     c("mixture_regression()", "coef.latents_mixture_regression()"), design, seed)
  }))
})
