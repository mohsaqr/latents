# multilpa() and multilca() with groups: profiles within group classes, on a
# simulated two-level set and the bundled course data, across numbers of
# profiles and group classes, both variance models, both covariance models,
# missing data, group weights, centring and scaling by person, and the types
# and order of the identifier. References: the two-level mixture likelihood
# (Vermunt 2003) rebuilt from the result tables with its local-maximum check,
# the single-level fit at one group class, and invariances.

sweep_module_twolevel <- function() {
  simulated <- sweep_two_level()
  course <- latents::course_engagement
  sets <- list(
    simulated = list(data = simulated, vars = c("y1", "y2", "y3"), id = "id"),
    course = list(data = course, vars = c("browse", "lectures", "forum_read", "forum_post"),
                  id = "student"))
  fit <- function(set, k, g, ..., data = set$data, id = set$id, n_starts = 10L) {
    multilpa(data, set$vars, id = id, n_profiles = k, n_group_classes = g, seed = 51L,
             n_starts = n_starts, tol = 1e-10, max_iter = 20000L, ...)
  }
  ll_of <- function(x) sweep_table(x, "information_criteria")$log_likelihood
  par_of <- function(x) sweep_table(x, "information_criteria")$n_parameters

  structure_cases <- lapply(
    sweep_grid(set = names(sets), k = 2:3, g = 1:3, variance = c("varying", "equal"),
               covariance = c("diagonal", "full")), function(option) {
      sweep_define(sprintf("structure[data=%s,k=%d,g=%d,%s,%s]", option$set, option$k, option$g,
                           option$variance, option$covariance), function() {
        set <- sets[[option$set]]
        model <- fit(set, option$k, option$g, variance_model = option$variance,
                     covariance_model = option$covariance)
        p <- length(set$vars)
        dispersion <- (if (identical(option$covariance, "full")) p * (p + 1L) / 2L else p) *
          (if (identical(option$variance, "varying")) option$k else 1L)
        criteria <- sweep_table(model, "information_criteria")
        groups <- length(unique(set$data[[set$id]]))
        rbind(
          sweep_likelihood_checks(model, set$data, group = set$data[[set$id]]),
          sweep_check("parameter count", criteria$n_parameters,
                      (option$g - 1L) + option$g * (option$k - 1L) + option$k * p + dispersion, 0),
          sweep_check("group BIC from the likelihood and parameter count", criteria$bic_groups,
                      -2 * criteria$log_likelihood + criteria$n_parameters * log(groups), 1e-10))
      })
    })

  missing_cases <- lapply(sweep_grid(set = names(sets), k = 2L, g = 1:2,
                                     covariance = c("diagonal", "full")), function(option) {
    sweep_define(sprintf("fiml[data=%s,k=%d,g=%d,%s]", option$set, option$k, option$g,
                         option$covariance), function() {
      set <- sets[[option$set]]
      data <- sweep_punch(set$data, 0.15, columns = set$vars)
      model <- fit(set, option$k, option$g, data = data, missing = "fiml",
                   covariance_model = option$covariance)
      observed <- rowSums(!is.na(data[set$vars])) > 0L
      rbind(sweep_likelihood_checks(model, data, group = data[[set$id]],
                                    n_bic = sum(observed), empty = !observed),
            sweep_check("drop_empty = FALSE keeps the likelihood", ll_of(model),
                        ll_of(fit(set, option$k, option$g, data = data, missing = "fiml",
                                  covariance_model = option$covariance, drop_empty = FALSE)),
                        1e-10))
    })
  })

  settings <- list(list(set = "simulated", k = 2L, g = 2L), list(set = "course", k = 2L, g = 2L),
                   list(set = "simulated", k = 3L, g = 2L))
  invariance <- function(name, code) {
    lapply(settings, function(s) {
      sweep_define(sprintf("%s[data=%s,k=%d,g=%d]", name, s$set, s$k, s$g),
                   function() code(sets[[s$set]], s$k, s$g))
    })
  }
  invariance_cases <- c(
    invariance("seed-reproducible", function(set, k, g) {
      sweep_check("two runs with one seed: every posterior", sweep_posteriors(fit(set, k, g)),
                  sweep_posteriors(fit(set, k, g)), 0)
    }),
    invariance("one-group-class-is-single-level", function(set, k, g) {
      nested <- fit(set, k, 1L)
      flat <- lpa(set$data, set$vars, n_profiles = k, seed = 51L, n_starts = 10L, tol = 1e-10,
                  max_iter = 20000L)
      sweep_check("one group class against lpa(): log likelihood", ll_of(nested), ll_of(flat),
                  1e-8)
    }),
    invariance("identifier-types", function(set, k, g) {
      base <- fit(set, k, g)
      relabel <- function(transform) {
        data <- set$data
        data$key <- transform(data[[set$id]])
        fit(set, k, g, data = data, id = "key")
      }
      rbind(
        sweep_same_solution("character identifier", base,
                            relabel(function(v) paste0("unit", v)), 1e-10),
        sweep_same_solution("factor identifier", base,
                            relabel(function(v) factor(v, levels = rev(unique(v)))), 1e-10),
        sweep_same_solution("double identifier", base, relabel(function(v) v + 0.5), 1e-10))
    }),
    invariance("row-order", function(set, k, g) {
      shuffled <- sweep_seeded(5L, function() set$data[sample(nrow(set$data)), , drop = FALSE])
      sweep_same_solution("rows shuffled across groups", fit(set, k, g),
                          fit(set, k, g, data = shuffled))
    }),
    invariance("cores", function(set, k, g) {
      sweep_same_solution("cores = 2 against cores = 1", fit(set, k, g),
                          fit(set, k, g, cores = 2L), 0)
    }),
    invariance("acceleration", function(set, k, g) {
      sweep_same_solution("acceleration = none against squarem", fit(set, k, g),
                          fit(set, k, g, acceleration = "none"), 1e-6)
    }),
    invariance("restart", function(set, k, g) {
      base <- fit(set, k, g)
      sweep_same_solution("refit from starting_values()", base,
                          fit(set, k, g, start = starting_values(base), n_starts = 1L), 1e-8)
    }),
    invariance("time-is-only-an-order", function(set, k, g) {
      data <- set$data
      data$occasion <- stats::ave(seq_len(nrow(data)), data[[set$id]], FUN = seq_along)
      sweep_same_solution("time given", fit(set, k, g),
                          fit(set, k, g, data = data, time = "occasion"), 0)
    }),
    invariance("unit-weights", function(set, k, g) {
      sweep_same_solution("group weights of one", fit(set, k, g),
                          fit(set, k, g, data = cbind(set$data, w = 1), weights = "w"), 1e-8)
    }),
    invariance("integer-group-weights", function(set, k, g) {
      ids <- unique(set$data[[set$id]])
      w <- sweep_seeded(9L, function() sample(1:3, length(ids), TRUE))
      data <- cbind(set$data, w = w[match(set$data[[set$id]], ids)])
      # Repeat each group as many times as its weight, each copy a group of its own.
      copies <- do.call(rbind, lapply(1:3, function(copy) {
        rows <- data[data$w >= copy, , drop = FALSE]
        rows$key <- paste(rows[[set$id]], copy)
        rows
      }))
      weighted <- fit(set, k, g, data = data, weights = "w")
      # From the weighted estimate, so that both fits climb the same hill.
      expanded <- fit(set, k, g, data = copies, id = "key", start = starting_values(weighted),
                      n_starts = 1L)
      rbind(
        sweep_check("integer group weights against repeated groups: pseudo log likelihood",
                    ll_of(weighted), length(ids) / sum(w) * ll_of(expanded), 1e-6),
        sweep_check("integer group weights against repeated groups: sorted means",
                    sort(sweep_parameters(weighted)$means), sort(sweep_parameters(expanded)$means),
                    1e-3),
        sweep_likelihood_checks(weighted, data, group = data[[set$id]],
                                weights = data$w * length(ids) / sum(w)))
    }),
    invariance("person-centering", function(set, k, g) {
      data <- set$data
      centred <- data
      centred[set$vars] <- lapply(set$vars, function(v) {
        data[[v]] - stats::ave(data[[v]], data[[set$id]])
      })
      scaled <- centred
      scaled[set$vars] <- lapply(set$vars, function(v) {
        centred[[v]] / stats::ave(data[[v]], data[[set$id]], FUN = stats::sd)
      })
      rbind(
        sweep_same_solution("centering = person against person-centred data",
                            fit(set, k, g, centering = "person"),
                            fit(set, k, g, data = centred), 1e-6),
        sweep_same_solution("centering and scaling = person against within-person z-scores",
                            fit(set, k, g, centering = "person", scaling = "person"),
                            fit(set, k, g, data = scaled), 1e-6))
    }))

  categorical <- sweep_two_level_categorical()
  items <- paste0("q", 1:4)
  fit_lca <- function(k, g, ..., data = categorical) {
    multilca(data, items, id = "id", n_profiles = k, n_group_classes = g, seed = 61L,
             n_starts = 10L, tol = 1e-10, max_iter = 20000L, ...)
  }
  categorical_cases <- c(
    lapply(sweep_grid(k = 2:3, g = 1:3, missing = c("complete", "fiml")), function(option) {
      sweep_define(sprintf("multilca[k=%d,g=%d,%s]", option$k, option$g, option$missing),
                   function() {
        holes <- identical(option$missing, "fiml")
        data <- if (holes) sweep_punch(categorical, 0.12, columns = items) else categorical
        model <- if (holes) fit_lca(option$k, option$g, data = data, missing = "fiml") else
          fit_lca(option$k, option$g)
        observed <- rowSums(!is.na(data[items])) > 0L
        rbind(sweep_likelihood_checks(model, data, group = data$id, n_bic = sum(observed),
                                      empty = !observed),
              sweep_check("parameter count", par_of(model),
                          (option$g - 1L) + option$g * (option$k - 1L) +
                            option$k * length(items), 0))
      })
    }),
    list(sweep_define("multilca-is-multilpa-with-categorical", function() {
      direct <- multilpa(categorical, items, id = "id", n_profiles = 2L, n_group_classes = 2L,
                         categorical = items, seed = 61L, n_starts = 10L, tol = 1e-10,
                         max_iter = 20000L)
      sweep_same_solution("multilca() against multilpa(categorical = vars)", fit_lca(2L, 2L),
                          direct, 1e-10)
    })))

  refusal_cases <- list(sweep_define("refusals", function() {
    set <- sets$simulated
    data <- set$data
    rbind(
      sweep_refuses("an unknown id", function() multilpa(data, set$vars, id = "nope", 2L, 2L)),
      sweep_refuses("a missing id",
                    function() multilpa(transform(data, id = replace(id, 3L, NA)), set$vars,
                                        id = "id", 2L, 2L)),
      sweep_refuses("the id among the indicators",
                    function() multilpa(data, c(set$vars, "id"), id = "id", 2L, 2L)),
      sweep_refuses("n_group_classes = 0",
                    function() multilpa(data, set$vars, id = "id", 2L, n_group_classes = 0L)),
      sweep_refuses("more group classes than groups",
                    function() multilpa(data[data$id <= 2L, ], set$vars, id = "id", 2L, 3L)),
      sweep_refuses("group classes with one profile",
                    function() multilpa(data, set$vars, id = "id", 1L, 2L)),
      sweep_refuses("group classes when every group has one row",
                    function() multilpa(transform(data, id = seq_along(id)), set$vars, id = "id",
                                        2L, 2L)),
      sweep_refuses("weights varying within a group",
                    function() multilpa(transform(data, w = seq_along(id)), set$vars, id = "id",
                                        2L, 2L, weights = "w")),
      sweep_refuses("time repeated within a group",
                    function() multilpa(transform(data, t = 1L), set$vars, id = "id", 2L, 2L,
                                        time = "t")),
      sweep_refuses("between_variance with the profiles family",
                    function() multilpa(data, set$vars, id = "id", 2L, 2L,
                                        between_variance = "equal")),
      sweep_refuses("an unknown family",
                    function() multilpa(data, set$vars, id = "id", 2L, 2L, family = "other")),
      sweep_refuses("an unknown covariate",
                    function() multilpa(data, set$vars, id = "id", 2L, 2L,
                                        profile_covariates = "nope")))
  }))

  c(structure_cases, missing_cases, invariance_cases, categorical_cases, refusal_cases)
}
