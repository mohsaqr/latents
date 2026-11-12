# lca() and the categorical arguments of multilpa(): one to three classes on
# four datasets, missing data, weights, local dependence, the probability
# floor and the input types. References: poLCA (its EM continued from the
# package's estimates, and its own search), observed proportions at one class,
# and the likelihood rebuilt from the result tables.

sweep_module_categorical <- function() {
  datasets <- sweep_categorical()
  fit_lca <- function(data, k, ..., vars = names(data), n_starts = 10L) {
    lca(data, vars, n_classes = k, seed = 21L, n_starts = n_starts, tol = 1e-10,
        max_iter = 20000L, ...)
  }
  as_codes <- function(data) as.data.frame(lapply(data, function(v) as.integer(factor(v))))
  item_levels <- function(data) lapply(data, function(v) levels(factor(v)))
  formula_of <- function(vars) {
    stats::as.formula(sprintf("cbind(%s) ~ 1", paste(vars, collapse = ", ")))
  }
  #' The fit's response probabilities as poLCA's list of class-by-category matrices
  probabilities_of <- function(fit, data) {
    table <- sweep_table(fit, "responses")
    levels_of <- item_levels(data)
    profile_levels <- sort(unique(table$profile))
    stats::setNames(lapply(names(data), function(v) {
      unname(sweep_wide(table[table$indicator == v, , drop = FALSE], "profile", "category",
                        "probability", profile_levels, levels_of[[v]]))
    }), names(data))
  }
  polca <- function(data, k, ...) {
    suppressWarnings(poLCA::poLCA(formula_of(names(data)), as_codes(data), nclass = k,
                                  maxiter = 50000L, tol = 1e-12, na.rm = FALSE,
                                  verbose = FALSE, calc.se = FALSE, ...))
  }
  #' Checks of one fit against poLCA and closed forms
  reference_checks <- function(fit, data, k) {
    criteria <- sweep_table(fit, "information_criteria")
    categories <- vapply(item_levels(data), length, integer(1))
    own <- probabilities_of(fit, data)
    rows <- sweep_check("parameter count", criteria$n_parameters,
                        (k - 1L) + k * sum(categories - 1L), 0)
    searched <- sweep_seeded(31L, function() polca(data, k, nrep = if (k == 1L) 1L else 15L))
    rows <- rbind(rows, sweep_at_least("optimum at least that of poLCA's search",
                                       criteria$log_likelihood, searched$llik, 1e-6))
    if (k == 1L) {
      proportions <- lapply(data, function(v) as.numeric(prop.table(table(factor(v)))))
      return(rbind(rows, sweep_check("one class: probabilities are the observed proportions",
                                     unlist(own), unlist(proportions), 1e-8),
                   sweep_check("one class: poLCA's log likelihood", criteria$log_likelihood,
                               searched$llik, 1e-8)))
    }
    continued <- polca(data, k, probs.start = own)
    rbind(rows,
          sweep_check("poLCA's EM from the package's probabilities: log likelihood",
                      criteria$log_likelihood, continued$llik, 1e-6),
          sweep_check("poLCA's EM from the package's probabilities: probabilities",
                      unlist(own), unlist(lapply(continued$probs, unname)), 2e-3),
          sweep_check("poLCA's EM from the package's probabilities: shares",
                      sweep_parameters(fit)$mixing, continued$P, 2e-3))
  }

  structure_cases <- lapply(
    sweep_grid(dataset = names(datasets), k = 1:3, missing = c("complete", "fiml")),
    function(option) {
      sweep_define(sprintf("classes[data=%s,k=%d,%s]", option$dataset, option$k, option$missing),
                   function() {
        holes <- identical(option$missing, "fiml")
        data <- if (holes) sweep_punch(datasets[[option$dataset]], 0.1) else datasets[[option$dataset]]
        fit <- if (holes) fit_lca(data, option$k, missing = "fiml") else fit_lca(data, option$k)
        observed <- rowSums(!is.na(data)) > 0L
        rbind(sweep_likelihood_checks(fit, data, n_bic = sum(observed), empty = !observed),
              reference_checks(fit, data, option$k))
      })
    })

  settings <- list(list(dataset = "carcinoma", k = 2L), list(dataset = "values", k = 2L),
                   list(dataset = "nominal", k = 3L))
  invariance <- function(name, code) {
    lapply(settings, function(s) {
      sweep_define(sprintf("%s[data=%s,k=%d]", name, s$dataset, s$k),
                   function() code(datasets[[s$dataset]], s$k))
    })
  }
  invariance_cases <- c(
    invariance("seed-reproducible", function(data, k) {
      sweep_check("two runs with one seed: every posterior", sweep_posteriors(fit_lca(data, k)),
                  sweep_posteriors(fit_lca(data, k)), 0)
    }),
    invariance("wrapper", function(data, k) {
      base <- fit_lca(data, k)
      direct <- multilpa(data, names(data), id = NULL, n_profiles = k, n_group_classes = 1L,
                         categorical = names(data), seed = 21L, n_starts = 10L, tol = 1e-10,
                         max_iter = 20000L)
      profiles <- lpa(data, names(data), n_profiles = k, categorical = names(data), seed = 21L,
                      n_starts = 10L, tol = 1e-10, max_iter = 20000L)
      grid <- enumerate_lca(data, names(data), n_classes = k, seed = 21L, n_starts = 10L,
                            tol = 1e-10, max_iter = 20000L)
      rbind(sweep_same_solution("lca() against multilpa(categorical = vars)", base, direct, 1e-10),
            sweep_same_solution("lca() against lpa(categorical = vars)", base, profiles, 1e-10),
            sweep_same_solution("lca() against the enumerate_lca() candidate", base,
                                get_candidate(grid, n_profiles = k), 1e-10))
    }),
    invariance("row-order", function(data, k) {
      shuffled <- sweep_seeded(5L, function() data[sample(nrow(data)), , drop = FALSE])
      sweep_same_solution("rows shuffled", fit_lca(data, k), fit_lca(shuffled, k))
    }),
    invariance("column-order", function(data, k) {
      sweep_same_solution("indicators reversed", fit_lca(data, k),
                          fit_lca(data, k, vars = rev(names(data))))
    }),
    invariance("cores", function(data, k) {
      sweep_same_solution("cores = 2 against cores = 1", fit_lca(data, k),
                          fit_lca(data, k, cores = 2L), 0)
    }),
    invariance("acceleration", function(data, k) {
      sweep_same_solution("acceleration = none against squarem", fit_lca(data, k),
                          fit_lca(data, k, acceleration = "none"), 1e-6)
    }),
    invariance("search", function(data, k) {
      base <- fit_lca(data, k, n_starts = 20L)
      rbind(sweep_same_solution("search = kmeans", base,
                                fit_lca(data, k, n_starts = 20L, search = "kmeans"), 1e-6),
            sweep_same_solution("search = emEM", base,
                                fit_lca(data, k, n_starts = 20L, search = "emEM"), 1e-6))
    }),
    invariance("select_start", function(data, k) {
      sweep_same_solution("select_start = converged", fit_lca(data, k),
                          fit_lca(data, k, select_start = "converged"), 1e-8)
    }),
    invariance("restart", function(data, k) {
      base <- fit_lca(data, k)
      again <- fit_lca(data, k, start = starting_values(base), n_starts = 1L)
      rbind(sweep_same_solution("refit from starting_values()", base, again, 1e-8),
            sweep_true("the refit needs few iterations",
                       sweep_table(again, "model")$iterations <= 5L,
                       sprintf("%d iterations", sweep_table(again, "model")$iterations)))
    }),
    invariance("fiml-on-complete-data", function(data, k) {
      sweep_same_solution("missing = fiml on complete data", fit_lca(data, k),
                          fit_lca(data, k, missing = "fiml"), 1e-8)
    }),
    invariance("min_probability", function(data, k) {
      # The floor moves a bounded probability by its own size, so the likelihood
      # moves by at most rows x items x floor.
      sweep_same_solution("min_probability = 1e-8 against 1e-10", fit_lca(data, k),
                          fit_lca(data, k, min_probability = 1e-8), 1e-5)
    }),
    invariance("unit-weights", function(data, k) {
      sweep_same_solution("weights of one", fit_lca(data, k),
                          fit_lca(cbind(data, w = 1), k, vars = names(data), weights = "w"), 1e-8)
    }),
    invariance("integer-weights", function(data, k) {
      w <- sweep_seeded(9L, function() sample(1:3, nrow(data), TRUE))
      weighted <- fit_lca(cbind(data, w = w), k, vars = names(data), weights = "w")
      expanded <- fit_lca(data[rep(seq_len(nrow(data)), w), , drop = FALSE], k,
                          start = starting_values(weighted), n_starts = 1L)
      rbind(
        sweep_check("integer weights against repeated rows: pseudo log likelihood",
                    sweep_table(weighted, "information_criteria")$log_likelihood,
                    nrow(data) / sum(w) *
                      sweep_table(expanded, "information_criteria")$log_likelihood, 1e-6),
        sweep_check("integer weights against repeated rows: sorted probabilities",
                    sort(sweep_table(weighted, "responses")$probability),
                    sort(sweep_table(expanded, "responses")$probability), 1e-3))
    }),
    invariance("pattern-weights", function(data, k) {
      # The response patterns with their frequencies as weights are the same data.
      key <- do.call(paste, c(data, sep = "\r"))
      patterns <- data[!duplicated(key), , drop = FALSE]
      frequency <- as.numeric(table(key)[key[!duplicated(key)]])
      weighted <- fit_lca(cbind(patterns, w = frequency), k, vars = names(data), weights = "w")
      base <- fit_lca(data, k)
      rbind(
        sweep_check("pattern frequencies as weights: pseudo log likelihood",
                    sweep_table(weighted, "information_criteria")$log_likelihood,
                    nrow(patterns) / nrow(data) *
                      sweep_table(base, "information_criteria")$log_likelihood, 1e-6),
        sweep_check("pattern frequencies as weights: sorted probabilities",
                    sort(sweep_table(weighted, "responses")$probability),
                    sort(sweep_table(base, "responses")$probability), 1e-3))
    }),
    invariance("input-types", function(data, k) {
      base <- fit_lca(data, k)
      rbind(
        sweep_same_solution("character columns", base,
                            fit_lca(as.data.frame(lapply(data, as.character)), k), 1e-10),
        sweep_same_solution("integer codes", base, fit_lca(as_codes(data), k), 1e-10),
        sweep_same_solution("ordered factors", base,
                            fit_lca(as.data.frame(lapply(data, function(v) factor(v, ordered = TRUE))),
                                    k), 1e-10),
        sweep_same_solution("an unused factor level", base,
                            fit_lca(as.data.frame(lapply(data, function(v) {
                              factor(v, levels = c(levels(factor(v)), "never"))
                            })), k), 1e-10))
    }),
    invariance("local-dependence", function(data, k) {
      pair <- names(data)[1:2]
      joined <- cbind(joint = interaction(data[[pair[[1L]]]], data[[pair[[2L]]]], drop = TRUE),
                      data[setdiff(names(data), pair)])
      dependent <- fit_lca(data, k, local_dependence = paste(pair, collapse = ":"))
      combined <- fit_lca(joined, k)
      complete_table <- nlevels(joined$joint) ==
        nlevels(factor(data[[pair[[1L]]]])) * nlevels(factor(data[[pair[[2L]]]]))
      rbind(
        sweep_check("local dependence against one item of the pair's combinations: log likelihood",
                    sweep_table(dependent, "information_criteria")$log_likelihood,
                    sweep_table(combined, "information_criteria")$log_likelihood, 1e-6),
        if (complete_table) {
          sweep_check("local dependence against one item of the pair's combinations: parameters",
                      sweep_table(dependent, "information_criteria")$n_parameters,
                      sweep_table(combined, "information_criteria")$n_parameters, 0)
        } else {
          sweep_no_reference("local dependence: parameters",
                             "an unobserved combination counts under local_dependence only")
        })
    }))

  logical_case <- list(sweep_define("input-types[logical]", function() {
    data <- datasets$carcinoma
    flags <- as.data.frame(lapply(data, function(v) as.integer(v) == 2L))
    sweep_same_solution("logical columns", fit_lca(data, 2L), fit_lca(flags, 2L), 1e-10)
  }))

  refusal_cases <- list(sweep_define("refusals", function() {
    data <- datasets$carcinoma
    vars <- names(data)
    rbind(
      sweep_refuses("n_classes = 0", function() lca(data, vars, n_classes = 0L)),
      sweep_refuses("n_classes = 2.5", function() lca(data, vars, n_classes = 2.5)),
      sweep_refuses("an unknown indicator", function() lca(data, c(vars, "nope"), 2L)),
      sweep_refuses("missing values without fiml",
                    function() lca(sweep_punch(data, 0.1), vars, 2L)),
      sweep_refuses("a single-category indicator",
                    function() lca(cbind(data, one = factor("x")), c(vars, "one"), 2L)),
      sweep_warns("a continuous column as an item warns",
                  function() lca(cbind(data, x = stats::rnorm(nrow(data))), c(vars, "x"), 2L,
                                 n_starts = 2L, seed = 1L),
                  "latents_suspect_categorical"),
      sweep_refuses("min_probability = 0.6", function() lca(data, vars, 2L, min_probability = 0.6)),
      sweep_refuses("min_probability = -1", function() lca(data, vars, 2L, min_probability = -1)),
      sweep_refuses("id passed to lca()", function() lca(data, vars, 2L, id = "A")),
      sweep_refuses("categorical naming a column outside vars",
                    function() multilpa(data, vars[1:3], id = NULL, n_profiles = 2L,
                                        n_group_classes = 1L, categorical = vars)),
      sweep_refuses("local dependence on an unknown item",
                    function() lca(data, vars, 2L, local_dependence = "A:nope")),
      sweep_refuses("local dependence of an item with itself",
                    function() lca(data, vars, 2L, local_dependence = "A:A")),
      sweep_refuses("more classes than rows", function() lca(data[1:2, ], vars, 3L)))
  }))

  c(structure_cases, invariance_cases, logical_case, refusal_cases)
}
