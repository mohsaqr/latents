# lpa() and the Gaussian arguments of multilpa(): the fourteen covariance
# structures by every spelling, missing data, weights, centring and scaling,
# the search and optimizer controls, and the input types, across four datasets
# and one to three profiles. References: mclust (its M-steps continued from
# the package's posteriors, its one-component fits, its parameter counts and
# its own search), closed forms, and the likelihood rebuilt from the tables.

#' The reference's warnings are not the package's: keep them out of the record
sweep_reference <- function(code) suppressWarnings(code)

sweep_module_gaussian <- function() {
  # mclust's fitting functions look up `mclustBIC` on the search path.
  suppressPackageStartupMessages(library(mclust))
  datasets <- sweep_continuous()
  models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE",
              "EEV", "VEV", "EVV", "VVV")
  fit_lpa <- function(data, k, ..., vars = names(data), n_starts = 10L) {
    lpa(data, vars, n_profiles = k, seed = 11L, n_starts = n_starts, tol = 1e-10,
        max_iter = 20000L, ...)
  }
  components <- function(model) {
    letters <- strsplit(model, "", fixed = TRUE)[[1L]]
    word <- c(E = "equal", V = "varying", I = "spherical")
    if (identical(letters[[2L]], "I")) {
      return(list(volume = word[[letters[[1L]]]], shape = "spherical"))
    }
    list(volume = word[[letters[[1L]]]], shape = word[[letters[[2L]]]],
         orientation = if (identical(letters[[3L]], "I")) "axis" else word[[letters[[3L]]]])
  }

  # ---- every structure, on every dataset, at one to three profiles ----------
  structure_cases <- lapply(
    sweep_grid(dataset = names(datasets), model = models, k = 1:3), function(option) {
      sweep_define(sprintf("structure[data=%s,model=%s,k=%d]", option$dataset, option$model,
                           option$k), function() {
        data <- datasets[[option$dataset]]
        x <- as.matrix(data)
        fit <- fit_lpa(data, option$k, model = option$model)
        criteria <- sweep_table(fit, "information_criteria")
        parameters <- sweep_parameters(fit)
        own <- rbind(
          sweep_likelihood_checks(fit, data),
          sweep_check("parameter count is mclust's", criteria$n_parameters,
                      mclust::nMclustParams(option$model, ncol(x), option$k), 0))
        searched <- sweep_reference(mclust::Mclust(x, G = option$k, modelNames = option$model,
                                                   verbose = FALSE))
        search_check <- if (is.null(searched) || !is.finite(searched$loglik)) {
          sweep_no_reference("optimum at least that of mclust's search", "mclust found no fit")
        } else {
          sweep_at_least("optimum at least that of mclust's search", criteria$log_likelihood,
                         searched$loglik, 1e-6)
        }
        if (option$k == 1L) return(rbind(own, search_check))
        continued <- sweep_reference(mclust::me(
          data = x, modelName = option$model, z = sweep_posteriors(fit),
          control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(20000L, 20000L))))
        if (!is.finite(continued$loglik)) {
          return(rbind(own, search_check, sweep_no_reference(
            "mclust's EM continued from the package's posteriors", "mclust stopped without a fit")))
        }
        # mclust's VVE M-step is inexact: continued from the package's solution it
        # lowers the likelihood (by 4e-4 to 1.2 on these data), and the package
        # restarted from where mclust stops returns to its own value
        # (diagnostics/vve.R). So under VVE mclust may only fail to improve.
        if (identical(option$model, "VVE")) {
          return(rbind(own, search_check, sweep_at_least(
            "mclust's EM from the package's posteriors does not improve the likelihood",
            criteria$log_likelihood, continued$loglik, 1e-8)))
        }
        rbind(own, search_check,
              sweep_check("mclust's EM from the package's posteriors: log likelihood",
                          criteria$log_likelihood, continued$loglik, 1e-6),
              sweep_check("mclust's EM from the package's posteriors: means",
                          parameters$means, t(continued$parameters$mean), 1e-3),
              sweep_check("mclust's EM from the package's posteriors: shares",
                          parameters$mixing, continued$parameters$pro, 1e-3),
              sweep_check("mclust's EM from the package's posteriors: covariances",
                          unlist(parameters$covariances),
                          continued$parameters$variance$sigma, 1e-3))
      })
    })

  # ---- the same structure by each spelling -----------------------------------
  spelling_cases <- lapply(models, function(model) {
    sweep_define(sprintf("spelling[model=%s]", model), function() {
      data <- datasets$blobs
      named <- fit_lpa(data, 2L, model = model)
      spelled <- do.call(fit_lpa, c(list(data, 2L), components(model)))
      rows <- sweep_same_solution("model code against volume/shape/orientation", named, spelled)
      pair <- list(EEI = c("equal", "diagonal"), VVI = c("varying", "diagonal"),
                   EEE = c("equal", "full"), VVV = c("varying", "full"))[[model]]
      if (is.null(pair)) return(rows)
      rbind(rows, sweep_same_solution(
        "model code against variance_model/covariance_model", named,
        fit_lpa(data, 2L, variance_model = pair[[1L]], covariance_model = pair[[2L]])))
    })
  })

  # ---- invariances that hold for any correct estimator ------------------------
  settings <- list(list(dataset = "srl", k = 2L, model = "VVI"),
                   list(dataset = "iris", k = 3L, model = "EEE"),
                   list(dataset = "faithful", k = 2L, model = "VVV"),
                   list(dataset = "blobs", k = 3L, model = "EEI"),
                   list(dataset = "blobs", k = 3L, model = "VEV"))
  invariance <- function(name, code) {
    lapply(settings, function(s) {
      sweep_define(sprintf("%s[data=%s,model=%s,k=%d]", name, s$dataset, s$model, s$k),
                   function() code(datasets[[s$dataset]], s$k, s$model))
    })
  }
  invariance_cases <- c(
    invariance("seed-reproducible", function(data, k, model) {
      first <- fit_lpa(data, k, model = model)
      second <- fit_lpa(data, k, model = model)
      sweep_check("two runs with one seed: every posterior", sweep_posteriors(first),
                  sweep_posteriors(second), 0)
    }),
    invariance("wrapper", function(data, k, model) {
      base <- fit_lpa(data, k, model = model)
      direct <- multilpa(data, names(data), id = NULL, n_profiles = k, n_group_classes = 1L,
                         model = model, seed = 11L, n_starts = 10L, tol = 1e-10,
                         max_iter = 20000L)
      grid <- enumerate_lpa(data, names(data), n_profiles = k, model = model, seed = 11L,
                            n_starts = 10L, tol = 1e-10, max_iter = 20000L)
      rbind(sweep_same_solution("lpa() against multilpa(id = NULL)", base, direct, 1e-10),
            sweep_same_solution("lpa() against the enumerate_lpa() candidate", base,
                                get_candidate(grid, n_profiles = k, model = model), 1e-10))
    }),
    invariance("row-order", function(data, k, model) {
      shuffled <- sweep_seeded(5L, function() data[sample(nrow(data)), , drop = FALSE])
      sweep_same_solution("rows shuffled", fit_lpa(data, k, model = model),
                          fit_lpa(shuffled, k, model = model))
    }),
    invariance("column-order", function(data, k, model) {
      sweep_same_solution("indicators reversed", fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, vars = rev(names(data))))
    }),
    invariance("cores", function(data, k, model) {
      sweep_same_solution("cores = 2 against cores = 1", fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, cores = 2L), 0)
    }),
    invariance("acceleration", function(data, k, model) {
      sweep_same_solution("acceleration = none against squarem",
                          fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, acceleration = "none"), 1e-6)
    }),
    invariance("search", function(data, k, model) {
      base <- fit_lpa(data, k, model = model, n_starts = 20L)
      rbind(sweep_same_solution("search = kmeans", base,
                                fit_lpa(data, k, model = model, n_starts = 20L, search = "kmeans"),
                                1e-6),
            sweep_same_solution("search = emEM", base,
                                fit_lpa(data, k, model = model, n_starts = 20L, search = "emEM"),
                                1e-6),
            sweep_same_solution("search = emEM, n_screen = 200", base,
                                fit_lpa(data, k, model = model, n_starts = 20L, search = "emEM",
                                        n_screen = 200L), 1e-6))
    }),
    invariance("select_start", function(data, k, model) {
      sweep_same_solution("select_start = converged", fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, select_start = "converged"), 1e-8)
    }),
    invariance("boundary", function(data, k, model) {
      sweep_same_solution("boundary = exclude on an interior solution",
                          fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, boundary = "exclude"), 1e-8)
    }),
    invariance("restart", function(data, k, model) {
      base <- fit_lpa(data, k, model = model)
      again <- fit_lpa(data, k, model = model, start = starting_values(base), n_starts = 1L)
      rbind(sweep_same_solution("refit from starting_values()", base, again, 1e-8),
            sweep_true("the refit needs few iterations",
                       sweep_table(again, "model")$iterations <= 5L,
                       sprintf("%d iterations", sweep_table(again, "model")$iterations)))
    }),
    invariance("fiml-on-complete-data", function(data, k, model) {
      sweep_same_solution("missing = fiml on complete data", fit_lpa(data, k, model = model),
                          fit_lpa(data, k, model = model, missing = "fiml"), 1e-8)
    }),
    invariance("unit-weights", function(data, k, model) {
      weighted <- cbind(data, w = 1)
      sweep_same_solution("weights of one", fit_lpa(data, k, model = model),
                          fit_lpa(weighted, k, model = model, vars = names(data), weights = "w"),
                          1e-8)
    }),
    invariance("integer-weights", function(data, k, model) {
      w <- sweep_seeded(9L, function() sample(1:3, nrow(data), TRUE))
      repeated <- data[rep(seq_len(nrow(data)), w), , drop = FALSE]
      weighted <- fit_lpa(cbind(data, w = w), k, model = model, vars = names(data), weights = "w")
      expanded <- fit_lpa(repeated, k, model = model, start = starting_values(weighted),
                          n_starts = 1L)
      scale <- nrow(data) / sum(w)
      rbind(
        sweep_check("integer weights against repeated rows: pseudo log likelihood",
                    sweep_table(weighted, "information_criteria")$log_likelihood,
                    scale * sweep_table(expanded, "information_criteria")$log_likelihood, 1e-6),
        sweep_check("integer weights against repeated rows: sorted means",
                    sort(sweep_parameters(weighted)$means), sort(sweep_parameters(expanded)$means),
                    1e-3),
        sweep_check("integer weights against repeated rows: sorted shares",
                    sort(sweep_parameters(weighted)$mixing),
                    sort(sweep_parameters(expanded)$mixing), 1e-3))
    }),
    invariance("zero-weights", function(data, k, model) {
      w <- rep(c(1, 0), c(nrow(data) - 20L, 20L))
      kept <- data[w > 0, , drop = FALSE]
      sweep_same_solution("zero-weight rows against their removal",
                          fit_lpa(cbind(data, w = w), k, model = model, vars = names(data),
                                  weights = "w"),
                          fit_lpa(kept, k, model = model), 1e-6)
    }),
    invariance("location-scale", function(data, k, model) {
      factors <- seq(0.5, 3, length.out = ncol(data))
      moved <- as.data.frame(Map(function(column, factor) 10 + factor * column, data, factors))
      base <- fit_lpa(data, k, model = model)
      other <- fit_lpa(moved, k, model = model)
      # Only structures closed under a rescaling of each indicator keep the solution.
      if (!model %in% c("EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVV", "VVV")) {
        return(sweep_no_reference("rescaled indicators: log likelihood",
                                  "this structure is not closed under rescaling"))
      }
      sweep_check("rescaled indicators: log likelihood shifts by the Jacobian",
                  sweep_table(other, "information_criteria")$log_likelihood,
                  sweep_table(base, "information_criteria")$log_likelihood -
                    nrow(data) * sum(log(factors)), 1e-6)
    }),
    invariance("grand-scaling", function(data, k, model) {
      population_sd <- function(v) sqrt(mean((v - mean(v))^2))
      standardized <- as.data.frame(lapply(data, function(v) (v - mean(v)) / population_sd(v)))
      rbind(
        sweep_same_solution("centering and scaling = grand against standardized data",
                            fit_lpa(data, k, model = model, centering = "grand",
                                    scaling = "grand"),
                            fit_lpa(standardized, k, model = model), 1e-6),
        sweep_same_solution("centering = grand against centred data",
                            fit_lpa(data, k, model = model, centering = "grand"),
                            fit_lpa(as.data.frame(lapply(data, function(v) v - mean(v))), k,
                                    model = model), 1e-6))
    }),
    invariance("input-types", function(data, k, model) {
      base <- fit_lpa(data, k, model = model)
      integers <- as.data.frame(lapply(data, function(v) round(v * 10)))
      rbind(
        sweep_same_solution("extra unused columns", base,
                            fit_lpa(cbind(data, unused = "text", other = NA), k, model = model,
                                    vars = names(data)), 0),
        sweep_same_solution("row names set", base,
                            fit_lpa(`rownames<-`(data, paste0("r", seq_len(nrow(data)))), k,
                                    model = model), 0),
        sweep_same_solution("integer columns against the same numbers as doubles",
                            fit_lpa(as.data.frame(lapply(integers, as.integer)), k, model = model),
                            fit_lpa(integers, k, model = model), 0))
    }))

  # ---- missing data ------------------------------------------------------------
  missing_cases <- lapply(
    sweep_grid(dataset = names(datasets), model = c("EEI", "VVI", "EEE", "VVV"), k = 1:3),
    function(option) {
      sweep_define(sprintf("fiml[data=%s,model=%s,k=%d]", option$dataset, option$model,
                           option$k), function() {
        data <- sweep_punch(datasets[[option$dataset]], 0.12)
        fit <- fit_lpa(data, option$k, model = option$model, missing = "fiml")
        observed <- rowSums(!is.na(data)) > 0L
        rows <- sweep_likelihood_checks(fit, data, n_bic = sum(observed), empty = !observed)
        refit <- fit_lpa(data, option$k, model = option$model, missing = "fiml",
                         start = starting_values(fit), n_starts = 1L)
        rows <- rbind(rows, sweep_same_solution("refit from starting_values()", fit, refit, 1e-8),
                      sweep_same_solution("drop_empty = FALSE keeps the estimates", fit,
                                          fit_lpa(data, option$k, model = option$model,
                                                  missing = "fiml", drop_empty = FALSE), 1e-10))
        if (option$k > 1L || !option$model %in% c("EEI", "VVI")) return(rows)
        # One profile, diagonal: each indicator's own observed mean and ML variance.
        parameters <- sweep_parameters(fit)
        rbind(rows,
              sweep_check("one profile: means are the observed means", parameters$means,
                          vapply(data, mean, numeric(1), na.rm = TRUE), 1e-8),
              sweep_check("one profile: variances are the observed ML variances",
                          diag(parameters$covariances[[1L]]),
                          vapply(data, function(v) mean((v - mean(v, na.rm = TRUE))^2, na.rm = TRUE),
                                 numeric(1)), 1e-8))
      })
    })

  # ---- one profile: closed forms under weights -----------------------------------
  weighted_cases <- lapply(names(datasets), function(dataset) {
    sweep_define(sprintf("weighted-one-profile[data=%s]", dataset), function() {
      data <- datasets[[dataset]]
      w <- sweep_seeded(13L, function() stats::runif(nrow(data), 0.2, 3))
      fit <- fit_lpa(cbind(data, w = w), 1L, model = "VVV", vars = names(data), weights = "w")
      moments <- stats::cov.wt(data, wt = w, method = "ML")
      parameters <- sweep_parameters(fit)
      scaled <- w * nrow(data) / sum(w)
      rbind(sweep_check("weighted means", parameters$means, moments$center, 1e-8),
            sweep_check("weighted ML covariance", parameters$covariances[[1L]], moments$cov, 1e-8),
            sweep_likelihood_checks(fit, data, weights = scaled))
    })
  })

  # ---- prior and noise: relations to the maximum likelihood fit ----------------------
  prior_noise_cases <- lapply(settings, function(s) {
    sweep_define(sprintf("prior-and-noise[data=%s,model=%s,k=%d]", s$dataset, s$model, s$k),
                 function() {
      data <- datasets[[s$dataset]]
      ll_of <- function(x) sweep_table(x, "information_criteria")$log_likelihood
      mclust_mode <- function(label, model, prior) {
        continued <- sweep_reference(mclust::me(
          as.matrix(data), modelName = s$model, z = sweep_posteriors(model), prior = prior,
          control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(20000L, 20000L))))
        quantity <- sprintf("%s: mclust's EM under the same prior, from the package's posteriors",
                            label)
        if (!is.finite(continued$loglik)) {
          return(sweep_no_reference(quantity, "mclust has no fit of this structure under a prior"))
        }
        sweep_check(quantity, ll_of(model), continued$loglik, 1e-6)
      }
      base <- fit_lpa(data, s$k, model = s$model)
      regularized <- fit_lpa(data, s$k, model = s$model, prior = prior_control())
      weaker <- fit_lpa(data, s$k, model = s$model, prior = prior_control(shrinkage = 1e-8))
      noisy <- fit_lpa(data, s$k, model = s$model, noise = TRUE)
      rbind(
        # The likelihood at a posterior mode cannot exceed the maximum likelihood.
        sweep_at_least("maximum likelihood at least the likelihood at the posterior mode",
                       ll_of(base), ll_of(regularized), 1e-8),
        sweep_check("two runs with one seed under a prior", ll_of(regularized),
                    ll_of(fit_lpa(data, s$k, model = s$model, prior = prior_control())), 0),
        # A noise component of share zero is the model without it.
        sweep_at_least("a noise component cannot lower the maximum", ll_of(noisy), ll_of(base),
                       1e-6),
        sweep_check("two runs with one seed with noise", ll_of(noisy),
                    ll_of(fit_lpa(data, s$k, model = s$model, noise = TRUE)), 0),
        # The gap to the maximum likelihood is not monotone in the prior's
        # strength (diagnostics/prior-strength.R), so each posterior mode is
        # checked against mclust's EM under the same prior instead.
        mclust_mode("default prior", regularized, mclust::priorControl()),
        mclust_mode("shrinkage 1e-8", weaker, mclust::priorControl(shrinkage = 1e-8)))
    })
  })

  # ---- arguments that must be refused with a classed condition -----------------------
  refusal_cases <- list(sweep_define("refusals", function() {
    data <- datasets$blobs
    vars <- names(data)
    holes <- sweep_punch(data, 0.1)
    rbind(
      sweep_refuses("n_profiles = 0", function() lpa(data, vars, n_profiles = 0L)),
      sweep_refuses("n_profiles = 1.5", function() lpa(data, vars, n_profiles = 1.5)),
      sweep_refuses("n_profiles = NA", function() lpa(data, vars, n_profiles = NA)),
      sweep_refuses("n_profiles of length two", function() lpa(data, vars, n_profiles = 1:2)),
      sweep_refuses("n_profiles above the number of rows",
                    function() lpa(data[1:3, ], vars, n_profiles = 5L)),
      sweep_refuses("an unknown indicator", function() lpa(data, c(vars, "nope"), 2L)),
      sweep_refuses("a repeated indicator", function() lpa(data, c(vars, "a"), 2L)),
      sweep_refuses("no indicator", function() lpa(data, character(), 2L)),
      sweep_refuses("a text indicator",
                    function() lpa(cbind(data, t = "x"), c(vars, "t"), 2L)),
      sweep_refuses("a constant indicator",
                    function() lpa(cbind(data, t = 1), c(vars, "t"), 2L)),
      sweep_refuses("an infinite value",
                    function() lpa(replace(data, 1L, list(a = replace(data$a, 2L, Inf))), vars, 2L)),
      sweep_refuses("missing values without fiml", function() lpa(holes, vars, 2L)),
      sweep_refuses("a matrix as data", function() lpa(as.matrix(data), vars, 2L)),
      sweep_refuses("an empty data frame", function() lpa(data[0L, ], vars, 2L)),
      sweep_refuses("an unknown model code", function() lpa(data, vars, 2L, model = "XYZ")),
      sweep_refuses("model with volume", function() lpa(data, vars, 2L, model = "VVI",
                                                        volume = "equal")),
      sweep_refuses("n_starts = 0", function() lpa(data, vars, 2L, n_starts = 0L)),
      sweep_refuses("max_iter = -1", function() lpa(data, vars, 2L, max_iter = -1L)),
      sweep_refuses("tol = -1", function() lpa(data, vars, 2L, tol = -1)),
      sweep_refuses("min_variance = -1", function() lpa(data, vars, 2L, min_variance = -1)),
      sweep_refuses("seed = 'a'", function() lpa(data, vars, 2L, seed = "a")),
      sweep_refuses("cores = 0", function() lpa(data, vars, 2L, cores = 0L)),
      sweep_refuses("an unknown missing option", function() lpa(data, vars, 2L, missing = "drop")),
      sweep_refuses("an unknown weights column", function() lpa(data, vars, 2L, weights = "w")),
      sweep_refuses("a prior that is not a prior_control()",
                    function() lpa(data, vars, 2L, prior = list(shrinkage = 1))),
      sweep_refuses("prior_control(shrinkage = -1)", function() prior_control(shrinkage = -1)),
      sweep_refuses("noise = NA", function() lpa(data, vars, 2L, noise = NA)),
      sweep_refuses("a prior with weights",
                    function() lpa(cbind(data, w = 1), vars, 2L, weights = "w",
                                   prior = prior_control())),
      sweep_refuses("negative weights",
                    function() lpa(cbind(data, w = -1), vars, 2L, weights = "w")),
      sweep_refuses("missing weights",
                    function() lpa(cbind(data, w = NA_real_), vars, 2L, weights = "w")),
      sweep_refuses("id passed to lpa()", function() lpa(data, vars, 2L, id = "a")),
      sweep_refuses("n_screen without emEM",
                    function() lpa(data, vars, 2L, search = "kmeans", n_screen = 50L)),
      sweep_refuses("a start of the wrong shape",
                    function() lpa(data, vars, 2L, start = starting_values(lpa(data, vars, 3L, seed = 1L)))))
  }))

  c(structure_cases, spelling_cases, invariance_cases, missing_cases, weighted_cases,
    prior_noise_cases, refusal_cases)
}
