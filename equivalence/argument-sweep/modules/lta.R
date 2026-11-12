# lta(): the transition arguments (`transitions`, `order`, `mover_stayer`,
# `n_group_classes`, `measurement`, `occasions`), the measurement arguments it
# shares with multilpa(), missing data, weights and the types of `id` and
# `time`, on six longitudinal datasets. References: the hidden Markov
# likelihood computed here by the forward recursion from the fit's result
# tables (first and second order, several chain classes, movers and stayers,
# occasion-specific transitions and measurement), smoothed posteriors by the
# backward recursion, nesting relations between the options, and invariances.

sweep_module_lta <- function() {
  sets <- sweep_longitudinal()
  fit <- function(set, k, ..., data = set$data, id = set$id, time = set$time, n_starts = 8L) {
    roles <- set[intersect(names(set), c("categorical", "count", "ordinal"))]
    do.call(lta, c(list(data, set$vars, id = id, n_profiles = k, time = time, seed = 71L,
                        n_starts = n_starts), roles, list(...)))
  }
  index_of <- function(label) {
    if (is.numeric(label)) as.integer(label) else as.integer(sub("^.*_", "", as.character(label)))
  }
  ll_of <- function(x) x$log_likelihood
  par_of <- function(x) x$n_parameters

  #' Emission parameters in the shape `sweep_log_density()` reads
  emissions <- function(tables, occasion = NULL) {
    of_occasion <- function(table) {
      if (is.null(occasion) || !"occasion" %in% names(table) || all(is.na(table$occasion))) {
        return(table)
      }
      table[table$occasion == occasion, , drop = FALSE]
    }
    profiles <- of_occasion(tables$profiles)
    profiles$profile <- index_of(profiles$profile)
    continuous <- unique(as.character(profiles$indicator))
    K <- max(index_of(tables$initial$profile))
    levels <- seq_len(K)
    covariances <- if (!length(continuous)) NULL else if (is.data.frame(tables$covariances)) {
      table <- of_occasion(tables$covariances)
      table$profile <- index_of(table$profile)
      lapply(levels, function(k) {
        sweep_wide(table[table$profile == k, , drop = FALSE], "indicator", "indicator_2",
                   "covariance", continuous, continuous)
      })
    } else {
      variances <- sweep_wide(profiles, "profile", "indicator", "variance", levels, continuous)
      lapply(levels, function(k) diag(variances[k, ], length(continuous)))
    }
    shared <- c("profile", "indicator", "category", "probability")
    discrete <- do.call(rbind, lapply(c("responses", "ordinal"), function(name) {
      table <- tables[[name]]
      if (!is.data.frame(table) || !nrow(table)) return(NULL)
      table <- of_occasion(table)[shared]
      table$profile <- index_of(table$profile)
      table
    }))
    counts <- if (is.data.frame(tables$count_means) && nrow(tables$count_means)) {
      table <- of_occasion(tables$count_means)
      table$profile <- index_of(table$profile)
      table
    }
    list(profile_levels = levels, continuous = continuous,
         means = if (length(continuous)) sweep_wide(profiles, "profile", "indicator", "mean",
                                                    levels, continuous),
         covariances = covariances, responses = discrete, counts = counts)
  }

  #' The chain classes of a fit: weight, initial vector and transition matrices
  chains <- function(tables, K, occasions) {
    initial <- tables$initial
    classes <- unique(as.character(initial$group_class))
    weights <- if (is.data.frame(tables$group_classes)) {
      tables$group_classes$weight[match(classes, as.character(tables$group_classes$group_class))]
    } else initial$group_class_probability[match(classes, as.character(initial$group_class))]
    transitions <- tables$transitions
    second <- tables$second_order_transitions
    lapply(seq_along(classes), function(c) {
      start <- initial[as.character(initial$group_class) == classes[[c]], , drop = FALSE]
      moves <- transitions[as.character(transitions$group_class) == classes[[c]], , drop = FALSE]
      matrix_of <- function(rows) {
        out <- matrix(0, K, K)
        out[cbind(index_of(rows$from), index_of(rows$to))] <- rows$probability
        out
      }
      by_occasion <- "occasion" %in% names(moves) && !all(is.na(moves$occasion))
      first_order <- lapply(occasions[-1L], function(occasion) {
        matrix_of(if (by_occasion) moves[moves$occasion == occasion, , drop = FALSE] else moves)
      })
      second_order <- if (is.data.frame(second)) {
        rows <- second[as.character(second$group_class) == classes[[c]], , drop = FALSE]
        out <- array(0, c(K, K, K))
        out[cbind(index_of(rows$previous), index_of(rows$from), index_of(rows$to))] <-
          rows$probability
        out
      }
      list(weight = weights[[c]],
           initial = start$probability[order(index_of(start$profile))],
           first_order = first_order, second_order = second_order)
    })
  }

  #' Log likelihood of one sequence in one chain class, with smoothed posteriors
  #' @param log_b Occasions x K emission log densities of the unit's rows.
  #' @param position Index of each row's occasion among the fitted occasions.
  forward <- function(log_b, position, chain) {
    b <- exp(log_b - apply(log_b, 1L, max))
    offset <- sum(apply(log_b, 1L, max))
    K <- ncol(b)
    steps <- nrow(b)
    if (!is.null(chain$second_order) && steps >= 2L) {
      # State pairs (previous, current); the first step uses the first-order matrix.
      alpha <- (chain$initial * b[1L, ]) * chain$first_order[[1L]] *
        matrix(b[2L, ], K, K, byrow = TRUE)
      scale <- sum(alpha)
      log_scale <- log(scale)
      alpha <- alpha / scale
      final <- Reduce(function(state, t) {
        moved <- vapply(seq_len(K), function(k) colSums(state$alpha * chain$second_order[, , k]),
                        numeric(K))
        moved <- moved * matrix(b[t, ], K, K, byrow = TRUE)
        list(alpha = moved / sum(moved), log_scale = state$log_scale + log(sum(moved)))
      }, seq_len(steps)[-(1:2)], list(alpha = alpha, log_scale = log_scale))
      return(list(log_likelihood = final$log_scale + offset, posteriors = NULL))
    }
    alphas <- matrix(0, steps, K)
    scales <- numeric(steps)
    alphas[1L, ] <- chain$initial * b[1L, ]
    scales[[1L]] <- sum(alphas[1L, ])
    alphas[1L, ] <- alphas[1L, ] / scales[[1L]]
    # A genuine recursion: each occasion needs the previous one.
    for (t in seq_len(steps)[-1L]) {
      moved <- as.numeric(alphas[t - 1L, ] %*% chain$first_order[[position[[t]] - 1L]]) * b[t, ]
      scales[[t]] <- sum(moved)
      alphas[t, ] <- moved / scales[[t]]
    }
    betas <- matrix(1, steps, K)
    for (t in rev(seq_len(steps))[-1L]) {
      betas[t, ] <- as.numeric(chain$first_order[[position[[t + 1L]] - 1L]] %*%
                                 (b[t + 1L, ] * betas[t + 1L, ])) / scales[[t + 1L]]
    }
    list(log_likelihood = sum(log(scales)) + offset, posteriors = alphas * betas)
  }

  #' The likelihood and posteriors of a transition fit rebuilt from its tables
  #' @param model The fit.
  #' @param set The dataset description.
  #' @param data The fitted data.
  #' @param weights Weight per unit, named by id, or `NULL`.
  #' @param gaps `TRUE` when a missing occasion is a step of the chain (`occasions = "grid"`).
  rebuild <- function(model, set, data, id = set$id, time = set$time, weights = NULL) {
    tables <- suppressMessages(get_results(model, "all"))
    data <- data[order(match(data[[id]], unique(data[[id]])), data[[time]]), , drop = FALSE]
    occasions <- sort(unique(data[[time]]))
    K <- max(index_of(tables$initial$profile))
    varies <- function(table) {
      is.data.frame(table) && "occasion" %in% names(table) && !all(is.na(table$occasion))
    }
    by_occasion <- any(vapply(tables[c("profiles", "responses", "ordinal", "count_means")],
                              varies, logical(1)))
    log_b <- if (by_occasion) {
      out <- matrix(0, nrow(data), K)
      filled <- lapply(occasions, function(occasion) {
        rows <- which(data[[time]] == occasion)
        list(rows = rows, values = sweep_log_density(
          data[rows, , drop = FALSE], emissions(tables, match(occasion, occasions))))
      })
      out[unlist(lapply(filled, `[[`, "rows")), ] <- do.call(rbind, lapply(filled, `[[`, "values"))
      out
    } else sweep_log_density(data, emissions(tables))
    classes <- chains(tables, K, seq_along(occasions))
    homogeneous <- !("occasion" %in% names(tables$transitions)) ||
      all(is.na(tables$transitions$occasion))
    units <- split(seq_len(nrow(data)), match(data[[id]], unique(data[[id]])))
    per_unit <- lapply(units, function(rows) {
      # Under homogeneous transitions each observed step uses the one matrix.
      position <- if (homogeneous) seq_along(rows) else match(data[[time]][rows], occasions)
      runs <- lapply(classes, function(chain) {
        forward(log_b[rows, , drop = FALSE], position,
                if (homogeneous) c(chain[c("weight", "initial", "second_order")],
                                   list(first_order = rep(chain$first_order[1L], length(rows))))
                else chain)
      })
      joint <- log(vapply(classes, `[[`, numeric(1), "weight")) +
        vapply(runs, `[[`, numeric(1), "log_likelihood")
      total <- max(joint) + log(sum(exp(joint - max(joint))))
      smoothed <- if (all(!vapply(runs, function(r) is.null(r$posteriors), logical(1)))) {
        Reduce(`+`, Map(function(run, share) share * run$posteriors, runs, exp(joint - total)))
      }
      list(log_likelihood = total, posteriors = smoothed)
    })
    unit_weights <- if (is.null(weights)) rep(1, length(units)) else weights
    posteriors <- if (all(!vapply(per_unit, function(u) is.null(u$posteriors), logical(1)))) {
      out <- matrix(NA_real_, nrow(data), K)
      out[unlist(units), ] <- do.call(rbind, lapply(per_unit, `[[`, "posteriors"))
      out[order(as.integer(rownames(data))), , drop = FALSE]
    }
    list(log_likelihood = sum(unit_weights * vapply(per_unit, `[[`, numeric(1), "log_likelihood")),
         posteriors = posteriors,
         reported = as.matrix(tables$assignments[grep("^posterior_profile_",
                                                      names(tables$assignments))]))
  }

  #' Checks of a transition fit against the rebuilt likelihood
  rebuilt_checks <- function(model, set, data = set$data, ...) {
    rownames(data) <- NULL
    rebuilt <- rebuild(model, set, data, ...)
    rows <- rbind(
      sweep_check("log likelihood by the forward recursion from the result tables",
                  ll_of(model), rebuilt$log_likelihood, 1e-8),
      sweep_check("AIC from the likelihood and parameter count", stats::AIC(model),
                  -2 * ll_of(model) + 2 * par_of(model), 1e-10),
      sweep_check("BIC from the likelihood and parameter count", stats::BIC(model),
                  -2 * ll_of(model) + par_of(model) * log(stats::nobs(model)), 1e-10))
    if (is.null(rebuilt$posteriors)) return(rows)
    kept <- stats::complete.cases(rebuilt$reported)
    rbind(rows, sweep_check("smoothed posteriors by the backward recursion",
                            rebuilt$reported[kept, , drop = FALSE],
                            rebuilt$posteriors[kept, , drop = FALSE], 1e-6))
  }

  option_cases <- lapply(
    c(sweep_grid(set = c("gaussian", "three_state", "binary", "mixed", "hours", "ragged"),
                 option = c("default", "equal-variance", "group-classes", "mover-stayer")),
      sweep_grid(set = c("gaussian", "three_state", "binary", "mixed", "hours"),
                 option = c("occasion-transitions", "second-order", "occasion-measurement",
                            "occasion-transitions+group-classes", "full-covariance"))),
    function(option) {
      sweep_define(sprintf("options[data=%s,%s]", option$set, option$option), function() {
        set <- sets[[option$set]]
        k <- if (identical(option$set, "three_state")) 3L else 2L
        gaussian <- length(setdiff(set$vars, unlist(set[c("categorical", "count", "ordinal")])))
        if (!gaussian && option$option %in% c("equal-variance", "full-covariance")) {
          return(sweep_no_reference("log likelihood", "no continuous indicator"))
        }
        if (gaussian < 2L && identical(option$option, "full-covariance")) {
          return(sweep_no_reference("log likelihood", "one continuous indicator"))
        }
        arguments <- switch(option$option,
          "default" = list(),
          "equal-variance" = list(variance_model = "equal"),
          "full-covariance" = list(covariance_model = "full"),
          "group-classes" = list(n_group_classes = 2L),
          "mover-stayer" = list(mover_stayer = TRUE),
          "occasion-transitions" = list(transitions = "occasion"),
          "second-order" = list(order = 2L),
          "occasion-measurement" = list(measurement = "occasion"),
          "occasion-transitions+group-classes" = list(transitions = "occasion",
                                                      n_group_classes = 2L))
        model <- do.call(fit, c(list(set, k), arguments))
        base <- if (identical(option$option, "default")) model else fit(set, k)
        rows <- rebuilt_checks(model, set)
        occasions <- length(unique(set$data[[set$time]]))
        if (identical(option$option, "default") && gaussian == length(set$vars)) {
          rows <- rbind(rows, sweep_check("parameter count", par_of(model),
                                          (k - 1L) + k * (k - 1L) + 2L * k * gaussian, 0))
        }
        if (identical(option$option, "occasion-transitions") && gaussian == length(set$vars)) {
          rows <- rbind(rows, sweep_check(
            "parameter count", par_of(model),
            (k - 1L) + (occasions - 1L) * k * (k - 1L) + 2L * k * gaussian, 0))
        }
        # Each of these options contains the default model, so its maximum is not lower.
        if (option$option %in% c("group-classes", "mover-stayer", "occasion-transitions",
                                 "second-order", "occasion-measurement", "full-covariance",
                                 "occasion-transitions+group-classes")) {
          rows <- rbind(rows, sweep_at_least("likelihood at least that of the default model",
                                             ll_of(model), ll_of(base), 1e-6))
        }
        if (identical(option$option, "equal-variance")) {
          rows <- rbind(rows, sweep_at_least("likelihood at most that of varying variances",
                                             ll_of(base), ll_of(model), 1e-6))
        }
        rows
      })
    })

  settings <- c("gaussian", "binary", "ragged")
  invariance <- function(name, code) {
    lapply(settings, function(s) {
      sweep_define(sprintf("%s[data=%s]", name, s), function() code(sets[[s]], 2L))
    })
  }
  same <- function(label, first, second, tolerance = 1e-7) {
    rbind(sweep_check(sprintf("%s: log likelihood", label), ll_of(first), ll_of(second), tolerance),
          sweep_check(sprintf("%s: parameter count", label), par_of(first), par_of(second), 0),
          sweep_check(sprintf("%s: sorted transition probabilities", label),
                      sort(sweep_table(first, "transitions")$probability),
                      sort(sweep_table(second, "transitions")$probability), sqrt(tolerance)))
  }
  invariance_cases <- c(
    invariance("seed-reproducible", function(set, k) {
      same("two runs with one seed", fit(set, k), fit(set, k), 0)
    }),
    invariance("row-order", function(set, k) {
      shuffled <- sweep_seeded(5L, function() set$data[sample(nrow(set$data)), , drop = FALSE])
      same("rows shuffled", fit(set, k), fit(set, k, data = shuffled))
    }),
    invariance("identifier-and-time-types", function(set, k) {
      base <- fit(set, k)
      data <- set$data
      data$key <- paste0("unit", data[[set$id]])
      data$when <- factor(data[[set$time]], levels = sort(unique(data[[set$time]])))
      data$day <- as.Date("2026-01-01") + 7L * data[[set$time]]
      data$stamp <- 10 * data[[set$time]] + 0.5
      rbind(same("character identifier", base, fit(set, k, data = data, id = "key"), 1e-10),
            same("factor time", base, fit(set, k, data = data, time = "when"), 1e-10),
            same("Date time", base, fit(set, k, data = data, time = "day"), 1e-10),
            same("unequally labelled numeric time", base, fit(set, k, data = data, time = "stamp"),
                 1e-10))
    }),
    invariance("cores", function(set, k) {
      same("cores = 2 against cores = 1", fit(set, k), fit(set, k, cores = 2L), 0)
    }),
    invariance("acceleration", function(set, k) {
      same("acceleration = none against squarem", fit(set, k), fit(set, k, acceleration = "none"),
           1e-6)
    }),
    invariance("search", function(set, k) {
      same("search = emEM against kmeans", fit(set, k, n_starts = 16L),
           fit(set, k, n_starts = 16L, search = "emEM"), 1e-6)
    }),
    invariance("select_start", function(set, k) {
      same("select_start = converged", fit(set, k), fit(set, k, select_start = "converged"), 1e-8)
    }),
    invariance("fiml-on-complete-data", function(set, k) {
      same("missing = fiml on complete data", fit(set, k), fit(set, k, missing = "fiml"), 1e-8)
    }),
    invariance("unit-weights", function(set, k) {
      same("weights of one", fit(set, k), fit(set, k, data = cbind(set$data, w = 1), weights = "w"),
           1e-8)
    }),
    invariance("integer-weights", function(set, k) {
      ids <- unique(set$data[[set$id]])
      w <- sweep_seeded(9L, function() sample(1:3, length(ids), TRUE))
      data <- cbind(set$data, w = w[match(set$data[[set$id]], ids)])
      copies <- do.call(rbind, lapply(1:3, function(copy) {
        rows <- data[data$w >= copy, , drop = FALSE]
        rows$key <- paste(rows[[set$id]], copy)
        rows
      }))
      weighted <- fit(set, k, data = data, weights = "w")
      expanded <- fit(set, k, data = copies, id = "key")
      rbind(
        sweep_check("integer weights against repeated units: pseudo log likelihood",
                    ll_of(weighted), length(ids) / sum(w) * ll_of(expanded), 1e-6),
        sweep_check("integer weights against repeated units: sorted transition probabilities",
                    sort(sweep_table(weighted, "transitions")$probability),
                    sort(sweep_table(expanded, "transitions")$probability), 1e-3),
        rebuilt_checks(weighted, set, data = data, weights = w * length(ids) / sum(w)))
    }),
    invariance("two-occasions", function(set, k) {
      # With two occasions there is one transition, so the options coincide.
      data <- set$data[set$data[[set$time]] <= 2L, , drop = FALSE]
      same("occasion-specific against homogeneous transitions", fit(set, k, data = data),
           fit(set, k, data = data, transitions = "occasion"), 1e-6)
    }))

  missing_cases <- lapply(c("gaussian", "three_state", "binary", "mixed", "hours"), function(s) {
    sweep_define(sprintf("fiml[data=%s]", s), function() {
      set <- sets[[s]]
      data <- sweep_punch(set$data, 0.12, columns = set$vars)
      k <- if (identical(s, "three_state")) 3L else 2L
      model <- fit(set, k, data = data, missing = "fiml")
      rbind(rebuilt_checks(model, set, data = data),
            sweep_check("drop_empty = FALSE keeps the likelihood", ll_of(model),
                        ll_of(fit(set, k, data = data, missing = "fiml", drop_empty = FALSE)),
                        1e-10))
    })
  })

  grid_cases <- list(
    sweep_define("grid-on-complete-occasions", function() {
      set <- sets$gaussian
      same("occasions = grid against observed on a complete panel", fit(set, 2L),
           fit(set, 2L, occasions = "grid"), 1e-8)
    }),
    sweep_define("grid-is-fiml-on-the-filled-panel", function() {
      # A skipped occasion on the grid is a row with nothing observed.
      set <- sets$ragged
      panel <- merge(expand.grid(id = unique(set$data$id), wave = 1:4), set$data, all.x = TRUE)
      on_grid <- fit(set, 2L, occasions = "grid")
      filled <- fit(set, 2L, data = panel, missing = "fiml")
      rbind(sweep_check("log likelihood", ll_of(on_grid), ll_of(filled), 1e-7),
            sweep_check("sorted transition probabilities",
                        sort(sweep_table(on_grid, "transitions")$probability),
                        sort(sweep_table(filled, "transitions")$probability), 1e-4))
    }))

  refusal_cases <- list(sweep_define("refusals", function() {
    set <- sets$gaussian
    data <- set$data
    call <- function(...) function() lta(data, set$vars, id = "id", time = "wave", ...)
    rbind(
      sweep_refuses("time = NULL", function() lta(data, set$vars, id = "id", n_profiles = 2L,
                                                  time = NULL)),
      sweep_refuses("an unknown time column",
                    function() lta(data, set$vars, id = "id", n_profiles = 2L, time = "nope")),
      sweep_refuses("time repeated within a unit",
                    function() lta(transform(data, wave = 1L), set$vars, id = "id",
                                   n_profiles = 2L, time = "wave")),
      sweep_refuses("missing time",
                    function() lta(transform(data, wave = replace(wave, 2L, NA)), set$vars,
                                   id = "id", n_profiles = 2L, time = "wave")),
      sweep_refuses("one profile", call(n_profiles = 1L)),
      sweep_refuses("order = 3", call(n_profiles = 2L, order = 3L)),
      sweep_refuses("order = 0", call(n_profiles = 2L, order = 0L)),
      sweep_refuses("mover_stayer = NA", call(n_profiles = 2L, mover_stayer = NA)),
      sweep_refuses("an unknown transitions option", call(n_profiles = 2L, transitions = "free")),
      sweep_refuses("an unknown measurement option", call(n_profiles = 2L, measurement = "free")),
      sweep_refuses("an unknown occasions option", call(n_profiles = 2L, occasions = "all")),
      sweep_refuses("an unknown transition covariate",
                    call(n_profiles = 2L, transition_covariates = "nope")),
      sweep_refuses("the time column among the indicators",
                    function() lta(data, c(set$vars, "wave"), id = "id", n_profiles = 2L,
                                   time = "wave")),
      sweep_refuses("every unit observed once",
                    function() lta(data[data$wave == 1L, ], set$vars, id = "id", n_profiles = 2L,
                                   time = "wave")))
  }))

  c(option_cases, invariance_cases, missing_cases, grid_cases, refusal_cases)
}
