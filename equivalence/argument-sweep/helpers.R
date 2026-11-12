# Shared plumbing of the argument sweep. A module returns a list of cases; a
# case is a name and a function that fits models and returns checks. Every
# check compares a package value with a reference value computed without the
# package's own likelihood code: a closed form, another package, or the
# mixture likelihood rebuilt here from the fit's public result tables.

#' One comparison of package values with reference values
#' @param quantity What is compared.
#' @param value,reference Numeric vectors of equal length.
#' @param tolerance Largest accepted difference, relative to `max(1, |reference|)`.
#' @param note Optional context.
#' @return One-row data.frame: quantity, difference, tolerance, n_values, status, note.
sweep_check <- function(quantity, value, reference, tolerance, note = NA_character_) {
  stopifnot("`quantity` must be one string" = is.character(quantity) && length(quantity) == 1L,
            "`tolerance` must be one nonnegative number" =
              is.numeric(tolerance) && length(tolerance) == 1L && tolerance >= 0)
  value <- as.numeric(value)
  reference <- as.numeric(reference)
  comparable <- length(value) == length(reference) && length(value) > 0L &&
    !anyNA(value) && !anyNA(reference)
  difference <- if (comparable) max(abs(value - reference) / pmax(1, abs(reference))) else Inf
  if (!comparable) {
    note <- sprintf("not comparable: %d package values (%d missing), %d reference values (%d missing)%s",
                    length(value), sum(is.na(value)), length(reference), sum(is.na(reference)),
                    if (is.na(note)) "" else paste0("; ", note))
  }
  data.frame(quantity = quantity, difference = difference, tolerance = tolerance,
             n_values = length(reference),
             status = if (difference <= tolerance) "agrees" else "disagrees",
             note = note, stringsAsFactors = FALSE)
}

#' A check of a logical condition
#' @param quantity What must hold.
#' @param condition One logical.
#' @param note Optional context.
#' @return One-row data.frame in the shape of `sweep_check()`.
sweep_true <- function(quantity, condition, note = NA_character_) {
  holds <- isTRUE(condition)
  data.frame(quantity = quantity, difference = if (holds) 0 else 1, tolerance = 0,
             n_values = 1L, status = if (holds) "agrees" else "disagrees",
             note = note, stringsAsFactors = FALSE)
}

#' A check that could not be made, recorded rather than dropped
#' @param quantity What would have been compared.
#' @param note Why there is no reference.
#' @return One-row data.frame with status "no reference".
sweep_no_reference <- function(quantity, note) {
  data.frame(quantity = quantity, difference = NA_real_, tolerance = NA_real_,
             n_values = 0L, status = "no reference", note = note, stringsAsFactors = FALSE)
}

#' The upper bound of a one-sided comparison: the package may not be worse
#' @param quantity What is compared.
#' @param value Package value (larger is better).
#' @param reference Reference value.
#' @param tolerance Accepted shortfall, relative to `max(1, |reference|)`.
#' @return One-row data.frame in the shape of `sweep_check()`.
sweep_at_least <- function(quantity, value, reference, tolerance, note = NA_character_) {
  shortfall <- max(0, (reference - value) / max(1, abs(reference)))
  if (!is.finite(shortfall)) shortfall <- Inf
  data.frame(quantity = quantity, difference = shortfall, tolerance = tolerance, n_values = 1L,
             status = if (shortfall <= tolerance) "agrees" else "disagrees",
             note = note, stringsAsFactors = FALSE)
}

#' Expect a classed error
#' @param quantity What is refused.
#' @param expression A function with no arguments that should fail.
#' @param class Accepted condition class prefix.
#' @return One-row data.frame: agrees when a `latents_` classed error is raised.
sweep_refuses <- function(quantity, expression, class = "^(latents|multilpa)_") {
  outcome <- tryCatch({
    expression()
    "no error"
  }, error = function(e) {
    classes <- class(e)
    if (any(grepl(class, classes))) "classed" else
      sprintf("unclassed error: %s", conditionMessage(e))
  })
  sweep_true(quantity, identical(outcome, "classed"),
             note = if (identical(outcome, "classed")) NA_character_ else outcome)
}

#' Expect a classed warning and a result
#' @param quantity What is accepted with a warning.
#' @param expression A function with no arguments.
#' @param class The warning class.
#' @return One-row data.frame: agrees when the call returns and signals `class`.
sweep_warns <- function(quantity, expression, class) {
  signalled <- FALSE
  outcome <- tryCatch(withCallingHandlers({
    expression()
    "returned"
  }, warning = function(w) {
    if (inherits(w, class)) signalled <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(e) sprintf("error: %s", conditionMessage(e)))
  sweep_true(quantity, identical(outcome, "returned") && signalled,
             note = if (!identical(outcome, "returned")) outcome else
               if (!signalled) sprintf("no `%s` warning", class) else NA_character_)
}

#' Define a case
#' @param case Name.
#' @param code Function with no arguments returning checks.
#' @return A list with `case` and `code`.
sweep_define <- function(case, code) {
  stopifnot(is.character(case), length(case) == 1L, is.function(code))
  list(case = case, code = code)
}

#' Checks as one data.frame, whether a case returned one or a list of them
sweep_bind <- function(checks) {
  if (is.data.frame(checks)) checks else do.call(rbind, checks)
}

#' Run one case, recording an error or warnings as data
#' @param module Module name.
#' @param definition A `sweep_define()` value.
#' @return data.frame of checks with module, case, conditions and seconds.
sweep_run_case <- function(module, definition) {
  seen <- character()
  started <- proc.time()[["elapsed"]]
  rows <- withCallingHandlers(
    tryCatch(sweep_bind(definition$code()), error = function(e) {
      data.frame(quantity = "(case ran)", difference = Inf, tolerance = 0, n_values = 0L,
                 status = "error",
                 note = sprintf("%s [%s]", conditionMessage(e), paste(class(e), collapse = "/")),
                 stringsAsFactors = FALSE)
    }),
    warning = function(w) {
      seen <<- c(seen, sprintf("%s: %s", class(w)[[1L]], conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  cbind(module = module, case = definition$case, rows,
        conditions = if (length(seen)) paste(unique(substr(seen, 1L, 160L)), collapse = " | ")
                     else NA_character_,
        seconds = round(proc.time()[["elapsed"]] - started, 2), stringsAsFactors = FALSE)
}

#' Expand named option vectors into a list of named option lists
#' @param ... Named vectors or lists of options.
#' @return A list with one named list per combination.
sweep_grid <- function(...) {
  grid <- expand.grid(..., stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, , drop = FALSE]))
}

# ---- the fit's public tables as arrays ------------------------------------

#' A long table as a matrix
#' @param table data.frame.
#' @param rows,columns,value Column names.
#' @param row_levels,column_levels Orders of the dimensions.
#' @return Numeric matrix.
sweep_wide <- function(table, rows, columns, value, row_levels = unique(table[[rows]]),
                       column_levels = unique(table[[columns]])) {
  out <- matrix(NA_real_, length(row_levels), length(column_levels),
                dimnames = list(as.character(row_levels), as.character(column_levels)))
  out[cbind(match(table[[rows]], row_levels), match(table[[columns]], column_levels))] <-
    table[[value]]
  out
}

#' One table of a fit, without the accessor's progress messages
sweep_table <- function(fit, what) suppressMessages(get_results(fit, what))

#' Parameters of a multilpa-family fit, read from its result tables
#' @param fit A fit of `multilpa()`, `lpa()`, `lca()` or `multilca()`.
#' @return List: `means` (K x P), `covariances` (list of K), `mixing` (G x K),
#'   `group_mixing` (G), `responses`, `counts`, and the indicator names.
sweep_parameters <- function(fit) {
  profiles <- sweep_table(fit, "profiles")
  # A fit with membership covariates has logits in place of a table of shares.
  mixing_table <- if (inherits(fit, "multilpa_covariates")) NULL else
    sweep_table(fit, "profile_probabilities")
  profile_levels <- if (is.null(mixing_table)) seq_len(fit$n_profiles) else
    sort(unique(mixing_table$profile))
  class_levels <- if (is.null(mixing_table)) seq_len(fit$n_group_classes) else
    sort(unique(mixing_table$group_class))
  mixing <- if (!is.null(mixing_table)) {
    sweep_wide(mixing_table, "group_class", "profile", "probability", class_levels, profile_levels)
  }
  group_mixing <- if (!is.null(mixing_table)) {
    mixing_table$group_class_probability[match(class_levels, mixing_table$group_class)]
  }
  continuous <- unique(as.character(profiles$indicator))
  means <- if (length(continuous)) {
    sweep_wide(profiles, "profile", "indicator", "mean", profile_levels, continuous)
  } else NULL
  covariances <- if (length(continuous)) {
    table <- sweep_table(fit, "covariances")
    lapply(profile_levels, function(k) {
      sweep_wide(table[table$profile == k, , drop = FALSE], "indicator", "indicator_2",
                 "covariance", continuous, continuous)
    })
  } else NULL
  shared <- c("profile", "indicator", "category", "probability")
  responses <- sweep_table(fit, "responses")
  responses <- if (is.data.frame(responses) && nrow(responses)) responses[shared] else NULL
  # An ordinal indicator's category probabilities enter the likelihood as any
  # categorical indicator's do; its structure is checked separately.
  if (length(fit$ordinal)) responses <- rbind(responses, sweep_table(fit, "ordinal")[shared])
  list(profile_levels = profile_levels, class_levels = class_levels, mixing = mixing,
       group_mixing = group_mixing, continuous = continuous, means = means,
       covariances = covariances, responses = responses,
       ordinal = as.character(fit$ordinal),
       counts = if (length(fit$count)) sweep_table(fit, "count_means") else NULL)
}

#' Log density of every row under every profile, from observed entries only
#' @param data The data the model was fitted to.
#' @param parameters Value of `sweep_parameters()`.
#' @param counts Optional count table (profile, indicator, mean, and dispersion for a
#'   negative binomial); defaults to the fit's own.
#' @return n x K matrix of log densities.
sweep_log_density <- function(data, parameters, counts = NULL) {
  K <- length(parameters$profile_levels)
  n <- nrow(data)
  gaussian <- if (length(parameters$continuous)) {
    x <- as.matrix(data[parameters$continuous])
    pattern <- apply(!is.na(x), 1L, function(r) paste(which(r), collapse = ","))
    out <- matrix(0, n, K)
    filled <- lapply(unique(pattern[nzchar(pattern)]), function(key) {
      rows <- which(pattern == key)
      seen <- as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
      vapply(seq_len(K), function(k) {
        root <- chol(parameters$covariances[[k]][seen, seen, drop = FALSE])
        centred <- sweep(x[rows, seen, drop = FALSE], 2L, parameters$means[k, seen])
        scaled <- backsolve(root, t(centred), transpose = TRUE)
        -0.5 * colSums(scaled^2) - sum(log(diag(root))) - 0.5 * length(seen) * log(2 * pi)
      }, numeric(length(rows)))
    })
    rows_by_pattern <- lapply(unique(pattern[nzchar(pattern)]), function(key) which(pattern == key))
    out[unlist(rows_by_pattern), ] <- do.call(rbind, lapply(filled, function(m) matrix(m, ncol = K)))
    out
  } else matrix(0, n, K)
  categorical <- if (!is.null(parameters$responses)) {
    table <- parameters$responses
    Reduce(`+`, lapply(unique(as.character(table$indicator)), function(v) {
      rows <- table[table$indicator == v, , drop = FALSE]
      observed <- as.character(data[[v]])
      vapply(parameters$profile_levels, function(k) {
        of_profile <- rows[rows$profile == k, , drop = FALSE]
        log_probability <- log(of_profile$probability[match(observed, as.character(of_profile$category))])
        ifelse(is.na(observed), 0, log_probability)
      }, numeric(n))
    }), matrix(0, n, K))
  } else matrix(0, n, K)
  counts <- counts %||% parameters$counts
  poisson <- if (!is.null(counts)) {
    Reduce(`+`, lapply(unique(as.character(counts$indicator)), function(v) {
      vapply(parameters$profile_levels, function(k) {
        row <- counts[counts$indicator == v & counts$profile == k, , drop = FALSE]
        # NB2 with variance mu + dispersion * mu^2; Poisson at dispersion zero.
        overdispersed <- "dispersion" %in% names(row) && isTRUE(row$dispersion > 0)
        log_mass <- if (overdispersed) {
          stats::dnbinom(data[[v]], mu = row$mean, size = 1 / row$dispersion, log = TRUE)
        } else stats::dpois(data[[v]], row$mean, log = TRUE)
        ifelse(is.na(data[[v]]), 0, log_mass)
      }, numeric(n))
    }), matrix(0, n, K))
  } else matrix(0, n, K)
  gaussian + categorical + poisson
}

#' Log of row sums of exponentials
sweep_log_sum_exp <- function(m) {
  top <- apply(m, 1L, max)
  top + log(rowSums(exp(m - top)))
}

#' The mixture likelihood and posteriors rebuilt from result tables
#'
#' Rows are conditionally independent given their profile; profiles are mixed
#' with class-specific shares within a group, and group classes are mixed
#' across groups (Vermunt 2003). A single-level model is the case of one
#' group class.
#' @param data Data the model was fitted to.
#' @param parameters Value of `sweep_parameters()`.
#' @param group Group of each row, or `NULL` for a single-level model.
#' @param weights Weight per row (single level), or `NULL`.
#' @param counts Optional Poisson means table.
#' @return List: `log_likelihood`, `posteriors` (n x K), `group_posteriors`.
sweep_mixture <- function(data, parameters, group = NULL, weights = NULL, counts = NULL) {
  log_density <- sweep_log_density(data, parameters, counts)
  n <- nrow(log_density)
  group <- if (is.null(group)) seq_len(n) else match(group, unique(group))
  G <- length(parameters$class_levels)
  # log p(row | group class g), then summed over the rows of each group.
  row_given_class <- vapply(seq_len(G), function(g) {
    sweep_log_sum_exp(sweep(log_density, 2L, log(parameters$mixing[g, ]), `+`))
  }, numeric(n))
  row_given_class <- matrix(row_given_class, n, G)
  group_given_class <- rowsum(row_given_class, group, reorder = TRUE)
  joint <- sweep(group_given_class, 2L, log(parameters$group_mixing), `+`)
  group_log_likelihood <- sweep_log_sum_exp(joint)
  group_posteriors <- exp(joint - group_log_likelihood)
  group_weights <- if (is.null(weights)) rep(1, length(group_log_likelihood)) else
    as.numeric(tapply(weights, group, function(w) w[[1L]]))
  posteriors <- Reduce(`+`, lapply(seq_len(G), function(g) {
    within_class <- sweep(log_density, 2L, log(parameters$mixing[g, ]), `+`)
    group_posteriors[group, g] * exp(within_class - row_given_class[, g])
  }))
  list(log_likelihood = sum(group_weights * group_log_likelihood), posteriors = posteriors,
       group_posteriors = group_posteriors)
}

#' Small moves of the free parameters, each a function of the parameters
#'
#' Every move stays inside the model: one mean, one count mean, one varying
#' dispersion, the common scale of all covariances, mass between two profile
#' shares or two group-class shares, mass between two categories of a
#' categorical indicator in one profile, and for an ordinal indicator one
#' profile's location or one shared category intercept.
#' @param parameters Value of `sweep_parameters()`.
#' @param step Relative size of a move.
#' @return List of functions mapping parameters to moved parameters.
sweep_moves <- function(parameters, step) {
  signs <- c(-1, 1)
  both <- function(move) lapply(signs, function(sign) function(p) move(p, sign * step))
  cells <- function(matrix) {
    if (is.null(matrix)) return(list())
    unname(split(which(!is.na(matrix), arr.ind = TRUE), seq_along(matrix)))
  }
  mean_moves <- unlist(lapply(cells(parameters$means), function(cell) {
    spread <- sqrt(mean(vapply(parameters$covariances, function(m) m[cell[[2L]], cell[[2L]]],
                               numeric(1))))
    both(function(p, h) {
      p$means[cell[[1L]], cell[[2L]]] <- p$means[cell[[1L]], cell[[2L]]] + h * spread
      p
    })
  }), recursive = FALSE)
  scale_moves <- if (is.null(parameters$covariances)) list() else both(function(p, h) {
    p$covariances <- lapply(p$covariances, function(m) m * (1 + h))
    p
  })
  transfer <- function(vector, from, to, h) {
    amount <- h * min(vector[[from]], vector[[to]])
    vector[[from]] <- vector[[from]] - amount
    vector[[to]] <- vector[[to]] + amount
    vector
  }
  K <- ncol(parameters$mixing)
  share_moves <- if (K < 2L) list() else unlist(lapply(seq_len(nrow(parameters$mixing)), function(g) {
    unlist(lapply(seq_len(K - 1L), function(k) both(function(p, h) {
      p$mixing[g, ] <- transfer(p$mixing[g, ], k, k + 1L, h)
      p
    })), recursive = FALSE)
  }), recursive = FALSE)
  G <- length(parameters$group_mixing)
  class_moves <- if (G < 2L) list() else unlist(lapply(seq_len(G - 1L), function(g) both(function(p, h) {
    p$group_mixing <- transfer(p$group_mixing, g, g + 1L, h)
    p
  })), recursive = FALSE)
  count_moves <- if (is.null(parameters$counts)) list() else {
    equal_dispersion <- "dispersion" %in% names(parameters$counts) &&
      all(tapply(parameters$counts$dispersion, parameters$counts$indicator,
                 function(d) diff(range(d)) < 1e-12))
    unlist(lapply(seq_len(nrow(parameters$counts)), function(i) {
      c(both(function(p, h) {
        p$counts$mean[[i]] <- p$counts$mean[[i]] * (1 + h)
        p
      }), if ("dispersion" %in% names(parameters$counts) && !equal_dispersion &&
              isTRUE(parameters$counts$dispersion[[i]] > 0)) both(function(p, h) {
        p$counts$dispersion[[i]] <- p$counts$dispersion[[i]] * (1 + h)
        p
      }))
    }), recursive = FALSE)
  }
  table <- parameters$responses
  response_moves <- if (is.null(table)) list() else {
    indicators <- unique(as.character(table$indicator))
    unlist(lapply(indicators, function(v) {
      if (v %in% parameters$ordinal) {
        categories <- unique(as.character(table$category[table$indicator == v]))
        tilt <- function(p, rows, weights) {
          tilted <- p$responses$probability[rows] * exp(weights)
          p$responses$probability[rows] <- tilted / sum(tilted)
          p
        }
        location <- unlist(lapply(parameters$profile_levels, function(k) both(function(p, h) {
          rows <- which(p$responses$indicator == v & p$responses$profile == k)
          tilt(p, rows, h * match(as.character(p$responses$category[rows]), categories))
        })), recursive = FALSE)
        intercept <- unlist(lapply(categories, function(category) both(function(p, h) {
          Reduce(function(q, k) {
            rows <- which(q$responses$indicator == v & q$responses$profile == k)
            tilt(q, rows, h * (as.character(q$responses$category[rows]) == category))
          }, parameters$profile_levels, p)
        })), recursive = FALSE)
        return(c(location, intercept))
      }
      unlist(lapply(parameters$profile_levels, function(k) {
        rows <- which(table$indicator == v & table$profile == k)
        unlist(lapply(seq_len(length(rows) - 1L), function(j) both(function(p, h) {
          p$responses$probability[rows] <- transfer(p$responses$probability[rows], j, j + 1L, h)
          p
        })), recursive = FALSE)
      }), recursive = FALSE)
    }), recursive = FALSE)
  }
  c(mean_moves, scale_moves, share_moves, class_moves, count_moves, response_moves)
}

#' The estimate is a local maximum of the rebuilt likelihood
#' @param data,parameters,group,weights As in `sweep_mixture()`.
#' @param step Relative size of each move.
#' @return One check: no move raises the likelihood by more than 1e-8 of its size.
sweep_local_maximum <- function(data, parameters, group = NULL, weights = NULL, step = 1e-3) {
  at_estimate <- sweep_mixture(data, parameters, group, weights)$log_likelihood
  moved <- vapply(sweep_moves(parameters, step), function(move) {
    sweep_mixture(data, move(parameters), group, weights)$log_likelihood
  }, numeric(1))
  sweep_at_least(sprintf("no small move of a free parameter raises the rebuilt likelihood (%d moves)",
                         length(moved)), at_estimate, max(moved), 1e-8)
}

#' The posterior table of a fit as an n x K matrix in row order
sweep_posteriors <- function(fit) {
  table <- sweep_table(fit, "posteriors")
  sweep_wide(table, "row", "profile", "posterior", sort(unique(table$row)),
             sort(unique(table$profile)))
}

#' Checks every multilpa-family fit must pass against the rebuilt likelihood
#' @param fit The fit.
#' @param data Data it was fitted to (after any scaling the caller applied).
#' @param group,weights,counts As in `sweep_mixture()`.
#' @param n_bic Sample size the fit's individual BIC should use.
#' @param empty Rows with nothing observed.
#' @param maximum Whether the estimate should be a local maximum of the likelihood
#'   (not under a prior or with held parameters).
#' @return data.frame of checks.
sweep_likelihood_checks <- function(fit, data, group = NULL, weights = NULL, counts = NULL,
                                    n_bic = nrow(data), tolerance = 1e-8,
                                    empty = rep(FALSE, nrow(data)), maximum = TRUE) {
  parameters <- sweep_parameters(fit)
  rebuilt <- sweep_mixture(data, parameters, group, weights, counts)
  criteria <- sweep_table(fit, "information_criteria")
  posteriors <- sweep_posteriors(fit)
  rows <- as.integer(rownames(posteriors))
  # A row with nothing observed keeps its place with missing posteriors (`drop_empty`).
  reported <- stats::complete.cases(posteriors)
  rbind(
    sweep_check("log likelihood rebuilt from the result tables", criteria$log_likelihood,
                rebuilt$log_likelihood, tolerance),
    sweep_true("posteriors are missing exactly for rows with nothing observed",
               identical(unname(!reported), unname(empty[rows]))),
    sweep_check("posteriors rebuilt from the result tables", posteriors[reported, , drop = FALSE],
                rebuilt$posteriors[rows[reported], , drop = FALSE], 1e-6),
    sweep_check("posterior rows sum to one", rowSums(posteriors[reported, , drop = FALSE]),
                rep(1, sum(reported)), 1e-8),
    sweep_check("profile shares sum to one", rowSums(parameters$mixing),
                rep(1, nrow(parameters$mixing)), 1e-8),
    if (!maximum) NULL else if (isTRUE(sweep_table(fit, "model")$boundary)) {
      # A fit with an estimate on its bound is a constrained maximum.
      sweep_no_reference("no small move of a free parameter raises the rebuilt likelihood",
                         "the fit is on a boundary")
    } else sweep_local_maximum(data, parameters, group, weights),
    sweep_check("AIC from the likelihood and parameter count", criteria$aic,
                -2 * criteria$log_likelihood + 2 * criteria$n_parameters, 1e-10),
    sweep_check("individual BIC from the likelihood and parameter count", criteria$bic_individual,
                -2 * criteria$log_likelihood + criteria$n_parameters * log(n_bic), 1e-10))
}

#' Checks that two fits are the same solution
#' @param label What makes the fits differ.
#' @param fit,other The two fits.
#' @param tolerance Tolerance on the likelihood; parameters use its square root.
#' @return data.frame of checks.
sweep_same_solution <- function(label, fit, other, tolerance = 1e-7) {
  first <- sweep_table(fit, "information_criteria")
  second <- sweep_table(other, "information_criteria")
  rbind(
    sweep_check(sprintf("%s: log likelihood", label), first$log_likelihood,
                second$log_likelihood, tolerance),
    sweep_check(sprintf("%s: parameter count", label), first$n_parameters,
                second$n_parameters, 0),
    sweep_check(sprintf("%s: sorted class shares", label),
                sort(sweep_parameters(fit)$mixing), sort(sweep_parameters(other)$mixing),
                sqrt(tolerance)))
}
