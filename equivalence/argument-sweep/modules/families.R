# multilpa(family = ): the three group-class families (additive, dispersion,
# additive_dispersion) and the two cross-level families, across
# `n_group_classes`, `between_variance`, `variance_model`, `weights`, `search`
# and the input types. References: the likelihood of each family written here
# from its definition and evaluated at the fit's result tables, with group
# posteriors and a local-maximum check; lme4's random-intercept model at one
# group class; parameter counts and criteria by formula; invariances; classed
# refusals.

sweep_module_families <- function() {
  simulated <- sweep_seeded(909L, function() {
    groups <- 50L
    size <- sample(5:10, groups, TRUE)
    class <- rep(1:2, each = groups / 2L)
    id <- rep(seq_len(groups), size)
    centre <- rbind(c(0, 0), c(2.5, -1.5))
    intercepts <- centre[class, ] + matrix(stats::rnorm(groups * 2L, sd = 0.6), groups)
    scale <- c(1, 1.8)[class[id]]
    values <- intercepts[id, ] + matrix(stats::rnorm(length(id) * 2L), length(id)) * scale
    cbind(id = id, stats::setNames(as.data.frame(values), c("y1", "y2")))
  })
  sets <- list(
    simulated = list(data = simulated, vars = c("y1", "y2"), id = "id"),
    two_level = list(data = sweep_two_level(), vars = c("y1", "y2", "y3"), id = "id"),
    engagement = list(data = subset(course_engagement, student <= 60),
                      vars = c("browse", "lectures", "forum_read"), id = "student"))
  fit <- function(set, family, ..., data = set$data, id = set$id, n_starts = 8L) {
    multilpa(data, set$vars, id, family = family, n_starts = n_starts, seed = 17L,
             tol = 1e-11, max_iter = 20000L, ...)
  }
  model_row <- function(model) sweep_table(model, "fit")
  normalised <- function(weights) weights * sum(weights > 0) / sum(weights)
  group_frames <- function(set, data = set$data, id = set$id) {
    lapply(split(data[set$vars], data[[id]]), as.matrix)
  }
  #' No one-at-a-time move of a parameter raises a likelihood
  #' @param parameters Named list of numeric arrays.
  #' @param likelihood Function of such a list.
  #' @param scales,simplexes Names of the positive and of the probability blocks.
  local_maximum <- function(parameters, likelihood, scales, simplexes, step = 1e-3) {
    at_estimate <- likelihood(parameters)
    plain <- setdiff(names(parameters), simplexes)
    # A column whose classes share one value is one parameter (an "equal"
    # constraint) and moves as one; otherwise each cell moves alone.
    units <- function(block) {
      cells <- matrix(seq_along(block), nrow = NROW(block))
      tied <- nrow(cells) > 1L & apply(matrix(block, nrow = nrow(cells)), 2L, function(column) {
        isTRUE(all.equal(min(column), max(column), tolerance = 1e-12))
      })
      c(lapply(which(tied), function(j) cells[, j]), as.list(cells[, !tied]))
    }
    moved <- unlist(lapply(plain, function(name) {
      block <- parameters[[name]]
      unlist(lapply(units(block), function(cells) {
        vapply(c(-1, 1), function(sign) {
          trial <- parameters
          size <- step * max(abs(block[cells]), 0.1)
          trial[[name]][cells] <- block[cells] + sign * size
          if (name %in% scales && any(trial[[name]][cells] <= 0)) return(-Inf)
          likelihood(trial)
        }, numeric(1))
      }))
    }))
    shifted <- unlist(lapply(simplexes, function(name) {
      block <- parameters[[name]]
      rows <- if (is.matrix(block)) seq_len(nrow(block)) else 1L
      unlist(lapply(rows, function(r) {
        shares <- if (is.matrix(block)) block[r, ] else block
        if (length(shares) < 2L) return(numeric())
        vapply(c(-1, 1), function(sign) {
          size <- sign * step * min(shares[1:2])
          shares[1:2] <- shares[1:2] + c(size, -size)
          trial <- parameters
          if (is.matrix(block)) trial[[name]][r, ] <- shares else trial[[name]] <- shares
          likelihood(trial)
        }, numeric(1))
      }))
    }))
    sweep_at_least(sprintf("no small move of a free parameter raises the rebuilt likelihood (%d moves)",
                           length(moved) + length(shifted)),
                   at_estimate, max(moved, shifted), 1e-8)
  }

  # ---- group-class families --------------------------------------------------
  #' Parameters of a group-class fit from its tables: classes x indicators
  additive_parameters <- function(model, vars) {
    table <- sweep_table(model, "parameters")
    classes <- sweep_table(model, "group_classes")
    block <- function(level, parameter) {
      rows <- table[table$level == level & table$parameter == parameter, ]
      matrix(vapply(vars, function(v) vapply(classes$group_class, function(h) {
        value <- rows$estimate[rows$indicator == v & rows$group_class %in% c(h, "shared")]
        stopifnot("one estimate per class and indicator" = length(value) == 1L)
        value
      }, numeric(1)), numeric(nrow(classes))), nrow = nrow(classes))
    }
    list(means = block("between", "mean"), between = block("between", "variance"),
         within = block("within", "variance"), shares = classes$weight)
  }
  #' A group's ratings of an indicator are exchangeable Gaussian: variance
  #' W + B, covariance B. Its log density from the group's size, mean and
  #' within-group sum of squares.
  additive_group_density <- function(frames, p) {
    t(vapply(frames, function(values) {
      n <- nrow(values)
      average <- colMeans(values)
      scatter <- colSums(sweep(values, 2L, average)^2)
      vapply(seq_along(p$shares), function(h) {
        total <- p$within[h, ] + n * p$between[h, ]
        sum(-n / 2 * log(2 * pi) - (n - 1) / 2 * log(p$within[h, ]) - log(total) / 2 -
              scatter / (2 * p$within[h, ]) - n * (average - p$means[h, ])^2 / (2 * total))
      }, numeric(1))
    }, numeric(length(p$shares))))
  }
  additive_mixture <- function(frames, p, weights = NULL) {
    density <- additive_group_density(frames, p)
    if (length(p$shares) == 1L) density <- matrix(density, ncol = 1L)
    joint <- sweep(density, 2L, log(p$shares), `+`)
    by_group <- sweep_log_sum_exp(joint)
    list(log_likelihood = sum((weights %||% 1) * by_group), posteriors = exp(joint - by_group))
  }
  additive_count <- function(model, classes, indicators, family) {
    varying <- identical(model_row(model)$between_variance, "varying")
    (classes - 1L) +
      indicators * (if (family == "dispersion") 1L else classes) +
      indicators * (if (family == "additive") 1L else classes) +
      indicators * (if (varying) classes else 1L)
  }
  additive_checks <- function(model, set, family, classes, data = set$data, id = set$id,
                              weights = NULL) {
    frames <- group_frames(set, data, id)
    p <- additive_parameters(model, set$vars)
    rebuilt <- additive_mixture(frames, p, weights)
    groups <- sweep_table(model, "groups")
    reported <- as.matrix(groups[grep("^probability_group_class_", names(groups))])
    order <- match(names(frames), as.character(groups[[id]]))
    row <- model_row(model)
    rbind(
      sweep_check("log likelihood rebuilt from the result tables", row$log_likelihood,
                  rebuilt$log_likelihood, 1e-8),
      sweep_check("group posteriors rebuilt from the result tables", reported[order, , drop = FALSE],
                  rebuilt$posteriors, 1e-6),
      sweep_check("group-class shares sum to one", sum(p$shares), 1, 1e-8),
      sweep_check("parameter count", row$n_parameters,
                  additive_count(model, classes, length(set$vars), family), 0),
      sweep_check("AIC from the likelihood and parameter count", row$aic,
                  -2 * row$log_likelihood + 2 * row$n_parameters, 1e-10),
      sweep_check("group BIC from the likelihood and parameter count", row$bic,
                  -2 * row$log_likelihood + row$n_parameters * log(row$n_groups), 1e-10),
      sweep_check("individual BIC from the likelihood and parameter count", row$bic_individual,
                  -2 * row$log_likelihood + row$n_parameters * log(row$n_obs), 1e-10),
      if (isTRUE(row$boundary)) {
        sweep_no_reference("no small move of a free parameter raises the rebuilt likelihood",
                           "the fit is on a boundary")
      } else {
        local_maximum(p, function(trial) additive_mixture(frames, trial, weights)$log_likelihood,
                      scales = c("between", "within"), simplexes = "shares")
      })
  }
  additive_grid <- sweep_grid(set = names(sets),
                              family = c("additive", "dispersion", "additive_dispersion"),
                              classes = 1:3, between = c("default", "equal", "varying"))
  additive_cases <- lapply(additive_grid, function(option) {
    sweep_define(sprintf("group-classes[data=%s,family=%s,classes=%d,between=%s]", option$set,
                         option$family, option$classes, option$between), function() {
      set <- sets[[option$set]]
      if (identical(option$family, "dispersion") && identical(option$between, "varying")) {
        # The dispersion family holds the between-group variances equal.
        return(sweep_refuses("between_variance = \"varying\" in the dispersion family",
                             function() fit(set, option$family,
                                            n_group_classes = option$classes,
                                            between_variance = "varying")))
      }
      model <- if (identical(option$between, "default")) {
        fit(set, option$family, n_group_classes = option$classes)
      } else {
        fit(set, option$family, n_group_classes = option$classes,
            between_variance = option$between)
      }
      rows <- additive_checks(model, set, option$family, option$classes)
      if (option$classes > 1L) return(rows)
      # One class: a random-intercept model per indicator.
      reference <- lapply(set$vars, function(v) {
        lme4::lmer(stats::as.formula(sprintf("%s ~ 1 + (1 | %s)", v, set$id)), set$data,
                   REML = FALSE, control = lme4::lmerControl(
                     optCtrl = list(xtol_abs = 1e-10, ftol_abs = 1e-12)))
      })
      p <- additive_parameters(model, set$vars)
      rbind(rows,
            sweep_check("one class: lme4's log likelihood", model_row(model)$log_likelihood,
                        sum(vapply(reference, function(m) as.numeric(stats::logLik(m)),
                                   numeric(1))), 1e-7),
            sweep_check("one class: lme4's means", as.numeric(p$means),
                        vapply(reference, function(m) unname(lme4::fixef(m)), numeric(1)), 1e-4),
            sweep_check("one class: lme4's within-group variances", as.numeric(p$within),
                        vapply(reference, function(m) stats::sigma(m)^2, numeric(1)), 1e-4),
            sweep_check("one class: lme4's between-group variances", as.numeric(p$between),
                        vapply(reference, function(m) as.numeric(lme4::VarCorr(m)[[1L]]),
                               numeric(1)), 1e-3))
    })
  })

  # ---- cross-level families --------------------------------------------------
  cross_parameters <- function(model, vars, family) {
    profiles <- sweep_table(model, "profiles")
    level <- function(name, value) {
      rows <- profiles[profiles$level == name, ]
      sweep_wide(rows, "profile", "indicator", value, unique(rows$profile), vars)
    }
    individual <- profiles[profiles$level == "individual", ]
    group <- profiles[profiles$level == "group", ]
    composition <- sweep_table(model, "composition")
    prevalence <- individual$prevalence[match(unique(individual$profile), individual$profile)]
    list(means = level("individual", "mean"), variances = level("individual", "variance"),
         group_means = level("group", "mean"), group_variances = level("group", "variance"),
         shares = group$prevalence[match(unique(group$profile), group$profile)],
         composition = if (identical(family, "full_cross_level")) {
           sweep_wide(composition, "group_class", "profile", "share",
                      unique(group$profile), unique(individual$profile))
         } else matrix(prevalence, nrow = 1L))
  }
  gaussian <- function(values, means, variances) {
    matrix(vapply(seq_len(nrow(means)), function(k) {
      rowSums(stats::dnorm(values, matrix(means[k, ], nrow(values), ncol(values), byrow = TRUE),
                           matrix(sqrt(variances[k, ]), nrow(values), ncol(values), byrow = TRUE),
                           log = TRUE))
    }, numeric(nrow(values))), nrow = nrow(values))
  }
  #' L_j = sum_h omega_h f(m_j | h) prod_i sum_k pi_{k|h} f(x_ij | k), with
  #' pi_{k|h} = pi_k in the restricted family (Houle, Morin & Harvey 2026).
  cross_mixture <- function(frames, p) {
    averages <- t(vapply(frames, colMeans, numeric(ncol(frames[[1L]]))))
    between <- gaussian(averages, p$group_means, p$group_variances)
    classes <- length(p$shares)
    ratings <- t(vapply(frames, function(values) {
      density <- gaussian(values, p$means, p$variances)
      vapply(seq_len(classes), function(h) {
        shares <- p$composition[min(h, nrow(p$composition)), ]
        sum(sweep_log_sum_exp(sweep(density, 2L, log(shares), `+`)))
      }, numeric(1))
    }, numeric(classes)))
    if (classes == 1L) ratings <- matrix(ratings, ncol = 1L)
    joint <- sweep(between + ratings, 2L, log(p$shares), `+`)
    by_group <- sweep_log_sum_exp(joint)
    list(log_likelihood = sum(by_group), posteriors = exp(joint - by_group))
  }
  cross_grid <- sweep_grid(set = names(sets),
                           family = c("restricted_cross_level", "full_cross_level"),
                           profiles = 2:3, classes = 2L,
                           variance = c("varying", "equal"), between = c("varying", "equal"))
  cross_cases <- lapply(cross_grid, function(option) {
    sweep_define(sprintf("cross-level[data=%s,family=%s,profiles=%d,variance=%s,between=%s]",
                         option$set, option$family, option$profiles, option$variance,
                         option$between), function() {
      set <- sets[[option$set]]
      model <- fit(set, option$family, n_profiles = option$profiles,
                   n_group_classes = option$classes, variance_model = option$variance,
                   between_variance = option$between)
      frames <- group_frames(set)
      p <- cross_parameters(model, set$vars, option$family)
      rebuilt <- cross_mixture(frames, p)
      groups <- sweep_table(model, "groups")
      reported <- as.matrix(groups[grep("^probability_group_class_", names(groups))])
      row <- model_row(model)
      d <- length(set$vars)
      full <- identical(option$family, "full_cross_level")
      expected <- option$profiles * d + d * (if (option$variance == "varying") option$profiles else 1L) +
        option$classes * d + d * (if (option$between == "varying") option$classes else 1L) +
        (option$classes - 1L) + (option$profiles - 1L) * (if (full) option$classes else 1L)
      rbind(
        sweep_check("log likelihood rebuilt from the result tables", row$log_likelihood,
                    rebuilt$log_likelihood, 1e-8),
        sweep_check("group posteriors rebuilt from the result tables",
                    reported[match(names(frames), as.character(groups[[set$id]])), , drop = FALSE],
                    rebuilt$posteriors, 1e-6),
        sweep_check("composition rows sum to one", rowSums(p$composition),
                    rep(1, nrow(p$composition)), 1e-8),
        sweep_check("parameter count", row$n_parameters, expected, 0),
        sweep_check("AIC from the likelihood and parameter count", row$aic,
                    -2 * row$log_likelihood + 2 * row$n_parameters, 1e-10),
        sweep_check("group BIC from the likelihood and parameter count", row$bic,
                    -2 * row$log_likelihood + row$n_parameters * log(row$n_groups), 1e-10),
        if (isTRUE(row$boundary)) {
          sweep_no_reference("no small move of a free parameter raises the rebuilt likelihood",
                             "the fit is on a boundary")
        } else {
          local_maximum(p, function(trial) cross_mixture(frames, trial)$log_likelihood,
                        scales = c("variances", "group_variances"),
                        simplexes = c("shares", "composition"))
        })
    })
  })

  # ---- invariances -----------------------------------------------------------
  families <- list(
    additive = list(n_group_classes = 2L),
    dispersion = list(n_group_classes = 2L),
    additive_dispersion = list(n_group_classes = 2L),
    restricted_cross_level = list(n_profiles = 2L, n_group_classes = 2L),
    full_cross_level = list(n_profiles = 2L, n_group_classes = 2L))
  refit <- function(set, family, ...) do.call(fit, c(list(set, family), families[[family]], list(...)))
  same <- function(label, first, second, tolerance = 1e-7) {
    rbind(
      sweep_check(sprintf("%s: log likelihood", label), model_row(first)$log_likelihood,
                  model_row(second)$log_likelihood, tolerance),
      sweep_check(sprintf("%s: parameter count", label), model_row(first)$n_parameters,
                  model_row(second)$n_parameters, 0),
      sweep_check(sprintf("%s: sorted group-class shares", label),
                  sort(sweep_table(first, "group_classes")$weight),
                  sort(sweep_table(second, "group_classes")$weight), sqrt(tolerance)))
  }
  invariance <- function(name, code, among = names(families)) {
    lapply(sweep_grid(set = c("simulated", "two_level"), family = among), function(option) {
      sweep_define(sprintf("%s[data=%s,family=%s]", name, option$set, option$family), function() {
        code(sets[[option$set]], option$family)
      })
    })
  }
  group_class_families <- c("additive", "dispersion", "additive_dispersion")
  invariance_cases <- c(
    invariance("row-order", function(set, family) {
      shuffled <- sweep_seeded(4L, function() set$data[sample(nrow(set$data)), ])
      same("rows shuffled", refit(set, family), refit(set, family, data = shuffled))
    }),
    invariance("identifier-type", function(set, family) {
      base <- refit(set, family)
      text <- transform(set$data, key = sprintf("g%03d", set$data[[set$id]]))
      rbind(same("character identifier", base, refit(set, family, data = text, id = "key")),
            same("factor identifier", base,
                 refit(set, family, data = transform(text, key = factor(key)), id = "key")))
    }),
    invariance("seed", function(set, family) {
      same("the same seed twice", refit(set, family), refit(set, family), 0)
    }),
    invariance("search", function(set, family) {
      same("search = emEM against kmeans, 20 starts each",
           refit(set, family, n_starts = 20L, search = "kmeans"),
           refit(set, family, n_starts = 20L, search = "emEM"), 1e-6)
    }),
    invariance("unit-weights", function(set, family) {
      same("group weights of one", refit(set, family),
           refit(set, family, data = cbind(set$data, w = 1), weights = "w"), 1e-8)
    }),
    invariance("units", function(set, family) {
      scaled <- set$data
      scaled[set$vars] <- scaled[set$vars] * 10
      base <- refit(set, family)
      other <- refit(set, family, data = scaled)
      # Cross-level families also count each group's mean as an observation.
      values <- nrow(set$data) * length(set$vars) +
        if (family %in% group_class_families) 0L else
          length(unique(set$data[[set$id]])) * length(set$vars)
      rbind(
        sweep_check("indicators times ten: log likelihood shifts by the Jacobian",
                    model_row(other)$log_likelihood,
                    model_row(base)$log_likelihood - values * log(10), 1e-6),
        sweep_check("indicators times ten: sorted group-class shares",
                    sort(sweep_table(other, "group_classes")$weight),
                    sort(sweep_table(base, "group_classes")$weight), 1e-3))
    }),
    invariance("integer-group-weights", function(set, family) {
      ids <- unique(set$data[[set$id]])
      w <- sweep_seeded(9L, function() sample(1:3, length(ids), TRUE))
      data <- cbind(set$data, w = w[match(set$data[[set$id]], ids)])
      copies <- do.call(rbind, lapply(1:3, function(copy) {
        rows <- data[data$w >= copy, , drop = FALSE]
        rows$key <- (copy - 1L) * max(ids) + rows[[set$id]]
        rows
      }))
      weighted <- refit(set, family, data = data, weights = "w", n_starts = 20L)
      expanded <- refit(set, family, data = copies, id = "key", n_starts = 20L)
      weights <- normalised(w[match(names(group_frames(set, data)), as.character(ids))])
      rbind(
        sweep_check("integer group weights against repeated groups: pseudo log likelihood",
                    model_row(weighted)$log_likelihood,
                    length(ids) / sum(w) * model_row(expanded)$log_likelihood, 1e-6),
        additive_checks(weighted, set, family, 2L, data = data, weights = weights))
    }, among = group_class_families))

  # ---- refusals --------------------------------------------------------------
  set <- sets$simulated
  refusal_cases <- list(sweep_define("refusals", function() {
    rbind(
      sweep_refuses("an unknown family",
                    function() multilpa(set$data, set$vars, "id", n_group_classes = 2L,
                                        family = "multiplicative")),
      sweep_refuses("a group-class family without `id`",
                    function() multilpa(set$data, set$vars, n_group_classes = 2L,
                                        family = "additive")),
      sweep_refuses("a group-class family with `n_profiles`",
                    function() multilpa(set$data, set$vars, "id", n_profiles = 2L,
                                        n_group_classes = 2L, family = "additive")),
      sweep_refuses("a group-class family with `categorical`",
                    function() multilpa(transform(set$data, c1 = factor(y1 > 0)),
                                        c(set$vars, "c1"), "id", n_group_classes = 2L,
                                        family = "dispersion", categorical = "c1")),
      sweep_refuses("an unknown `between_variance`",
                    function() multilpa(set$data, set$vars, "id", n_group_classes = 2L,
                                        family = "additive", between_variance = "free")),
      sweep_refuses("`between_variance` in the profile model",
                    function() multilpa(set$data, set$vars, "id", n_profiles = 2L,
                                        n_group_classes = 2L, between_variance = "equal")),
      sweep_refuses("a missing indicator value in a group-class family",
                    function() multilpa(sweep_punch(set$data, 0.05, columns = set$vars), set$vars,
                                        "id", n_group_classes = 2L, family = "additive")),
      sweep_refuses("zero group classes",
                    function() multilpa(set$data, set$vars, "id", n_group_classes = 0L,
                                        family = "additive")),
      sweep_refuses("a cross-level family with `missing`",
                    function() multilpa(set$data, set$vars, "id", n_profiles = 2L,
                                        n_group_classes = 2L, family = "full_cross_level",
                                        missing = "fiml")),
      sweep_refuses("an unknown cross-level `variance_model`",
                    function() multilpa(set$data, set$vars, "id", n_profiles = 2L,
                                        n_group_classes = 2L, family = "full_cross_level",
                                        variance_model = "free")),
      sweep_refuses("weights varying within a group",
                    function() multilpa(transform(set$data, w = seq_along(id)), set$vars, "id",
                                        n_group_classes = 2L, family = "additive",
                                        weights = "w")),
      sweep_refuses("an unknown search",
                    function() multilpa(set$data, set$vars, "id", n_group_classes = 2L,
                                        family = "additive", search = "random")))
  }))
  c(additive_cases, cross_cases, invariance_cases, refusal_cases)
}
