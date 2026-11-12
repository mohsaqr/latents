# pool_imputations(): `imputed` (a list or a mice `mids`), the number of
# imputations, `vcov_type`, `level`, `adjust`, single-level and two-level
# models. References: Rubin's rules written here and applied to the
# per-imputation table, mitools::MIcombine on the same estimates, the
# single-fit inference when every imputation is the same data, the order of
# the imputations, and classed refusals.

sweep_module_pooling <- function() {
  blobs <- sweep_continuous()$blobs
  nested <- sweep_two_level()
  #' Completed data sets: the data with noise added to one indicator
  completed <- function(data, column, m, seed = 61L) {
    sweep_seeded(seed, function() lapply(seq_len(m), function(i) {
      data[[column]] <- data[[column]] + stats::rnorm(nrow(data), sd = 0.25)
      data
    }))
  }
  models <- list(
    single = list(data = blobs, vars = c("a", "b"), id = NULL, column = "a",
                  arguments = list(n_profiles = 2L)),
    two_level = list(data = nested, vars = c("y1", "y2"), id = "id", column = "y1",
                     arguments = list(n_profiles = 2L, n_group_classes = 2L)))
  pool <- function(model, imputed, ...) {
    do.call(pool_imputations, c(list(imputed, model$vars, model$id), model$arguments,
                                list(n_starts = 4L, seed = 5L, ...)))
  }
  key <- function(table) paste(table$level, table$outcome, table$term, table$parameter)
  #' Rubin (1987): the pooled estimate, its variance components and the
  #' reference t distribution, per parameter
  rubin <- function(imputations, level) {
    parts <- split(imputations, key(imputations))
    out <- do.call(rbind, lapply(names(parts), function(name) {
      part <- parts[[name]]
      m <- nrow(part)
      within <- mean(part$standard_error^2)
      between <- stats::var(part$estimate)
      total <- within + (1 + 1 / m) * between
      riv <- if (within > 0) (1 + 1 / m) * between / within else 0
      df <- if (riv > 0) (m - 1) * (1 + 1 / riv)^2 else Inf
      data.frame(key = name, estimate = mean(part$estimate), standard_error = sqrt(total),
                 within = within, between = between, riv = riv, df = df,
                 fmi = if (riv > 0) (riv + 2 / (df + 3)) / (riv + 1) else 0)
    }))
    out$conf_low <- out$estimate - stats::qt(1 - (1 - level) / 2, out$df) * out$standard_error
    out$conf_high <- out$estimate + stats::qt(1 - (1 - level) / 2, out$df) * out$standard_error
    out
  }
  pooled_checks <- function(pooled, level, adjust) {
    estimates <- as.data.frame(pooled)
    imputations <- sweep_table(pooled, "imputations")
    expected <- rubin(imputations, level)
    expected <- expected[match(key(estimates), expected$key), ]
    finite <- is.finite(estimates$df)
    tested <- !is.na(estimates$p_value)
    combined <- lapply(split(imputations, key(imputations)), function(part) {
      mitools::MIcombine(as.list(part$estimate), as.list(part$standard_error^2))
    })[key(estimates)]
    rbind(
      sweep_true("one pooled row per parameter of the imputations",
                 setequal(key(estimates), key(imputations)) && !anyDuplicated(key(estimates))),
      sweep_check("Rubin's rules: estimates", estimates$estimate, expected$estimate, 1e-10),
      sweep_check("Rubin's rules: standard errors", estimates$standard_error,
                  expected$standard_error, 1e-10),
      sweep_check("Rubin's rules: within and between variances",
                  c(estimates$within, estimates$between), c(expected$within, expected$between),
                  1e-10),
      sweep_check("Rubin's rules: relative increase in variance", estimates$riv, expected$riv,
                  1e-8),
      sweep_check("Rubin's rules: degrees of freedom", estimates$df[finite], expected$df[finite],
                  1e-8),
      sweep_true("no between-imputation variance gives infinite degrees of freedom",
                 identical(finite, is.finite(expected$df))),
      sweep_check("Rubin's rules: fraction of missing information", estimates$fmi, expected$fmi,
                  1e-8),
      sweep_check("intervals from the t reference distribution",
                  c(estimates$conf_low, estimates$conf_high),
                  c(expected$conf_low, expected$conf_high), 1e-8),
      sweep_check("p-values from the t reference distribution", estimates$p_value[tested],
                  2 * stats::pt(-abs(estimates$estimate / estimates$standard_error),
                                estimates$df)[tested], 1e-8),
      sweep_check(sprintf("adjusted p-values are p.adjust(method = \"%s\") of the tested ones",
                          adjust), estimates$p_adjusted[tested],
                  stats::p.adjust(estimates$p_value[tested], adjust), 1e-10),
      sweep_check("mitools::MIcombine: estimates", estimates$estimate,
                  vapply(combined, function(x) unname(x$coefficients), numeric(1)), 1e-10),
      sweep_check("mitools::MIcombine: standard errors", estimates$standard_error,
                  vapply(combined, function(x) sqrt(as.numeric(x$variance)), numeric(1)), 1e-10),
      sweep_check("mitools::MIcombine: degrees of freedom", estimates$df[finite],
                  vapply(combined, function(x) unname(x$df), numeric(1))[finite], 1e-8))
  }

  grid <- sweep_grid(model = names(models), m = c(2L, 5L),
                     vcov_type = c("observed", "robust", "opg"), level = c(0.9, 0.95),
                     adjust = c("none", "BH"))
  rule_cases <- lapply(grid, function(option) {
    sweep_define(sprintf("rubin[model=%s,m=%d,vcov=%s,level=%s,adjust=%s]", option$model, option$m,
                         option$vcov_type, format(option$level), option$adjust), function() {
      model <- models[[option$model]]
      pooled <- pool(model, completed(model$data, model$column, option$m),
                     vcov_type = option$vcov_type, level = option$level, adjust = option$adjust)
      fits <- sweep_table(pooled, "fits")
      rbind(pooled_checks(pooled, option$level, option$adjust),
            sweep_check("one fit row per imputation", nrow(fits), option$m, 0),
            sweep_true("every imputation converged", all(fits$converged)))
    })
  })

  invariance_cases <- lapply(names(models), function(name) {
    sweep_define(sprintf("invariances[model=%s]", name), function() {
      model <- models[[name]]
      sets <- completed(model$data, model$column, 4L)
      base <- as.data.frame(pool(model, sets))
      reversed <- as.data.frame(pool(model, rev(sets)))
      copies <- as.data.frame(pool(model, list(model$data, model$data, model$data)))
      single <- do.call(multilpa, c(list(model$data, model$vars, model$id), model$arguments,
                                    list(n_starts = 4L, seed = 5L)))
      inference <- as.data.frame(parameter_inference(single))
      holm <- as.data.frame(pool(model, sets, adjust = "holm"))
      tested <- !is.na(holm$p_value)
      rbind(
        sweep_check("imputations in reverse order: estimates",
                    reversed$estimate[match(key(base), key(reversed))], base$estimate, 1e-8),
        sweep_check("imputations in reverse order: standard errors",
                    reversed$standard_error[match(key(base), key(reversed))],
                    base$standard_error, 1e-8),
        sweep_check("three copies of one data set: no between-imputation variance",
                    copies$between, rep(0, nrow(copies)), 1e-12),
        sweep_check("three copies of one data set: the single fit's estimates",
                    copies$estimate, inference$estimate[match(key(copies), key(inference))],
                    1e-8),
        sweep_check("three copies of one data set: the single fit's standard errors",
                    copies$standard_error,
                    inference$standard_error[match(key(copies), key(inference))], 1e-8),
        sweep_check("adjust = \"holm\" is p.adjust()", holm$p_adjusted[tested],
                    stats::p.adjust(holm$p_value[tested], "holm"), 1e-10))
    })
  })

  mids_case <- sweep_define("mice-mids-input", function() {
    model <- models$single
    holed <- sweep_punch(model$data[model$vars], 0.1)
    imputed <- sweep_seeded(8L, function() mice::mice(holed, m = 3L, printFlag = FALSE))
    from_mids <- as.data.frame(pool(model, imputed))
    from_list <- as.data.frame(pool(model, lapply(1:3, function(i) mice::complete(imputed, i))))
    rbind(
      sweep_check("a mids object against its completed data sets: estimates",
                  from_mids$estimate, from_list$estimate, 1e-10),
      sweep_check("a mids object against its completed data sets: standard errors",
                  from_mids$standard_error, from_list$standard_error, 1e-10),
      pooled_checks(pool(model, imputed), 0.95, "none"))
  })

  refusal_case <- sweep_define("refusals", function() {
    model <- models$single
    sets <- completed(model$data, model$column, 2L)
    rbind(
      sweep_refuses("one data frame, not a list of them", function() pool(model, model$data)),
      sweep_refuses("a single completed data set", function() pool(model, sets[1L])),
      sweep_refuses("completed data sets of different sizes",
                    function() pool(model, list(sets[[1L]], sets[[2L]][-1L, ]))),
      sweep_refuses("completed data sets with different columns",
                    function() pool(model, list(sets[[1L]], sets[[2L]]["a"]))),
      sweep_refuses("a completed data set with a missing value",
                    function() pool(model, list(sets[[1L]], sweep_punch(sets[[2L]], 0.05)))),
      sweep_refuses("`level` outside (0, 1)", function() pool(model, sets, level = 1.2)),
      sweep_refuses("an unknown `vcov_type`", function() pool(model, sets, vcov_type = "hc3")),
      sweep_refuses("an unknown `adjust`", function() pool(model, sets, adjust = "sidak")))
  })
  c(rule_cases, invariance_cases, list(mids_case, refusal_case))
}
