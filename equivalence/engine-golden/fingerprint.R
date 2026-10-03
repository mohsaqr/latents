# Fingerprints of fitted models: everything a user can observe about a fit,
# reduced to plain data that can be saved and compared across sessions.
#
# A fingerprint never holds a function, an environment or a formula with its
# environment: those differ between sessions even when nothing has changed.
# Formulas and calls are stored deparsed.

# Fit classes, in the order the package defines them. Each gets the full set of
# fit probes (print, summary, every table, logLik, coef, vcov, inference ...).
.golden_fit_classes <- c(
  "multilpa", "multilpa_covariates", "multilpa_transitions", "multilpa_lta",
  "multilpa_additive", "multilpa_cross_level", "latents_mixture_regression",
  "latents_growth_mixture"
)

# Classes whose get_results() takes `vcov_type`: their tables are recorded
# again under the robust sandwich, which is a second inference path.
.golden_vcov_table_classes <- c(
  "multilpa_additive", "multilpa_lta", "latents_mixture_regression",
  "latents_growth_mixture"
)

#' Reduce any object to plain, session-independent data
#'
#' Recurses through lists and attributes. Functions, environments, external
#' pointers and plot objects become placeholders; formulas and calls become
#' their deparsed text.
#' @param x Any R object.
#' @return A plain object of atomic vectors and lists.
golden_clean <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, c("ggplot", "gg", "gtable", "grob"))) return("<plot>")
  if (is.function(x)) return("<function>")
  if (is.environment(x)) return("<environment>")
  if (inherits(x, "formula") || is.language(x)) {
    return(paste0("<language> ", paste(deparse(unclass(x), width.cutoff = 500L),
                                       collapse = " ")))
  }
  if (typeof(x) %in% c("externalptr", "S4", "weakref", "bytecode", "promise",
                       "char", "...", "any")) {
    return(sprintf("<%s>", typeof(x)))
  }
  kept <- attributes(x)
  value <- x
  if (is.list(value)) {
    value <- lapply(unclass(value), golden_clean)
  }
  attributes(value) <- NULL
  # `names` and `row.names` are plain; any other attribute may itself hold
  # matrices, formulas or environments.
  if (!is.null(kept)) attributes(value) <- lapply(kept, golden_clean)
  value
}

#' Is this condition one of the package's own (classed) refusals?
#' @param condition A condition object.
#' @return A single logical.
golden_is_refusal <- function(condition) {
  any(grepl("^(latents|multilpa)_", class(condition)))
}

#' Evaluate one probe, recording its value, its conditions and any error
#'
#' Warnings and messages are recorded and muffled; an error ends the probe and
#' is recorded as `list(error_class =, message =, expected =)`, where
#' `expected` is `TRUE` for one of the package's classed refusals and `FALSE`
#' for anything else (an unclassed error, an internal failure), which the run
#' report lists for a human to look at.
#' @param expr The probe, evaluated lazily.
#' @param clean Whether to reduce the value with `golden_clean()`.
#' @return `list(value =, conditions =)`; `conditions` is a list of
#'   `list(type, class, message)`.
golden_probe <- function(expr, clean = TRUE) {
  conditions <- list()
  record <- function(type, condition) {
    conditions[[length(conditions) + 1L]] <<- list(
      type = type, class = class(condition),
      message = conditionMessage(condition))
  }
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(list(error_class = class(e), message = conditionMessage(e),
                     expected = golden_is_refusal(e)),
                class = "golden_error")
    }),
    warning = function(w) {
      record("warning", w)
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      record("message", m)
      invokeRestart("muffleMessage")
    }
  )
  error <- inherits(value, "golden_error")
  if (error) value <- unclass(value)
  if (clean && !error) value <- golden_clean(value)
  list(value = value, error = error, conditions = conditions)
}

#' Does a package method exist for this generic and any class of `x`?
#' @param generic Generic name.
#' @param x Object.
#' @return A single logical.
golden_has_method <- function(generic, x) {
  namespace <- asNamespace("latents")
  any(vapply(class(x), function(cls) {
    method <- paste(generic, cls, sep = ".")
    exists(method, envir = namespace, inherits = FALSE) ||
      !is.null(utils::getS3method(generic, cls, optional = TRUE,
                                  envir = namespace))
  }, logical(1L)))
}

#' Which `what` values to record for a non-fit object with a get_results()
#' @param x Object.
#' @return Character vector of table names.
golden_whats <- function(x) {
  if (inherits(x, c("latents_family_enumeration",
                    "latents_transition_enumeration"))) {
    return(c("candidates", "best"))
  }
  if (inherits(x, "latents_regression_enumeration")) return("fit")
  if (inherits(x, "latents_pooled")) {
    return(c("estimates", "imputations", "fits"))
  }
  "all"
}

#' The table names a fit's get_results() offers (without "all")
#' @param x A fit.
#' @return Character vector.
golden_table_names <- function(x) {
  namespace <- asNamespace("latents")
  if (inherits(x, c("multilpa", "multilpa_covariates",
                    "multilpa_transitions"))) {
    return(names(namespace$.multilpa_catalogue(x)))
  }
  method <- utils::getS3method("get_results", class(x)[1L], envir = namespace)
  setdiff(eval(formals(method)$what), c("all", "recovery"))
}

#' Run every probe of a fitted model
#' @param x A fit of one of `.golden_fit_classes`.
#' @param seed Seed set before each probe, so no probe's random draws depend
#'   on what earlier probes consumed.
#' @return A named list of probe results.
golden_fit_probes <- function(x, seed) {
  probe <- function(expr) {
    set.seed(seed)
    golden_probe(expr)
  }
  has <- function(generic) golden_has_method(generic, x)
  result <- list(
    class = class(x),
    print = probe(utils::capture.output(print(x))),
    summary = probe(utils::capture.output(print(summary(x)))),
    tables = probe(get_results(x, "all")),
    as_data_frame = probe(as.data.frame(x))
  )
  # When "all" fails on one table, the rest would go unrecorded; record each
  # table of the catalogue on its own instead, so the failure stays confined
  # to the table that has it.
  if (isTRUE(result$tables$error)) {
    names_available <- golden_table_names(x)
    result$tables_each <- stats::setNames(lapply(names_available, function(what) {
      probe(get_results(x, what))
    }), names_available)
  }
  if (inherits(x, .golden_vcov_table_classes)) {
    result$tables_robust <- probe(get_results(x, "all", vcov_type = "robust"))
  }
  if (has("logLik")) {
    result$logLik <- probe({
      value <- stats::logLik(x)
      list(log_likelihood = as.numeric(value), df = attr(value, "df"),
           nobs = attr(value, "nobs"))
    })
  }
  if (has("nobs")) result$nobs <- probe(stats::nobs(x))
  if (has("coef")) result$coef <- probe(stats::coef(x))
  if (has("vcov")) result$vcov <- probe(stats::vcov(x))
  if (has("confint")) result$confint <- probe(stats::confint(x))
  if (has("parameter_inference")) {
    result$inference_observed <- probe(
      parameter_inference(x, vcov_type = "observed"))
    result$inference_robust <- probe(
      parameter_inference(x, vcov_type = "robust"))
  }
  if (has("simulate")) result$simulate <- probe(stats::simulate(x, seed = 1))
  result
}

#' Run the probes of any other object a verb returns
#' @param x Object (enumeration, bootstrap test, diagnostics, table, ...).
#' @param seed Seed set before each probe.
#' @return A named list of probe results.
golden_object_probes <- function(x, seed) {
  if (inherits(x, .golden_fit_classes)) return(golden_fit_probes(x, seed))
  probe <- function(expr) {
    set.seed(seed)
    golden_probe(expr)
  }
  result <- list(class = class(x))
  if (is.data.frame(x) || is.atomic(x)) {
    result$value <- probe(x)
    result$print <- probe(utils::capture.output(print(x)))
    return(result)
  }
  result$print <- probe(utils::capture.output(print(x)))
  if (golden_has_method("summary", x)) {
    result$summary <- probe(utils::capture.output(print(summary(x))))
  }
  if (golden_has_method("get_results", x)) {
    whats <- golden_whats(x)
    result$tables <- stats::setNames(
      lapply(whats, function(what) probe(get_results(x, what))), whats)
  }
  # Fits carried inside a result (an enumeration's grid, a pooled object's
  # imputations) are reduced to their likelihoods: the tables above already
  # carry what the user reads, and the full fits would repeat them.
  fits <- if (is.list(x)) x[["fits", exact = TRUE]] else NULL
  if (is.list(fits) && !is.data.frame(fits)) {
    result$fit_likelihoods <- probe(vapply(fits, function(fit) {
      value <- if (is.list(fit)) fit[["log_likelihood", exact = TRUE]] else NULL
      if (is.numeric(value)) value[1L] else NA_real_
    }, numeric(1L)))
  }
  result
}

#' The extra service probes a case may name in `services`
#'
#' Each takes the fitted object and the case and returns whatever the service
#' returns; the result is then fingerprinted by `golden_object_probes()`.
#' Cases supply `outcome`, `covariates`, `null_fit` (a function building the
#' nested model for `bootstrap_lrt`), `other_fits` (functions building more fits
#' for `compare_models`) and `newdata` where a service needs them.
#' @return A named list of functions `function(fit, case)`.
golden_services <- function() {
  list(
    diagnostics = function(fit, case) diagnostics(fit, plots = FALSE),
    descriptives = function(fit, case) descriptives(fit),
    report = function(fit, case) {
      utils::capture.output(report(fit, plots = FALSE))
    },
    predict = function(fit, case) {
      stats::predict(fit, newdata = case$newdata %||% case$data)
    },
    starting_values = function(fit, case) starting_values(fit),
    starting_values_measurement = function(fit, case) {
      starting_values(fit, what = "measurement")
    },
    sensitivity = function(fit, case) sensitivity(fit, seeds = 1:2),
    three_step_bch = function(fit, case) {
      three_step(fit, data = case$data, outcome = case$outcome, method = "bch")
    },
    three_step_proportional = function(fit, case) {
      three_step(fit, data = case$data, outcome = case$outcome,
                 method = "proportional")
    },
    three_step_modal = function(fit, case) {
      three_step(fit, data = case$data, outcome = case$outcome, method = "modal")
    },
    three_step_pairs = function(fit, case) {
      three_step(fit, data = case$data, outcome = case$outcome, method = "bch",
                 contrast = "pairs")
    },
    three_step_groups = function(fit, case) {
      three_step(fit, data = case$data, outcome = case$group_outcome,
                 level = "groups")
    },
    r3step = function(fit, case) {
      r3step(fit, data = case$data, covariates = case$covariates)
    },
    r3step_groups = function(fit, case) {
      r3step(fit, data = case$data, covariates = case$group_covariates,
             level = "groups")
    },
    inference_opg = function(fit, case) {
      parameter_inference(fit, vcov_type = "opg")
    },
    inference_fix = function(fit, case) {
      parameter_inference(fit, boundary = "fix")
    },
    inference_bootstrap = function(fit, case) {
      parameter_inference(fit, method = "bootstrap", iter = 5L, n_starts = 2L,
                          seed = 1L)
    },
    bootstrap_lrt = function(fit, case) {
      bootstrap_lrt(case$null_fit(), fit, iter = 5L, n_starts = 2L, seed = 1L)
    },
    # Named, so the comparison labels its rows by name rather than by a
    # deparse of each fit object, which would change with any internal field.
    compare_models = function(fit, case) {
      others <- lapply(case$other_fits, function(make) make())
      do.call(compare_models, c(list(fit = fit), stats::setNames(
        others, paste0("other_", seq_along(others)))))
    },
    simulate = function(fit, case) stats::simulate(fit, seed = 1L),
    get_tna = function(fit, case) {
      if (!requireNamespace("tna", quietly = TRUE)) return("<tna not installed>")
      get_tna(fit)
    },
    get_group_tna = function(fit, case) {
      if (!requireNamespace("tna", quietly = TRUE)) return("<tna not installed>")
      get_group_tna(fit)
    }
  )
}

#' Fingerprint a fitted object (or a refusal) for one catalogue case
#'
#' @param object The value of `case$fit()`, as returned by `golden_probe(...,
#'   clean = FALSE)`: a list with `value`, `error` and `conditions`.
#' @param case The catalogue case.
#' @param seed The case's seed.
#' @return A plain list: `case`, `engine`, `fit_conditions`, and either
#'   `refusal` (the fit was refused) or `probes` and `services`.
golden_fingerprint <- function(object, case, seed = golden_case_seed(case$name)) {
  result <- list(
    case = case$name,
    engine = case$engine,
    fit_conditions = object$conditions
  )
  if (isTRUE(object$error)) {
    result$refusal <- object$value
    return(result)
  }
  fit <- object$value
  result$probes <- golden_object_probes(fit, seed)
  registry <- golden_services()
  unknown <- setdiff(case$services, names(registry))
  if (length(unknown) > 0L) {
    stop(sprintf("Case `%s` names unknown services: %s", case$name,
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  result$services <- stats::setNames(lapply(case$services, function(name) {
    set.seed(seed)
    service <- golden_probe(registry[[name]](fit, case), clean = FALSE)
    if (isTRUE(service$error)) return(service)
    list(value = golden_object_probes(service$value, seed), error = FALSE,
         conditions = service$conditions)
  }), case$services)
  result
}

#' A fixed seed per case, derived from its name
#' @param name Case name.
#' @return A single integer.
golden_case_seed <- function(name) {
  codes <- utf8ToInt(name)
  as.integer(sum(codes * seq_along(codes)) %% 99991L + 1L)
}

#' Every unexpected (non-refusal) error recorded anywhere in a fingerprint
#' @param fingerprint A fingerprint.
#' @return A data.frame with `path`, `error_class`, `message`.
golden_unexpected_errors <- function(fingerprint) {
  walk <- function(x, path) {
    if (!is.list(x)) return(list())
    if (!is.null(x[["error_class", exact = TRUE]]) &&
        identical(x[["expected", exact = TRUE]], FALSE)) {
      return(list(data.frame(path = path,
                             error_class = paste(x$error_class, collapse = ","),
                             message = x$message, stringsAsFactors = FALSE)))
    }
    children <- names(x) %||% rep("", length(x))
    unlist(Map(function(child, name) {
      walk(child, paste(path, if (nzchar(name)) name else "?", sep = "/"))
    }, x, children), recursive = FALSE)
  }
  found <- walk(fingerprint, fingerprint$case)
  if (length(found) == 0L) {
    return(data.frame(path = character(), error_class = character(),
                      message = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, found)
}

# ---------------------------------------------------------------------------
# Comparing two fingerprints.

# Tolerances: a path naming a log likelihood is held to 1e-10 relative,
# everything else to 1e-8; an absolute difference of at most 1e-10 always
# passes (values at or near zero, where a relative difference is meaningless).
.golden_tolerance <- list(log_likelihood = 1e-10, default = 1e-8,
                          absolute = 1e-10)

#' Is `path` a log likelihood (held to the tighter tolerance)?
#' @param path A fingerprint path.
#' @return A single logical.
golden_is_loglik_path <- function(path) {
  grepl("log_?lik|loglik", path, ignore.case = TRUE)
}

#' One comparison outcome
#' @param status 0 identical, 1 within tolerance, 2 different.
#' @param max_rel Largest relative difference seen.
#' @param path First differing path (or the path of `max_rel`).
#' @param detail A short description of the first difference.
#' @return A list.
golden_outcome <- function(status = 0L, max_rel = 0, path = NA_character_,
                           detail = NA_character_) {
  list(status = status, max_rel = max_rel, path = path, detail = detail)
}

#' Combine several outcomes: the worst status, the largest difference, and the
#' first path that carries the worst status
#' @param outcomes A list of outcomes.
#' @return One outcome.
golden_combine <- function(outcomes) {
  if (length(outcomes) == 0L) return(golden_outcome())
  status <- vapply(outcomes, function(o) o$status, integer(1L))
  max_rel <- vapply(outcomes, function(o) o$max_rel, numeric(1L))
  worst <- max(status)
  if (worst == 0L) return(golden_outcome())
  pick <- if (worst == 2L) which(status == 2L)[1L] else which.max(max_rel)
  golden_outcome(worst, max(max_rel), outcomes[[pick]]$path,
                 outcomes[[pick]]$detail)
}

#' Short text for a value, for the detail column
#' @param x A value.
#' @return A single string.
golden_brief <- function(x) {
  text <- paste(utils::capture.output(utils::str(x, give.attr = FALSE,
                                                 vec.len = 2L)),
                collapse = " ")
  substr(gsub("\\s+", " ", text), 1L, 70L)
}

#' Compare two doubles vectors (attributes already handled)
#' @param a,b Numeric vectors without attributes.
#' @param path Their path.
#' @return An outcome.
golden_compare_numbers <- function(a, b, path) {
  if (length(a) != length(b)) {
    return(golden_outcome(2L, Inf, path, sprintf("length %d vs %d", length(a),
                                                 length(b))))
  }
  if (!identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b))) {
    return(golden_outcome(2L, Inf, path, "NA/NaN pattern differs"))
  }
  keep <- !is.na(a)
  a <- a[keep]
  b <- b[keep]
  infinite <- is.infinite(a) | is.infinite(b)
  if (!identical(a[infinite], b[infinite])) {
    return(golden_outcome(2L, Inf, path, "infinite values differ"))
  }
  a <- a[!infinite]
  b <- b[!infinite]
  difference <- abs(a - b)
  if (length(difference) == 0L || all(difference == 0)) return(golden_outcome())
  relative <- ifelse(difference <= .golden_tolerance$absolute, 0,
                     difference / pmax(abs(a), abs(b)))
  worst <- which.max(relative)
  max_rel <- relative[worst]
  tolerance <- if (golden_is_loglik_path(path)) {
    .golden_tolerance$log_likelihood
  } else {
    .golden_tolerance$default
  }
  detail <- sprintf("baseline %.17g, now %.17g", a[worst], b[worst])
  golden_outcome(if (max_rel <= tolerance) 1L else 2L, max_rel, path, detail)
}

#' Compare a baseline fingerprint (or any part of one) with a fresh one
#'
#' Characters, logicals, integers, names and every attribute must be
#' identical; doubles are compared with `.golden_tolerance`.
#' @param baseline,current The two values.
#' @param path Path of this value, for the report.
#' @return An outcome: `status` (0 identical, 1 within tolerance, 2
#'   different), `max_rel`, `path` and `detail`.
golden_compare <- function(baseline, current, path = "") {
  if (identical(baseline, current)) return(golden_outcome())
  if (!identical(typeof(baseline), typeof(current))) {
    return(golden_outcome(2L, Inf, path, sprintf("type %s vs %s",
                                                 typeof(baseline),
                                                 typeof(current))))
  }
  attributes_a <- attributes(baseline)
  attributes_b <- attributes(current)
  names_a <- sort(names(attributes_a))
  names_b <- sort(names(attributes_b))
  if (!identical(names_a, names_b)) {
    return(golden_outcome(2L, Inf, path, sprintf(
      "attributes [%s] vs [%s]", paste(names_a, collapse = ","),
      paste(names_b, collapse = ","))))
  }
  attribute_outcomes <- lapply(names_a, function(name) {
    golden_compare(attributes_a[[name]], attributes_b[[name]],
                   sprintf("%s@%s", path, name))
  })
  bare_a <- baseline
  bare_b <- current
  attributes(bare_a) <- NULL
  attributes(bare_b) <- NULL
  value_outcomes <- if (is.list(bare_a)) {
    if (length(bare_a) != length(bare_b)) {
      list(golden_outcome(2L, Inf, path, sprintf("length %d vs %d",
                                                 length(bare_a),
                                                 length(bare_b))))
    } else {
      labels <- names(baseline) %||% rep("", length(bare_a))
      labels <- ifelse(nzchar(labels), labels,
                       sprintf("[[%d]]", seq_along(bare_a)))
      Map(function(a, b, label) golden_compare(a, b, paste0(path, "/", label)),
          bare_a, bare_b, labels)
    }
  } else if (is.double(bare_a)) {
    list(golden_compare_numbers(bare_a, bare_b, path))
  } else if (is.complex(bare_a)) {
    list(golden_compare_numbers(Re(bare_a), Re(bare_b), path),
         golden_compare_numbers(Im(bare_a), Im(bare_b), path))
  } else if (identical(bare_a, bare_b)) {
    list(golden_outcome())
  } else {
    first <- if (length(bare_a) == length(bare_b)) {
      which(!(bare_a == bare_b) | xor(is.na(bare_a), is.na(bare_b)))[1L]
    } else {
      NA_integer_
    }
    detail <- if (is.na(first)) {
      sprintf("%s vs %s", golden_brief(bare_a), golden_brief(bare_b))
    } else {
      sprintf("[%d] %s vs %s", first, golden_brief(bare_a[first]),
              golden_brief(bare_b[first]))
    }
    list(golden_outcome(2L, Inf, path, detail))
  }
  golden_combine(c(attribute_outcomes, value_outcomes))
}
