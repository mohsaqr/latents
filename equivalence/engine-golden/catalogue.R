# The golden catalogue: one case per engine and documented option.
#
# A case is
#   list(name, engine, covers, fit = function() <the fit>, data = <data the fit
#        used, for verbs that take data>, services = <extra probes>, ...)
# where `...` carries whatever a service needs (`outcome`, `covariates`,
# `group_outcome`, `group_covariates`, `null_fit`, `other_fits`, `newdata`).
#
# `fit()` must be self-contained and seeded (`seed =` wherever the verb takes
# one), small and fast. A fit the package refuses with a classed error is a
# case too: the refusal is behaviour to preserve.
#
# The cases live in section files (catalogue-<section>.R) next to this one,
# each defining `golden_cases_<section>()`; `golden_catalogue()` joins them.

#' Build one catalogue case
#' @param name Unique case name (also the RDS file name).
#' @param engine Which estimation engine the case exercises.
#' @param covers One line on what the case covers, for the README.
#' @param fit A function of no arguments returning the fitted object.
#' @param data The data the fit used, or `NULL`.
#' @param services Names of extra service probes (see `golden_services()`).
#' @param ... Inputs those services need.
#' @return A case list.
golden_case <- function(name, engine, covers, fit, data = NULL,
                        services = character(), ...) {
  stopifnot(
    "`name` must be a single string" = is.character(name) && length(name) == 1L,
    "`name` must be usable as a file name" = grepl("^[A-Za-z0-9_.-]+$", name),
    "`fit` must be a function" = is.function(fit),
    "`services` must be a character vector" = is.character(services)
  )
  c(list(name = name, engine = engine, covers = covers, fit = fit, data = data,
         services = services), list(...))
}

# ---------------------------------------------------------------------------
# Small data shared by the sections. Every fit uses a subset small enough to
# fit in about a second.

#' The first `n` students of `course_engagement` (about 14 rows each)
golden_engagement <- function(n = 20L) {
  data <- subset(course_engagement, student <= n)
  row.names(data) <- NULL
  data
}

#' The first `n` students of `student_esm`
golden_esm <- function(n = 30L) {
  data <- subset(student_esm, student <= n)
  row.names(data) <- NULL
  data
}

#' The first `n` students of `growth_scores`
golden_growth <- function(n = 60L) {
  data <- subset(growth_scores, student <= n)
  row.names(data) <- NULL
  data
}

#' The first `n` students of `study_hours`
golden_hours <- function(n = 80L) {
  data <- subset(study_hours, student <= n)
  row.names(data) <- NULL
  data
}

#' Continuous `course_engagement` indicators
golden_engagement_vars <- c("browse", "lectures", "forum_read")

#' Every catalogue case, in section order
#' @param pattern Optional regular expression; only matching case names.
#' @return A named list of cases.
golden_catalogue <- function(pattern = NULL) {
  sections <- c("profiles", "covariates_families", "transitions",
                "regression", "services")
  cases <- unlist(lapply(sections, function(section) {
    builder <- get0(paste0("golden_cases_", section), mode = "function")
    if (is.null(builder)) {
      stop(sprintf("Section `%s` is not loaded.", section), call. = FALSE)
    }
    builder()
  }), recursive = FALSE)
  names(cases) <- vapply(cases, function(case) case$name, character(1L))
  duplicated_names <- unique(names(cases)[duplicated(names(cases))])
  if (length(duplicated_names) > 0L) {
    stop(sprintf("Duplicate case names: %s",
                 paste(duplicated_names, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(pattern)) cases <- cases[grepl(pattern, names(cases))]
  cases
}
