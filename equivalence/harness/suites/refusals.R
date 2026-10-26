# Refusal parity: one degenerate input, every entry point, one answer.
#
# The 2026-10-26 review found the same defect in the data accepted by one verb
# and refused by another: exactly constant indicators passed where their
# computed standard deviation was 1e-17 rather than zero (R1), and
# `enumerate_lpa()` returned a grid of failed cells, with no condition, for a
# column that did not exist while `multilpa()` refused it. A refusal that one
# entry point raises and another omits is invisible to every test that uses
# only the first.
#
# Each row puts one degenerate input through one verb and passes when the verb
# refuses it with the expected condition class (or, for input that is legal
# but suspicious, fits and warns with the expected class). An unclassed error,
# a different class, or a silent fit fails.

suite_refusals <- function() {
  cases <- c(.refusal_continuous_cases(), .refusal_categorical_cases())
  rows <- do.call(rbind, lapply(cases, .refusal_row))
  rows$source <- "refusal parity"
  rows
}

#' What a call did: the condition class it raised, or that it fitted
#'
#' Warnings are recorded, not muffled silently: the classes of every warning
#' are returned so a case can require one.
#'
#' @param expression An unevaluated call.
#' @param envir Where to evaluate it.
#' @return A list with `outcome` (`"fit"`, the first condition class of the
#'   package, or `"unclassed"`) and `warnings` (their classes).
.refusal_outcome <- function(expression, envir) {
  warnings <- character()
  outcome <- withCallingHandlers(
    tryCatch({
      eval(expression, envir)
      "fit"
    }, error = function(condition) {
      own <- grep("^(latents|multilpa)_", class(condition), value = TRUE)
      if (length(own) == 0L) "unclassed" else own[[1L]]
    }),
    warning = function(condition) {
      warnings <<- c(warnings, class(condition)[[1L]])
      invokeRestart("muffleWarning")
    })
  list(outcome = outcome, warnings = unique(warnings))
}

#' One comparison row for one case
#' @param case A list with `defect`, `verb`, `call`, `data` and `expect`
#'   (a condition class, or `"fit"` with `warning` naming the class required).
#' @return A comparison row.
.refusal_row <- function(case) {
  envir <- new.env(parent = asNamespace("latents"))
  envir$data <- case$data
  seen <- .refusal_outcome(case$call, envir)
  passed <- if (identical(case$expect, "fit")) {
    identical(seen$outcome, "fit") && case$warning %in% seen$warnings
  } else identical(seen$outcome, case$expect)
  expected <- if (identical(case$expect, "fit")) sprintf("fits, warning %s", case$warning) else
    sprintf("refuses with %s", case$expect)
  got <- if (identical(seen$outcome, "fit")) {
    sprintf("fitted%s", if (length(seen$warnings)) sprintf(", warning %s",
                                                           paste(seen$warnings, collapse = "/")) else "")
  } else seen$outcome
  compare_values(sprintf("%s: %s %s (got: %s)", case$defect, case$verb, expected, got),
                 reference = 1, obtained = as.numeric(passed), tolerance = 0.5)
}

#' Data for the refusal cases: 30 groups of 4 occasions, two Gaussian and two
#' categorical indicators and a predictor
.refusal_data <- function() {
  set.seed(20261026L)
  n <- 120L
  data.frame(g = rep(seq_len(30L), each = 4L), t = rep(seq_len(4L), 30L),
             y1 = stats::rnorm(n), y2 = stats::rnorm(n),
             c1 = factor(sample(c("a", "b", "c"), n, TRUE)),
             c2 = factor(sample(c("p", "q"), n, TRUE)), x = stats::rnorm(n))
}

#' The Gaussian entry points, each as a call on `data` with indicators `y1`
#' and `y2` (or the outcome `y1` for the regression)
.refusal_gaussian_verbs <- function(vars = c("y1", "y2")) {
  list(
    multilpa = bquote(multilpa(data, vars = .(vars), id = "g", n_profiles = 2,
                               n_group_classes = 1, n_starts = 1, seed = 1)),
    lpa = bquote(lpa(data, .(vars), n_profiles = 2, n_starts = 1, seed = 1)),
    lta = bquote(lta(data, .(vars), "g", n_profiles = 2, time = "t", n_starts = 1,
                     seed = 1)),
    enumerate_lpa = bquote(enumerate_lpa(data, .(vars), n_profiles = 1:2,
                                         n_starts = 1, seed = 1)),
    enumerate_classes = bquote(enumerate_classes(data, .(vars), id = "g",
                                                 n_profiles = 1:2, n_starts = 1,
                                                 seed = 1)),
    mixture_regression = bquote(mixture_regression(.(as.name(vars[[1L]])) ~ x, data,
                                                   n_classes = 2, n_starts = 1,
                                                   seed = 1)))
}

#' Gaussian defects, put through every Gaussian entry point they apply to
.refusal_continuous_cases <- function() {
  base <- .refusal_data()
  defects <- list(
    list(defect = "an exactly constant indicator (0.1 in every row)",
         data = within(base, y1 <- 0.1), expect = "latents_bad_data"),
    list(defect = "an indicator missing in every row",
         data = within(base, y1 <- NA_real_), expect = "latents_bad_data",
         # The regression refuses any missing outcome under its default
         # `missing = "error"`, with the class that contract documents.
         override = c(mixture_regression = "latents_missing_data")),
    list(defect = "an infinite value", data = within(base, y1[5L] <- Inf),
         expect = "latents_bad_data"),
    list(defect = "a text indicator", data = within(base, y1 <- as.character(y1)),
         expect = "latents_bad_data"),
    list(defect = "no rows", data = base[0L, ], expect = "latents_bad_data"),
    list(defect = "a missing group identifier", data = within(base, g[3L] <- NA),
         expect = "latents_bad_data",
         # Single-level verbs do not read the group column.
         skip = c("lpa", "enumerate_lpa", "mixture_regression")))
  cases <- unlist(lapply(defects, function(defect) {
    verbs <- .refusal_gaussian_verbs()
    verbs <- verbs[setdiff(names(verbs), defect$skip)]
    lapply(names(verbs), function(verb) {
      list(defect = defect$defect, verb = verb, call = verbs[[verb]], data = defect$data,
           expect = if (verb %in% names(defect$override)) defect$override[[verb]] else
             defect$expect)
    })
  }), recursive = FALSE)
  # A misspelled indicator, through every verb that takes `vars`.
  misspelled <- .refusal_gaussian_verbs(c("y1", "nope"))
  misspelled$mixture_regression <- NULL
  c(cases, lapply(names(misspelled), function(verb) {
    list(defect = "a misspelled indicator", verb = verb, call = misspelled[[verb]],
         data = base, expect = "latents_bad_data")
  }))
}

#' The categorical entry points
.refusal_categorical_verbs <- function(vars = c("c1", "c2")) {
  list(
    multilpa = bquote(multilpa(data, vars = .(vars), id = "g", n_profiles = 2,
                               n_group_classes = 1, categorical = .(vars),
                               n_starts = 1, seed = 1)),
    lca = bquote(lca(data, .(vars), n_classes = 2, n_starts = 1, seed = 1)),
    lta = bquote(lta(data, .(vars), "g", n_profiles = 2, time = "t",
                     categorical = .(vars), n_starts = 1, seed = 1)),
    enumerate_lca = bquote(enumerate_lca(data, .(vars), n_classes = 1:2, n_starts = 1,
                                         seed = 1)))
}

#' Categorical defects, put through every categorical entry point
.refusal_categorical_cases <- function() {
  base <- .refusal_data()
  defects <- list(
    list(defect = "a categorical indicator with one category",
         data = within(base, c1 <- factor(rep("a", nrow(base)))),
         expect = "latents_bad_categorical"),
    list(defect = "a categorical indicator missing in every row",
         data = within(base, c1 <- factor(rep(NA, nrow(base)), levels = c("a", "b"))),
         expect = "latents_bad_categorical"),
    # Legal but almost certainly a mistake: fitted, with a classed warning.
    list(defect = "a continuous variable declared categorical",
         data = within(base, c1 <- stats::rnorm(nrow(base))), expect = "fit",
         warning = "latents_suspect_categorical"))
  unlist(lapply(defects, function(defect) {
    verbs <- .refusal_categorical_verbs()
    lapply(names(verbs), function(verb) {
      list(defect = defect$defect, verb = verb, call = verbs[[verb]], data = defect$data,
           expect = defect$expect, warning = defect$warning)
    })
  }), recursive = FALSE)
}
