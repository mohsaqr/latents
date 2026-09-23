# Equivalence harness for multilpa.
#
# Every external comparison in this package reports through one tidy frame, so
# that the evidence can be regenerated and read as a single table rather than
# recovered by reading a pile of scripts. A suite is a file in suites/ that
# defines `suite_<name>()` returning the frame built by `compare_values()`.
#
# Run every suite from the project root:
#   Rscript equivalence/harness/run.R
#
# Three things the harness guarantees, each added after the 2026-09-20 review
# found the saved report had been attributed to a version it was never run
# against:
#
#   1. `check_validation_api()` parses every validation script and checks each
#      call to a multilpa verb against that verb's current formals, before any
#      suite is loaded. A script left on a renamed argument is named with its
#      file, line and the replacement argument, rather than dying at suite load
#      with a bare "unused arguments".
#   2. A suite whose external package is not installed raises
#      `multilpa_suite_skipped` through `require_suite_packages()` and is
#      recorded as skipped with its reason. It never silently disappears, and
#      it never aborts the rest of the run.
#   3. The returned object carries the multilpa version, the R version and
#      every external package version the run actually used, and `run.R` writes
#      them beside the numbers. A report can no longer be read as evidence for
#      a version it was not produced under.

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

#' The tolerances this harness compares at
#'
#' Pinned here rather than left implicit at each call site, so that loosening
#' one is a visible edit. `compare_values()` has no default tolerance: every
#' comparison states the tolerance it claims.
#'
#' @format A named numeric vector.
#' - `machine`: `sqrt(.Machine$double.eps)`, the default for an agreement that
#'   nothing but the two implementations bounds.
#' - `likelihood`: 1e-5, for a maximised log-likelihood two different
#'   multi-start schemes reached independently.
#' - `parameter`: 1e-4, for a probability or a mean read off two independent
#'   optimisers stopped at their own convergence tolerances.
#' - `count`: 0.5, for an integer-valued quantity such as a free-parameter
#'   count or a sample size, where anything inside half a unit is the same
#'   integer.
.equivalence_tolerances <- c(machine = sqrt(.Machine$double.eps),
                             likelihood = 1e-5, parameter = 1e-4, count = 0.5)

#' Compare reference values against values this package produced
#'
#' @param quantity Character vector naming each compared quantity.
#' @param reference Numeric vector of reference values from the external source.
#' @param obtained Numeric vector of values multilpa produced.
#' @param tolerance Numeric tolerance, recycled across quantities. Required:
#'   see [.equivalence_tolerances] for the pinned values this harness uses.
#' @param scale `"absolute"` compares `|reference - obtained|`; `"relative"`
#'   divides that difference by `max(|reference|, 1)`, for quantities whose
#'   magnitude makes an absolute tolerance meaningless.
#' @return A `data.frame` with one row per quantity and columns `quantity`,
#'   `reference`, `obtained`, `difference`, `tolerance`, `agrees`.
compare_values <- function(quantity, reference, obtained, tolerance,
                           scale = c("absolute", "relative")) {
  scale <- match.arg(scale)
  stopifnot(
    "`quantity` must be character" = is.character(quantity),
    "`reference` must be numeric" = is.numeric(reference),
    "`obtained` must be numeric" = is.numeric(obtained),
    "`reference` and `obtained` must have one value per quantity" =
      length(reference) == length(quantity) && length(obtained) == length(quantity),
    "`tolerance` must be stated explicitly and be positive" =
      !missing(tolerance) && is.numeric(tolerance) && all(tolerance > 0)
  )
  divisor <- if (identical(scale, "relative")) pmax(abs(reference), 1) else 1
  difference <- abs(reference - obtained) / divisor
  data.frame(quantity = quantity, reference = reference, obtained = obtained,
             difference = difference, tolerance = tolerance,
             agrees = difference <= tolerance, stringsAsFactors = FALSE)
}

#' Compare against a value read from published digits
#'
#' A value printed to `digits` decimals states its quantity only to within half
#' a unit in the last place. Reporting the raw gap against such a reference
#' measures how the source was typeset, not how the implementation behaves: a
#' G-squared printed as `15.3` sits up to 0.05 from whatever was computed, so an
#' exact reproduction still shows a difference of that order.
#'
#' This expresses the gap as a fraction of the rounding the printed value
#' admits. A value of 0.77 means the computed quantity is well inside the
#' printed precision and therefore reproduces the published digits; anything
#' above 1 means it does not, and no amount of rounding explains it.
#'
#' @param quantity Character vector naming each compared quantity.
#' @param reference Numeric vector of published values.
#' @param obtained Numeric vector of values multilpa produced.
#' @param digits Number of decimals each reference was printed to, recycled
#'   across quantities.
#' @return A `data.frame` shaped like [compare_values()], whose `difference` is
#'   in units of the half-unit rounding bound and whose `tolerance` is 1.
compare_printed <- function(quantity, reference, obtained, digits) {
  stopifnot(
    "`digits` must be a nonnegative whole number" =
      is.numeric(digits) && all(digits >= 0) && all(digits == floor(digits)),
    "`reference` and `obtained` must have one value per quantity" =
      length(reference) == length(quantity) && length(obtained) == length(quantity))
  bound <- 0.5 * 10^(-digits)
  data.frame(quantity = quantity, reference = reference, obtained = obtained,
             difference = abs(reference - obtained) / bound,
             tolerance = 1, agrees = abs(reference - obtained) <= bound,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Guarding external dependencies
# ---------------------------------------------------------------------------

#' Require the external packages a suite compares against
#'
#' Raises `multilpa_suite_skipped` rather than an error when a package is
#' missing, so [run_equivalence()] records the suite as skipped with its reason
#' instead of aborting the run or dropping the suite silently.
#'
#' @param ... Package names the calling suite needs.
#' @param reason Optional sentence explaining what the packages are needed for.
#' @return Invisibly `TRUE` when every package is installed; otherwise a
#'   `multilpa_suite_skipped` condition is signalled.
require_suite_packages <- function(..., reason = NULL) {
  packages <- c(...)
  stopifnot("`...` must name packages" = is.character(packages) && length(packages) > 0L)
  present <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (all(present)) return(invisible(TRUE))
  missing_packages <- paste(packages[!present], collapse = ", ")
  message_text <- sprintf("not installed: %s%s", missing_packages,
                          if (is.null(reason)) "" else sprintf(" (%s)", reason))
  stop(errorCondition(message_text, class = "multilpa_suite_skipped",
                      packages = packages[!present], call = NULL))
}

#' Require the retained program output a suite reads its reference from
#'
#' Same contract as [require_suite_packages()], for a suite whose reference is
#' a file rather than a package: a missing artifact is a skip with a reason,
#' not a cryptic connection error halfway through a run.
#'
#' @param ... Paths the calling suite needs.
#' @param reason Optional sentence explaining what the files are.
#' @return Invisibly `TRUE` when every file exists; otherwise a
#'   `multilpa_suite_skipped` condition is signalled.
require_suite_files <- function(..., reason = NULL) {
  paths <- c(...)
  stopifnot("`...` must name files" = is.character(paths) && length(paths) > 0L)
  absent <- paths[!file.exists(paths)]
  if (length(absent) == 0L) return(invisible(TRUE))
  stop(errorCondition(
    sprintf("missing artifact(s): %s%s", paste(absent, collapse = ", "),
            if (is.null(reason)) "" else sprintf(" (%s)", reason)),
    class = "multilpa_suite_skipped", files = absent, call = NULL))
}

# ---------------------------------------------------------------------------
# Stale-API detection
# ---------------------------------------------------------------------------

#' The 0.8.0 argument renames, as a lookup from the old spelling
#'
#' multilpa 0.8.0 renamed public arguments to match the sibling package
#' Nestimate and shipped no deprecation shim, so a script left on the old
#' spelling fails with `unused argument`. The review of 2026-09-20 found the
#' whole validation tree still on the old names. This table turns that failure
#' into a named replacement.
#'
#' @format A named character vector: names are the retired spellings, values
#'   the current ones.
.multilpa_api_renames <- c(object = "x", indicators = "vars", group = "id",
                           level_ci = "ci_level", p_adjust = "adjust",
                           n_boot = "iter", profiles = "n_profiles",
                           group_classes = "n_group_classes")

#' Verbs that left the exported API, as a lookup to what replaced them
#'
#' The argument check below only inspects calls to verbs the package still
#' exports, so a call to a removed verb was never looked at. `get_data()` became
#' `get_results()` in 0.4.0 and seven validation scripts kept calling the old
#' name unnoticed. A removed verb is listed here with its replacement, or with
#' where it went if it was deferred rather than renamed.
#'
#' @format A named character vector: names are the removed verbs, values the
#'   replacement to report.
.multilpa_removed_verbs <- c(
  get_data = "get_results",
  fit_transitions = "lta",
  fit_random_intercept = "deferred: deferred_verb(\"random-intercept.R\", \"fit_random_intercept\")",
  lmr_lrt = "deferred: deferred_verb(\"lmr.R\", \"lmr_lrt\")")

#' Calls to removed verbs in one file
#'
#' A file that assigns the name itself -- loading a deferred verb with
#' `deferred_verb()`, say -- is calling its own binding, not the removed export,
#' and is not reported.
#'
#' @param path Path to an R script.
#' @return A `data.frame` with columns `file`, `line`, `call`, `argument` and
#'   `replacement`, one row per offending call; `argument` is `NA`.
.validation_removed_calls <- function(path) {
  parsed <- parse(path, keep.source = TRUE)
  data <- utils::getParseData(parsed)
  empty <- data.frame(file = character(), line = integer(), call = character(),
                      argument = character(), replacement = character(),
                      stringsAsFactors = FALSE)
  if (is.null(data) || nrow(data) == 0L) return(empty)
  assigned <- vapply(as.list(parsed), function(expression) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.name(expression[[2L]])) as.character(expression[[2L]]) else NA_character_
  }, character(1))
  removed <- setdiff(names(.multilpa_removed_verbs), assigned)
  called <- data[data$token == "SYMBOL_FUNCTION_CALL" & data$text %in% removed, ,
                 drop = FALSE]
  if (nrow(called) == 0L) return(empty)
  data.frame(file = path, line = called$line1, call = called$text,
             argument = NA_character_,
             replacement = unname(.multilpa_removed_verbs[called$text]),
             stringsAsFactors = FALSE)
}

#' Every named argument passed to a multilpa verb in one file
#'
#' Parses the file and reads argument names off the parse data, so the line
#' number reported is the line the argument is written on rather than the line
#' the enclosing expression starts on.
#'
#' @param path Path to an R script.
#' @param verbs Character vector of function names to look for.
#' @return A `data.frame` with columns `file`, `line`, `call` and `argument`,
#'   one row per named argument of a matching call.
.validation_call_arguments <- function(path, verbs) {
  parsed <- parse(path, keep.source = TRUE)
  data <- utils::getParseData(parsed)
  if (is.null(data) || nrow(data) == 0L) {
    return(data.frame(file = character(), line = integer(), call = character(),
                      argument = character(), stringsAsFactors = FALSE))
  }
  called <- data[data$token == "SYMBOL_FUNCTION_CALL" & data$text %in% verbs, , drop = FALSE]
  if (nrow(called) == 0L) {
    return(data.frame(file = character(), line = integer(), call = character(),
                      argument = character(), stringsAsFactors = FALSE))
  }
  # The function symbol sits inside its own `expr`, whose parent is the call.
  call_ids <- data$parent[match(called$parent, data$id)]
  rows <- lapply(seq_len(nrow(called)), function(index) {
    arguments <- data[data$token == "SYMBOL_SUB" & data$parent == call_ids[[index]], , drop = FALSE]
    if (nrow(arguments) == 0L) return(NULL)
    data.frame(file = path, line = arguments$line1, call = called$text[[index]],
               argument = arguments$text, stringsAsFactors = FALSE)
  })
  stacked <- do.call(rbind, rows)
  if (is.null(stacked)) {
    data.frame(file = character(), line = integer(), call = character(),
               argument = character(), stringsAsFactors = FALSE)
  } else {
    stacked
  }
}

#' Check every validation script against multilpa's current signatures
#'
#' The defect this exists to catch: a public argument is renamed, the package's
#' own tests are updated, and the validation scripts are not. The harness then
#' fails at suite load with `unused argument`, which names neither the file nor
#' the replacement, and the saved report keeps claiming a pass.
#'
#' Every named argument of every call to an exported multilpa verb is checked
#' against that verb's formals. Partial matching is honoured, because R honours
#' it, and a verb with `...` only has its arguments checked against the formals
#' before the dots.
#'
#' @param root Project root directory.
#' @param paths Character vector of scripts to check; defaults to every `.R`
#'   file under `equivalence/` and `validation/`.
#' @return A `data.frame` with one row per offending call site and columns
#'   `file`, `line`, `call`, `argument` and `replacement`. Zero rows means
#'   every call matches the loaded package's signatures.
check_validation_api <- function(root = ".", paths = NULL) {
  stopifnot("multilpa must be loaded before its signatures can be checked" =
              "multilpa" %in% loadedNamespaces())
  verbs <- getNamespaceExports("multilpa")
  verbs <- verbs[vapply(verbs, function(name) {
    is.function(get(name, envir = asNamespace("multilpa")))
  }, logical(1))]
  if (is.null(paths)) {
    # Both trees call the package's verbs: the equivalence material and what
    # is left under validation/ (the recovery study and the demo). A tree
    # missing from this list would stop being checked without any error.
    paths <- list.files(file.path(root, c("equivalence", "validation")),
                        pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  }
  removed <- do.call(rbind, lapply(paths, .validation_removed_calls))
  used <- do.call(rbind, lapply(paths, .validation_call_arguments, verbs = verbs))
  empty <- data.frame(file = character(), line = integer(), call = character(),
                      argument = character(), replacement = character(),
                      stringsAsFactors = FALSE)
  if (is.null(used) || nrow(used) == 0L) return(rbind(empty, removed))
  known <- vapply(seq_len(nrow(used)), function(index) {
    formal_names <- names(formals(get(used$call[[index]], envir = asNamespace("multilpa"))))
    argument <- used$argument[[index]]
    dots <- match("...", formal_names)
    # Arguments before `...` may be abbreviated, exactly as R matches them;
    # those after it must be given in full.
    matchable <- if (is.na(dots)) formal_names else formal_names[seq_len(dots - 1L)]
    if (argument %in% formal_names) return(TRUE)
    # A retired spelling is stale even when it partially matches a current
    # argument. `group` abbreviates `group_covariates`, so an old
    # `multilpa(group = "school")` does not fail: R silently passes the ID
    # column as a group covariate instead.
    if (argument %in% names(.multilpa_api_renames)) return(FALSE)
    if (sum(startsWith(matchable, argument)) == 1L) {
      return(TRUE)
    }
    # A verb with `...` forwards what it does not recognise, so an unmatched
    # argument there is legitimate unless it is one of the spellings 0.8.0
    # retired -- which would be forwarded to a verb that no longer accepts it.
    !is.na(dots) && !(argument %in% names(.multilpa_api_renames))
  }, logical(1))
  stale <- used[!known, , drop = FALSE]
  replacement <- unname(.multilpa_api_renames[stale$argument])
  stale$replacement <- ifelse(is.na(replacement), "(no known replacement)", replacement)
  stale <- rbind(empty, stale, removed)
  row.names(stale) <- NULL
  stale
}

#' Stop with a readable account of every stale call site
#' @param stale The frame returned by [check_validation_api()].
#' @return Invisibly `TRUE` when `stale` has no rows; otherwise a
#'   `multilpa_stale_validation_api` condition is signalled.
.stop_if_stale_api <- function(stale) {
  if (nrow(stale) == 0L) return(invisible(TRUE))
  # A removed verb has no argument to name; the call itself is what is stale.
  site <- ifelse(is.na(stale$argument), sprintf("%s()", stale$call),
                 sprintf("%s(%s = )", stale$call, stale$argument))
  lines <- sprintf("  %s:%d  %s -> %s", stale$file, stale$line, site,
                   stale$replacement)
  stop(errorCondition(
    paste(c(sprintf("%d validation call site(s) do not match multilpa %s:",
                    nrow(stale), utils::packageVersion("multilpa")),
            lines), collapse = "\n"),
    class = "multilpa_stale_validation_api", stale = stale, call = NULL))
}

# ---------------------------------------------------------------------------
# Running the suites
# ---------------------------------------------------------------------------

#' Load and run every registered equivalence suite
#'
#' Checks every validation script against the loaded package's signatures
#' before any suite is sourced, so a script left on a renamed argument is
#' reported with its file and line rather than failing at suite load.
#'
#' @param suites Character vector of suite names, or `NULL` for every suite
#'   found in `equivalence/harness/suites`.
#' @param root Project root directory.
#' @param check_api Whether to check every validation script against the loaded
#'   package's signatures first. Only turn this off to inspect a known-stale
#'   tree.
#' @return An object of class `multilpa_equivalence`: a list with
#'   `comparisons` (one row per compared quantity), `suites` (one row per
#'   suite, with its status and, when skipped or failed, the reason) and
#'   `session` (the versions the run was produced under). Use
#'   `as.data.frame(x, what = )` to reach any of the three.
run_equivalence <- function(suites = NULL, root = ".", check_api = TRUE) {
  # Read the version and fingerprint the source before anything is fitted. A
  # long run can outlive the tree it started against -- the run of 2026-09-20
  # began at 0.10.0 and ended at 0.11.0, because the package was being edited
  # while it ran -- and a report that records only the closing version claims
  # evidence the loaded code never produced.
  opened <- .source_state(root)
  if (isTRUE(check_api)) .stop_if_stale_api(check_validation_api(root = root))
  suite_dir <- file.path(root, "equivalence", "harness", "suites")
  files <- sort(list.files(suite_dir, pattern = "[.]R$", full.names = TRUE))
  names(files) <- sub("[.]R$", "", basename(files))
  available <- names(files)
  if (!is.null(suites)) {
    unknown <- setdiff(suites, names(files))
    if (length(unknown) > 0L) {
      stop(errorCondition(sprintf("Unknown suite(s): %s", paste(unknown, collapse = ", ")),
                          class = "multilpa_unknown_suite", call = NULL))
    }
    files <- files[suites]
  }
  coverage <- if (is.null(suites)) sprintf("all %d suites", length(files)) else
    sprintf("PARTIAL: %d of %d suites (%s); %s not run", length(files),
            length(available), paste(names(files), collapse = ", "),
            paste(setdiff(available, names(files)), collapse = ", "))
  outcomes <- lapply(names(files), .run_one_suite, files = files)
  names(outcomes) <- names(files)
  comparisons <- do.call(rbind, lapply(outcomes, `[[`, "comparisons"))
  if (is.null(comparisons)) comparisons <- .empty_comparisons()
  row.names(comparisons) <- NULL
  status <- do.call(rbind, lapply(outcomes, `[[`, "status"))
  row.names(status) <- NULL
  structure(list(comparisons = comparisons, suites = status,
                 session = .equivalence_session(status, opened, .source_state(root),
                                                coverage)),
            class = "multilpa_equivalence")
}

#' The columns every comparison frame carries
#' @return A zero-row `data.frame` with the comparison columns.
.empty_comparisons <- function() {
  data.frame(suite = character(), source = character(), precision = character(),
             dataset = character(), quantity = character(), reference = numeric(),
             obtained = numeric(), difference = numeric(), tolerance = numeric(),
             agrees = logical(), seconds = numeric(), stringsAsFactors = FALSE)
}

#' Run one suite and record what happened to it
#'
#' A suite that raises `multilpa_suite_skipped` is recorded as skipped with its
#' reason. Any other error is recorded as a failure, with the condition's own
#' message and class kept, and the run continues so that one broken suite does
#' not hide the state of the other six. Nothing is discarded: `run.R` reports
#' every skip and every failure and exits non-zero when any suite failed.
#'
#' @param name Suite name.
#' @param files Named character vector of suite file paths.
#' @return A list with `comparisons` (possibly `NULL`) and `status`, a one-row
#'   `data.frame`.
.run_one_suite <- function(name, files) {
  message(sprintf("running suite: %s", name))
  warnings_seen <- character()
  started <- proc.time()[["elapsed"]]
  outcome <- withCallingHandlers(
    tryCatch(
      list(state = "ran", value = .evaluate_suite(name, files[[name]]), reason = NA_character_),
      multilpa_suite_skipped = function(condition) {
        list(state = "skipped", value = NULL, reason = conditionMessage(condition))
      },
      error = function(condition) {
        list(state = "failed", value = NULL,
             reason = sprintf("%s: %s", paste(class(condition), collapse = "/"),
                              conditionMessage(condition)))
      }),
    warning = function(condition) {
      warnings_seen <<- c(warnings_seen, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  elapsed <- round(proc.time()[["elapsed"]] - started, 2)
  if (identical(outcome$state, "skipped")) {
    message(sprintf("  skipped: %s", outcome$reason))
  }
  if (identical(outcome$state, "failed")) {
    message(sprintf("  FAILED: %s", outcome$reason))
  }
  if (length(warnings_seen) > 0L) {
    message(sprintf("  %d warning(s), first: %s", length(warnings_seen), warnings_seen[[1L]]))
  }
  comparisons <- if (is.null(outcome$value)) NULL else
    data.frame(outcome$value, seconds = elapsed, stringsAsFactors = FALSE)
  status <- data.frame(
    suite = name, status = outcome$state,
    compared = if (is.null(comparisons)) 0L else nrow(comparisons),
    agreed = if (is.null(comparisons)) 0L else sum(comparisons$agrees),
    warnings = length(warnings_seen), seconds = elapsed,
    reason = outcome$reason %||% NA_character_,
    first_warning = if (length(warnings_seen) == 0L) NA_character_ else warnings_seen[[1L]],
    stringsAsFactors = FALSE)
  list(comparisons = comparisons, status = status)
}

# The package Depends on R (>= 4.1.0), where base `%||%` does not exist yet.
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Source a suite file and call its entry point
#' @param name Suite name.
#' @param path Path to the suite file.
#' @return The suite's comparison frame, with `suite`, `source`, `precision`
#'   and `dataset` columns filled in.
.evaluate_suite <- function(name, path) {
  environment_for_suite <- new.env(parent = globalenv())
  sys.source(path, envir = environment_for_suite)
  entry_name <- paste0("suite_", gsub("-", "_", name))
  if (!exists(entry_name, envir = environment_for_suite, inherits = FALSE)) {
    stop(errorCondition(sprintf("%s does not define %s()", path, entry_name),
                        class = "multilpa_bad_suite", call = NULL))
  }
  result <- get(entry_name, envir = environment_for_suite)()
  stopifnot("a suite must return a data.frame with a `source` column" =
              is.data.frame(result) && "source" %in% names(result))
  if (is.null(result$dataset)) result$dataset <- NA_character_
  if (is.null(result$precision)) result$precision <- "machine"
  data.frame(suite = name, source = result$source, precision = result$precision,
             dataset = result$dataset,
             result[, c("quantity", "reference", "obtained", "difference",
                        "tolerance", "agrees")],
             stringsAsFactors = FALSE)
}

#' Record the versions a run was produced under
#'
#' The 2026-09-20 review found a saved report of 1,012 agreements being read as
#' evidence for version 0.10.0, when the scripts that produced it could not run
#' against 0.10.0 at all. A report that carries its own version cannot be
#' misread that way.
#'
#' @param status The per-suite status frame.
#' @param opened The source state recorded before the first suite ran.
#' @param closed The source state recorded after the last suite finished.
#' @param coverage A sentence saying whether every suite ran or only some.
#' @return A `data.frame` with one row per recorded component and columns
#'   `component` and `value`.
.equivalence_session <- function(status, opened, closed, coverage) {
  external <- c("mclust", "glca", "poLCA", "multilevLCA", "depmixS4", "tidySEM")
  version_of <- vapply(external, function(name) {
    if (requireNamespace(name, quietly = TRUE)) as.character(utils::packageVersion(name)) else "not installed"
  }, character(1))
  unchanged <- identical(opened$version, closed$version) &&
    identical(opened$fingerprint, closed$fingerprint)
  data.frame(
    component = c("multilpa", "multilpa at end", "R/ fingerprint",
                  "R/ fingerprint at end", "source unchanged during run",
                  "coverage", "R", "platform", "generated", "suites run",
                  "suites skipped", "suites failed", external),
    value = c(opened$version, closed$version, opened$fingerprint,
              closed$fingerprint, if (unchanged) "yes" else "NO", coverage,
              R.version.string, R.version$platform,
              format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
              sum(status$status == "ran"), sum(status$status == "skipped"),
              sum(status$status == "failed"), unname(version_of)),
    stringsAsFactors = FALSE)
}

#' The declared version and a fingerprint of the package source
#'
#' The fingerprint is one md5 over the per-file md5 sums of `R/`, so a single
#' string states whether the estimation code changed between two moments. It
#' uses only base R and `tools`.
#'
#' @param root Project root directory.
#' @return A list with `version` and `fingerprint`, both single strings.
.source_state <- function(root = ".") {
  sources <- sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))
  fingerprint <- if (length(sources) == 0L) "no R/ directory" else {
    manifest <- tempfile("multilpa-source-")
    on.exit(unlink(manifest), add = TRUE, after = FALSE)
    writeLines(paste(basename(sources), tools::md5sum(sources)), manifest)
    unname(tools::md5sum(manifest))
  }
  version <- tryCatch(as.character(utils::packageVersion("multilpa")),
                      error = function(condition) {
                        # The DESCRIPTION is the only place the version lives;
                        # if it cannot be read, say so rather than invent one.
                        sprintf("unreadable (%s)", conditionMessage(condition))
                      })
  list(version = version, fingerprint = fingerprint)
}

# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------

#' Reach a table of an equivalence run
#'
#' @param x An object returned by [run_equivalence()].
#' @param row.names Passed through for consistency with the generic; ignored.
#' @param optional Passed through for consistency with the generic; ignored.
#' @param what Which table: `"comparisons"` is one row per compared quantity,
#'   `"suites"` one row per suite with its status and reason, `"session"` the
#'   versions the run was produced under, and `"summary"` the by-source
#'   summary.
#' @param ... Unused.
#' @return A `data.frame`.
as.data.frame.multilpa_equivalence <- function(x, row.names = NULL, optional = FALSE,
                                               what = c("comparisons", "suites",
                                                        "session", "summary"), ...) {
  what <- match.arg(what)
  switch(what,
         comparisons = x$comparisons,
         suites = x$suites,
         session = x$session,
         summary = equivalence_summary(x))
}

#' One recorded provenance value of an equivalence run
#'
#' The version a report is evidence for is needed as a scalar, in a filename, a
#' header line and a column. Reaching into the session table for it at each
#' call site is the ritual this verb exists to remove.
#'
#' @param run An object returned by [run_equivalence()].
#' @param component Which component to report: `"multilpa"`, `"R"`,
#'   `"platform"`, `"generated"`, or any external package the run recorded.
#' @return A single character string.
equivalence_provenance <- function(run, component = "multilpa") {
  stopifnot(
    "`run` must come from run_equivalence()" = inherits(run, "multilpa_equivalence"),
    "`component` must be a single string" =
      is.character(component) && length(component) == 1L)
  matched <- match(component, run$session$component)
  if (is.na(matched)) {
    stop(errorCondition(
      sprintf("`component` must be one of: %s",
              paste(run$session$component, collapse = ", ")),
      class = "multilpa_bad_argument", call = NULL))
  }
  run$session$value[[matched]]
}

#' Print an equivalence run
#' @param x An object returned by [run_equivalence()].
#' @param ... Unused.
#' @return `x`, invisibly.
print.multilpa_equivalence <- function(x, ...) {
  cat(sprintf("Equivalence run against multilpa %s under %s\n",
              equivalence_provenance(x, "multilpa"), R.version.string))
  print(x$suites[, c("suite", "status", "compared", "agreed", "warnings", "seconds")])
  cat(sprintf("\n%d quantities compared, %d agreed, %d disagreed\n",
              nrow(x$comparisons), sum(x$comparisons$agrees),
              sum(!x$comparisons$agrees)))
  skipped <- x$suites[x$suites$status != "ran", , drop = FALSE]
  if (nrow(skipped) > 0L) {
    cat("\nSuites not run:\n")
    print(skipped[, c("suite", "status", "reason")])
  }
  invisible(x)
}

#' Summarise an equivalence report by source
#'
#' @param report An object returned by [run_equivalence()], or its comparisons
#'   frame.
#' @return A `data.frame` with one row per suite, source and reference
#'   precision, giving the number of quantities compared, how many agreed, the
#'   largest difference seen and the loosest tolerance claimed.
equivalence_summary <- function(report) {
  if (inherits(report, "multilpa_equivalence")) report <- report$comparisons
  stopifnot("`report` must be an equivalence report" =
              is.data.frame(report) && all(c("suite", "source", "agrees") %in% names(report)))
  if (nrow(report) == 0L) {
    return(data.frame(suite = character(), source = character(),
                      precision = character(), compared = integer(),
                      agreed = integer(), worst_difference = numeric(),
                      loosest_tolerance = numeric(), stringsAsFactors = FALSE))
  }
  parts <- split(report, list(report$suite, report$source, report$precision),
                 drop = TRUE)
  summaries <- lapply(parts, function(part) {
    data.frame(suite = part$suite[[1L]], source = part$source[[1L]],
               precision = part$precision[[1L]],
               compared = nrow(part), agreed = sum(part$agrees),
               worst_difference = max(part$difference),
               loosest_tolerance = max(part$tolerance),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, summaries)
  row.names(out) <- NULL
  out[order(out$precision, out$suite, out$source), ]
}

# ---------------------------------------------------------------------------
# Helpers the suites share
# ---------------------------------------------------------------------------

#' Stack a per-indicator list of class-by-category matrices into one matrix
#'
#' @param blocks List of matrices sharing a row count, one per indicator.
#' @return A matrix with one row per latent class and every category in columns.
flatten_blocks <- function(blocks) {
  stopifnot("`blocks` must be a list of matrices" =
              is.list(blocks) && all(vapply(blocks, is.matrix, logical(1))))
  do.call(cbind, lapply(blocks, function(block) unname(as.matrix(block))))
}

#' Align two labelings of the same latent classes
#'
#' Latent class labels are arbitrary, so two programs fitting the same model
#' recover the same classes under some permutation. This finds the permutation
#' minimising total absolute distance, which is exact for the class counts used
#' here rather than a greedy match.
#'
#' @param target Matrix with one row per class, in the order to match.
#' @param candidate Matrix of the same dimensions in an unknown class order.
#' @return An integer permutation `p` such that `candidate[p, ]` aligns with
#'   `target`.
align_classes <- function(target, candidate) {
  target <- unname(as.matrix(target))
  candidate <- unname(as.matrix(candidate))
  stopifnot("`target` and `candidate` must have the same dimensions" =
              identical(dim(target), dim(candidate)))
  orders <- .class_permutations(nrow(target))
  distances <- vapply(seq_len(nrow(orders)), function(row) {
    sum(abs(target - candidate[orders[row, ], , drop = FALSE]))
  }, numeric(1))
  orders[which.min(distances), ]
}

#' Every permutation of seq_len(n), one per row
#' @param n Number of elements; kept small because this enumerates n! rows.
#' @return An integer matrix with `factorial(n)` rows.
.class_permutations <- function(n) {
  stopifnot("permutation search is only used for small class counts" = n <= 8L)
  if (n == 1L) return(matrix(1L, 1L, 1L))
  smaller <- .class_permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(value) {
    cbind(value, matrix(setdiff(seq_len(n), value)[smaller], nrow(smaller)),
          deparse.level = 0)
  }))
}

#' Refit a model across several seeds and report what varied
#'
#' A single-seed agreement is not a finding: every suite here starts its
#' optimiser from random values, and a comparison that passes under one seed
#' and fails under another is measuring the seed. This refits the same model
#' under each of `seeds` and reports two quantities: the spread of the
#' maximised log-likelihood, which should be zero if every seed reaches the
#' same mode, and the number of distinct optima reached.
#'
#' @param label A name for the model, used in the reported quantities.
#' @param fit_one A function of one argument, a seed, returning a fitted
#'   `multilpa` model.
#' @param seeds Integer vector of seeds.
#' @param tolerance How far apart two maximised log-likelihoods may be and
#'   still count as the same mode.
#' @return A `data.frame` shaped like [compare_values()], with one row for the
#'   spread and one for the number of distinct optima.
seed_stability <- function(label, fit_one, seeds, tolerance = .equivalence_tolerances[["likelihood"]]) {
  stopifnot(
    "`label` must be a single string" = is.character(label) && length(label) == 1L,
    "`fit_one` must be a function of a seed" = is.function(fit_one),
    "`seeds` must be at least three seeds" = is.numeric(seeds) && length(seeds) >= 3L)
  optima <- vapply(seeds, function(seed) as.numeric(stats::logLik(fit_one(seed))), numeric(1))
  spread <- max(optima) - min(optima)
  distinct <- length(unique(round(optima / tolerance)))
  rbind(
    compare_values(sprintf("%s: spread of the maximised log-likelihood across %d seeds",
                           label, length(seeds)),
                   reference = 0, obtained = spread, tolerance = tolerance),
    compare_values(sprintf("%s: distinct optima reached across %d seeds",
                           label, length(seeds)),
                   reference = 1, obtained = distinct,
                   tolerance = .equivalence_tolerances[["count"]]))
}

# ---------------------------------------------------------------------------
# Writing the report
# ---------------------------------------------------------------------------

#' Render an equivalence run as markdown
#'
#' @param run The object returned by [run_equivalence()].
#' @return A character vector of markdown lines.
.equivalence_markdown <- function(run) {
  stopifnot("`run` must come from run_equivalence()" =
              inherits(run, "multilpa_equivalence"))
  report <- run$comparisons
  summary_table <- equivalence_summary(run)
  row_text <- function(values) paste0("| ", paste(values, collapse = " | "), " |")
  version <- equivalence_provenance(run, "multilpa")
  header <- c(
    "# Equivalence report",
    "",
    sprintf("Generated by `Rscript equivalence/harness/run.R` against **multilpa %s** under %s on %s.",
            version, R.version.string,
            equivalence_provenance(run, "generated")),
    "",
    "The version above is the version this report is evidence for. It is",
    "recorded by the harness, not typed in, so the report cannot be attributed",
    "to a release it was not produced under.",
    "",
    if (startsWith(equivalence_provenance(run, "coverage"), "PARTIAL")) {
      sprintf("**This report covers only part of the harness: %s.** Run `Rscript equivalence/harness/run.R` with no arguments to regenerate complete evidence.",
              equivalence_provenance(run, "coverage"))
    } else {
      sprintf("This report covers %s.", equivalence_provenance(run, "coverage"))
    },
    "",
    if (identical(equivalence_provenance(run, "source unchanged during run"), "yes")) {
      "The package source was unchanged for the whole run."
    } else {
      sprintf(paste("**The package source changed while this run was in progress.**",
                    "It opened at version %s with `R/` fingerprint `%s` and closed at",
                    "version %s with fingerprint `%s`. The estimation code used is the",
                    "snapshot loaded at the start; re-run the harness once the tree is",
                    "settled before citing these numbers."),
              equivalence_provenance(run, "multilpa"),
              equivalence_provenance(run, "R/ fingerprint"),
              equivalence_provenance(run, "multilpa at end"),
              equivalence_provenance(run, "R/ fingerprint at end"))
    },
    "",
    sprintf("**%d quantities compared across %d sources; %d agreed within tolerance.**",
            nrow(report), length(unique(report$source)), sum(report$agrees)),
    "",
    sprintf("Largest difference against a reference computed here at full double precision: **%.3e**.",
            max(c(report$difference[report$precision == "machine"], 0))),
    sprintf("Largest gap against a reference read from published digits: **%.3g** of the rounding that value's printed precision admits, where 1 would be the whole of it. Published references are compared this way because the raw gap against a statistic printed to one decimal measures the typesetting, not the implementation.",
            max(c(report$difference[report$precision == "printed"], 0))),
    "",
    "## Suites",
    "",
    row_text(c("Suite", "Status", "Compared", "Agreed", "Warnings", "Seconds", "Reason")),
    row_text(rep("---", 7L)),
    vapply(seq_len(nrow(run$suites)), function(row) {
      row_text(c(run$suites$suite[[row]], run$suites$status[[row]],
                 run$suites$compared[[row]], run$suites$agreed[[row]],
                 run$suites$warnings[[row]], run$suites$seconds[[row]],
                 if (is.na(run$suites$reason[[row]])) "" else run$suites$reason[[row]]))
    }, character(1)),
    "",
    "## By source",
    "",
    row_text(c("Suite", "Reference source", "Reference precision", "Dataset(s)",
               "Compared", "Agreed", "Largest difference", "Loosest tolerance")),
    row_text(rep("---", 8L)))
  datasets <- vapply(seq_len(nrow(summary_table)), function(row) {
    part <- report[report$suite == summary_table$suite[[row]] &
                     report$source == summary_table$source[[row]] &
                     report$precision == summary_table$precision[[row]], , drop = FALSE]
    paste(sort(unique(stats::na.omit(part$dataset))), collapse = ", ")
  }, character(1))
  body <- vapply(seq_len(nrow(summary_table)), function(row) {
    row_text(c(summary_table$suite[[row]], summary_table$source[[row]],
               summary_table$precision[[row]], datasets[[row]],
               summary_table$compared[[row]], summary_table$agreed[[row]],
               sprintf("%.2e", summary_table$worst_difference[[row]]),
               sprintf("%.2e", summary_table$loosest_tolerance[[row]])))
  }, character(1))
  failures <- report[!report$agrees, , drop = FALSE]
  tail_lines <- if (nrow(failures) == 0L) {
    c("", "Every compared quantity agreed within its stated tolerance.")
  } else {
    c("", "## Disagreements", "",
      row_text(c("Suite", "Quantity", "Reference", "Obtained", "Difference",
                 "Tolerance")),
      row_text(rep("---", 6L)),
      vapply(seq_len(nrow(failures)), function(row) {
        row_text(c(failures$suite[[row]], failures$quantity[[row]],
                   sprintf("%.6g", failures$reference[[row]]),
                   sprintf("%.6g", failures$obtained[[row]]),
                   sprintf("%.2e", failures$difference[[row]]),
                   sprintf("%.2e", failures$tolerance[[row]])))
      }, character(1)))
  }
  session_lines <- c(
    "", "## Produced under", "",
    row_text(c("Component", "Value")), row_text(rep("---", 2L)),
    vapply(seq_len(nrow(run$session)), function(row) {
      row_text(c(run$session$component[[row]], run$session$value[[row]]))
    }, character(1)),
    "", "The full `sessionInfo()` is in `SESSION.txt`.")
  c(header, body, tail_lines, session_lines, "",
    "Per-quantity detail, including every tolerance, is in `report.csv`;",
    "per-suite status, including anything skipped, is in `suites.csv`.")
}
