# Report test and equivalence results to CVS, the validation service in
# ../validation (docs/AUTOMATION.md there). Each run becomes one JSON file in
# the vitest layout CVS ingests: one entry per source file, one assertion per
# test or compared quantity, with the difference and tolerance where the run
# recorded them. Each file declares its CVS tier (`_cvs.tier`): CVS classifies
# by file name otherwise, and an R file name says nothing about the tier.
#
# The file goes to the CVS inbox when reporting is on (mini-ci sets
# CVS_INBOX_DIR for jobs of a repo whose .mini-ci says `report cvs`; it sets
# CVS_DISABLED=1 otherwise), and a copy always goes to tmp/cvs/ so a run that
# could not report keeps its file. The file name starts with "latents-", which
# is how the CVS server assigns it to the latents target.
#
# Every assertion is tagged with the reference it is evidence for (`_cvs`
# axis, authority, reference), from the verification matrix's own attribution
# rules (verification-matrix/attribution.R), so CVS's per-target Verification
# tab and the matrix count the same evidence the same way
# (../validation/docs/TARGET-VERIFICATION.md).
#
# Used by equivalence/run.R (testthat equivalence files), equivalence/cvs/
# test-local.R (the package's own tests), equivalence/harness/run.R and
# equivalence/argument-sweep/run.R. Needs jsonlite
# (equivalence/install-reference-packages.R).

#' The verification-matrix folder, from the project root or from equivalence/
#' (where testthat runs this folder's tests)
#' @noRd
.cvs_matrix_dir <- function() {
  candidates <- c(file.path("equivalence", "verification-matrix"), "verification-matrix")
  found <- candidates[file.exists(file.path(candidates, "attribution.R"))]
  if (!length(found)) {
    stop(errorCondition(
      "verification-matrix/attribution.R not found: run from the project root or equivalence/",
      class = "latents_cvs_no_attribution", call = NULL))
  }
  found[[1L]]
}
source(file.path(.cvs_matrix_dir(), "attribution.R"), local = TRUE)

#' The report entry (file name) of a harness suite and of a sweep module
#' @noRd
.cvs_harness_entry <- function(suite) file.path("equivalence", "harness", "suites", paste0(suite, ".R"))
.cvs_sweep_entry <- function(module) file.path("equivalence", "argument-sweep", "modules", paste0(module, ".R"))

# CVS's reference axes (TARGET-VERIFICATION.md, "Tag every assertion") for
# each kind of reference in the matrix's rule tables. A commercial program is
# told apart by the program its reference names ("Mplus 9 Demo (JStats
# fixtures)", "IMMERSE lca-bookdown (Mplus 8.11)" are Mplus). Provenance and infrastructure
# checks (kit completeness, the reporter's own tests) are latents' internal
# checks: tagged internal, so CVS shows them and never counts them as a source,
# rather than left untagged, which CVS would count as R. "Python package" has
# no CVS axis in latents' configuration and stays untagged.
.cvs_axis_map <- data.frame(
  kind = c("R package", "commercial program", "commercial program", "commercial program",
           "published", "independent code", "internal", "component check", "provenance",
           "infrastructure"),
  program = c(NA, "Latent GOLD", "Mplus", "Stata", NA, NA, NA, NA, NA, NA),
  axis = c("r-equivalence", "latentgold", "mplus", "stata", "paper", "independent",
           "internal", "internal", "internal", "internal"),
  authority = c("r", "latentgold-only", "mplus-only", "stata-only", "paper", "independent",
                "internal", "internal", "internal", "internal"),
  software = c(NA, "latentgold", "Mplus", "stata", "paper", "independent", NA, NA, NA, NA),
  stringsAsFactors = FALSE)

#' The CVS axis of each reference
#'
#' @param reference,kind A reference string and its kind from the rule tables,
#'   of equal length.
#' @return A data.frame (one row per reference): `reference`, `kind`, `axis`,
#'   `authority`, `software` (the manifest's `cross_software_expected` name, NA
#'   for R and internal). `axis` and `authority` are NA for a kind CVS has no
#'   axis for.
#' @section Errors: a commercial program whose reference names none, or more
#'   than one, of Latent GOLD, Mplus and Stata raises `latents_cvs_unknown_program`.
#'   It is not left untagged: the rule tables are fixed files, and the
#'   reporter's tests check that every commercial program in them maps, so a
#'   new program fails those tests (seconds) rather than a run; untagged, its
#'   evidence would be counted by CVS as R.
#' @noRd
.cvs_axis <- function(reference, kind) {
  stopifnot("`reference` and `kind` are character vectors of equal length" =
              is.character(reference) && is.character(kind) && length(reference) == length(kind))
  row <- vapply(seq_along(kind), function(i) {
    candidates <- which(.cvs_axis_map$kind == kind[[i]])
    programs <- .cvs_axis_map$program[candidates]
    named <- is.na(programs) | vapply(as.character(programs), function(program) {
      !is.na(program) && grepl(program, reference[[i]], fixed = TRUE)
    }, logical(1))
    hit <- candidates[named]
    if (identical(kind[[i]], "commercial program") && length(hit) != 1L) {
      stop(errorCondition(sprintf(
        "commercial program \"%s\" names no single program with a CVS axis: extend .cvs_axis_map",
        reference[[i]]),
        class = "latents_cvs_unknown_program", call = NULL))
    }
    if (length(hit)) hit[[1L]] else NA_integer_
  }, integer(1))
  data.frame(reference = reference, kind = kind, axis = .cvs_axis_map$axis[row],
             authority = .cvs_axis_map$authority[row], software = .cvs_axis_map$software[row],
             stringsAsFactors = FALSE)
}

#' The attribution rule tables the reporter reads
#' @param dir The verification-matrix folder.
#' @return A list: `harness` (harness-references.csv), `tests` (test-attribution.csv).
#' @noRd
cvs_rules <- function(dir = .cvs_matrix_dir()) {
  list(harness = matrix_read_rules("harness-references.csv", dir),
       tests = matrix_read_rules("test-attribution.csv", dir))
}

#' The references each assertion of one report entry compares with
#'
#' Harness quantities by their source (`_cvs$rPackage`), testthat blocks by the
#' rule of test-attribution.csv their file and name match (a numbered repeat,
#' "name [2]", by its name). Sweep checks, the harness's coverage and
#' known-open rows, and tests no rule names have none.
#' @return A list with one element per assertion: a data.frame of
#'   `reference` and `kind` (zero rows when unattributed).
#' @noRd
.cvs_entry_references <- function(entry, rules) {
  none <- data.frame(reference = character(), kind = character(), stringsAsFactors = FALSE)
  assertions <- entry$assertionResults
  if (startsWith(entry$name, file.path("equivalence", "argument-sweep", ""))) {
    return(rep(list(none), length(assertions)))
  }
  if (startsWith(entry$name, file.path("equivalence", "harness", "suites", ""))) {
    sources <- vapply(assertions, function(a) a[["_cvs"]]$rPackage %||% NA_character_, character(1))
    hit <- rep(NA_integer_, length(sources))
    hit[!is.na(sources)] <- matrix_harness_reference(sources[!is.na(sources)], rules$harness)
    return(lapply(hit, function(h) {
      if (is.na(h)) none else data.frame(reference = rules$harness$reference[[h]],
                                         kind = rules$harness$kind[[h]], stringsAsFactors = FALSE)
    }))
  }
  titles <- sub(" \\[[0-9]+\\]$", "", vapply(assertions, `[[`, character(1), "title"))
  hit <- matrix_test_rule(dirname(entry$name), basename(entry$name), titles, rules$tests)
  lapply(hit, function(h) {
    if (is.na(h)) none else matrix_rule_references(rules$tests$reference[[h]], rules$tests$kind[[h]])
  })
}

#' Tag every assertion with the reference it is evidence for
#'
#' Adds `axis`, `authority` and `reference` to each assertion's `_cvs` when
#' its references all fall on one CVS axis. An assertion comparing with
#' references on several axes (one block checking Stata and a published table)
#' gets `reference` only: CVS takes one axis per assertion, and choosing one
#' would be a guess. An assertion with no attributed reference is left as it is.
#' @param suites A list of report entries from the `.cvs_from_*()` builders.
#' @param rules The rule tables (`cvs_rules()`).
#' @return `suites`, tagged.
#' @noRd
cvs_tag_suites <- function(suites, rules = cvs_rules()) {
  stopifnot("`suites` must be a list of report entries" = is.list(suites),
            "`rules` must hold the harness and test tables" =
              all(c("harness", "tests") %in% names(rules)))
  lapply(suites, function(entry) {
    references <- .cvs_entry_references(entry, rules)
    entry$assertionResults <- Map(function(assertion, refs) {
      if (!nrow(refs)) return(assertion)
      axes <- .cvs_axis(refs$reference, refs$kind)
      tag <- list(reference = paste(refs$reference, collapse = "; "))
      if (!anyNA(axes$axis) && length(unique(axes$axis)) == 1L) {
        tag <- c(list(axis = axes$axis[[1L]], authority = axes$authority[[1L]]), tag)
      }
      assertion[["_cvs"]] <- utils::modifyList(assertion[["_cvs"]] %||% list(), tag)
      assertion
    }, entry$assertionResults, references)
    entry
  })
}

#' Assertions per CVS axis
#' @param suites Tagged report entries.
#' @return A data.frame, one row per axis: `axis` and `assertions`. Assertions
#'   without an axis are counted as "(several axes)" when they name references
#'   and "(unattributed)" when they do not. Sorted by axis.
#' @noRd
cvs_axis_counts <- function(suites) {
  axis <- unlist(lapply(suites, function(entry) {
    vapply(entry$assertionResults, function(assertion) {
      tag <- assertion[["_cvs"]]
      tag$axis %||% if (is.null(tag$reference)) "(unattributed)" else "(several axes)"
    }, character(1))
  }))
  counts <- table(factor(axis, levels = sort(unique(axis))))
  data.frame(axis = names(counts), assertions = as.integer(counts), stringsAsFactors = FALSE)
}

#' One CVS assertion
#' @param group Title of the enclosing block (the suite or test file).
#' @param title Assertion title; unique within its file.
#' @param status One of "passed", "failed", "skipped".
#' @param failure Failure messages (character, possibly empty).
#' @param meta Named list for `_cvs` (delta, tolerance, rFunction, rPackage,
#'   module), or NULL. Non-finite numbers are dropped, since JSON has none.
#' @return A list in vitest's assertionResults layout.
#' @noRd
.cvs_assertion <- function(group, title, status, failure = character(), meta = NULL) {
  stopifnot(
    "`status` must be passed, failed or skipped" =
      length(status) == 1L && status %in% c("passed", "failed", "skipped"))
  out <- list(ancestorTitles = list(group), title = title,
              fullName = sprintf("%s > %s", group, title), status = status,
              duration = 0, failureMessages = as.list(failure))
  if (!is.null(meta)) {
    keep <- vapply(meta, function(value) {
      length(value) == 1L && !is.na(value) && (!is.numeric(value) || is.finite(value))
    }, logical(1))
    out[["_cvs"]] <- meta[keep]
  }
  out
}

#' Make titles unique within a file by numbering repeats
#' @noRd
.cvs_unique <- function(titles) {
  repeat_number <- stats::ave(seq_along(titles), titles, FUN = seq_along)
  ifelse(repeat_number > 1L, sprintf("%s [%d]", titles, repeat_number), titles)
}

#' Suites from a testthat run
#'
#' One suite per test file, one assertion per `test_that()` block. A block is
#' failed when any expectation failed or errored, skipped when it skipped
#' (even after passing expectations, so a partial check stays visible), and
#' passed otherwise.
#' @param results The `testthat_results` returned by `testthat::test_dir()`.
#' @param directory The directory the files were run from, prefixed to names.
#' @param tier CVS tier of these files: "equiv" for comparisons with references,
#'   "unit" for the package's own tests. CVS cannot tell from an R file name.
#' @return A list of suites.
#' @noRd
.cvs_from_testthat <- function(results, directory, tier) {
  # No default: a defaulted tier once filed the package's own tests as equivalence.
  stopifnot(
    "`results` must be testthat results" = inherits(results, "testthat_results"),
    "`tier` must be \"equiv\" or \"unit\"" = is.character(tier) && length(tier) == 1L &&
      tier %in% c("equiv", "unit"))
  tests <- lapply(results, function(test) {
    kinds <- vapply(test$results, function(expectation) class(expectation)[[1L]], character(1))
    problems <- test$results[kinds %in% c("expectation_failure", "expectation_error")]
    skips <- test$results[kinds == "expectation_skip"]
    list(file = test$file, test = test$test,
         status = if (length(problems)) "failed" else if (length(skips)) "skipped" else "passed",
         messages = vapply(c(problems, skips), conditionMessage, character(1)))
  })
  files <- vapply(tests, `[[`, character(1), "file")
  unname(lapply(split(tests, factor(files, levels = unique(files))), function(file_tests) {
    file <- file_tests[[1L]]$file
    # The title is the test's name alone: CVS follows an assertion across runs by it.
    names <- .cvs_unique(vapply(file_tests, `[[`, character(1), "test"))
    # A test_that() block has no single difference to report; checkKind says what the
    # row is (CVS otherwise reads a run with no numeric metadata as a capture outage).
    module <- gsub("^test-|[.][Rr]$", "", file)
    assertions <- Map(function(test, name) {
      .cvs_assertion(file, name, test$status, test$messages,
                     list(checkKind = "testthat", module = module))
    }, file_tests, names)
    list(name = file.path(directory, file), `_cvs` = list(tier = tier),
         assertionResults = unname(assertions))
  }))
}

#' Suites from the harness comparisons (equivalence/harness/run.R)
#'
#' One suite per harness suite, one assertion per compared quantity. A
#' quantity that agrees passes; one that does not fails, also when a known
#' open defect explains it, since the defect is real (the message names it).
#' The other reasons the harness fails a run become failed assertions of
#' their suite too: an output it produced but did not completely compare
#' ("output compared: <fit>: <output>", from coverage.csv) and a known-open
#' entry that matches nothing ("known-open entry <id> matches", from
#' known-open-status.csv). A suite that was skipped or failed before
#' comparing anything is one assertion, "(suite ran)", with the reason.
#' @param comparisons The comparisons table the harness writes as report.csv.
#' @param suites The suites table (suites.csv).
#' @param coverage The coverage table (coverage.csv), or NULL when the run
#'   kept none (runs before the coverage ledger).
#' @param known_open The known-open status table (known-open-status.csv), or NULL.
#' @return A list of suites.
#' @noRd
.cvs_from_harness <- function(comparisons, suites, coverage = NULL, known_open = NULL) {
  needed <- c("suite", "source", "quantity", "difference", "tolerance", "agrees",
              "output", "n_values", "known_open")
  stopifnot(
    "`comparisons` lacks harness columns" = all(needed %in% names(comparisons)),
    "`suites` lacks harness columns" = all(c("suite", "status", "reason") %in% names(suites)),
    "`coverage` lacks harness columns" = is.null(coverage) ||
      all(c("suite", "fit", "output", "n_produced", "n_compared", "status", "reason") %in% names(coverage)),
    "`known_open` lacks harness columns" = is.null(known_open) ||
      all(c("id", "suite", "status") %in% names(known_open)))
  # Each source gives a data.frame of (suite, assertion); they are grouped by suite at the end.
  quantity_rows <- lapply(split(comparisons, factor(comparisons$suite, levels = unique(comparisons$suite))),
                          function(rows) {
    suite <- rows$suite[[1L]]
    titles <- .cvs_unique(rows$quantity)
    lapply(seq_len(nrow(rows)), function(i) {
      agrees <- isTRUE(rows$agrees[[i]])
      failure <- if (agrees) character() else sprintf(
        "difference %.3g exceeds tolerance %.3g against %s%s", rows$difference[[i]],
        rows$tolerance[[i]], rows$source[[i]],
        if (is.na(rows$known_open[[i]])) "" else sprintf(" (known open defect %s)", rows$known_open[[i]]))
      list(suite = suite, assertion = .cvs_assertion(
        suite, titles[[i]], if (agrees) "passed" else "failed", failure,
        list(delta = abs(rows$difference[[i]]), tolerance = rows$tolerance[[i]],
             rFunction = rows$output[[i]], rPackage = rows$source[[i]],
             module = suite, values = rows$n_values[[i]])))
    })
  })
  gaps <- if (is.null(coverage)) coverage else coverage[coverage$status %in% c("partial", "not covered"), ]
  gap_rows <- lapply(seq_len(NROW(gaps)), function(i) {
    list(suite = gaps$suite[[i]], assertion = .cvs_assertion(
      gaps$suite[[i]], sprintf("output compared: %s: %s", gaps$fit[[i]], gaps$output[[i]]), "failed",
      sprintf("%s: %d of %d values compared%s", gaps$status[[i]], gaps$n_compared[[i]],
              gaps$n_produced[[i]], if (is.na(gaps$reason[[i]])) "" else paste0("; ", gaps$reason[[i]]))))
  })
  stale <- if (is.null(known_open)) known_open else known_open[known_open$status %in% "stale", ]
  stale_rows <- lapply(seq_len(NROW(stale)), function(i) {
    list(suite = stale$suite[[i]], assertion = .cvs_assertion(
      stale$suite[[i]], sprintf("known-open entry %s matches", stale$id[[i]]), "failed",
      "the register explains a disagreement that no longer occurs; remove the entry"))
  })
  not_run <- suites[suites$status != "ran", ]
  not_run_rows <- lapply(seq_len(nrow(not_run)), function(i) {
    list(suite = not_run$suite[[i]], assertion = .cvs_assertion(
      not_run$suite[[i]], "(suite ran)",
      if (identical(not_run$status[[i]], "skipped")) "skipped" else "failed",
      sprintf("%s: %s", not_run$status[[i]], not_run$reason[[i]])))
  })
  rows <- c(unlist(quantity_rows, recursive = FALSE), gap_rows, stale_rows, not_run_rows)
  suite_of <- vapply(rows, `[[`, character(1), "suite")
  unname(lapply(split(rows, factor(suite_of, levels = unique(suite_of))), function(suite_rows) {
    list(name = .cvs_harness_entry(suite_rows[[1L]]$suite),
         `_cvs` = list(tier = "equiv"),
         assertionResults = unname(lapply(suite_rows, `[[`, "assertion")))
  }))
}

#' Suites from the argument sweep (equivalence/argument-sweep/run.R)
#'
#' One suite per module, one assertion per check. "agrees" passes;
#' "disagrees" and "error" fail; "known difference" (registered in
#' known-differences.csv) and "no reference" are skipped with their reason.
#' @param results The results table the sweep writes as results.csv.
#' @return A list of suites.
#' @noRd
.cvs_from_sweep <- function(results) {
  needed <- c("module", "case", "quantity", "difference", "tolerance", "status", "note")
  stopifnot("`results` lacks sweep columns" = all(needed %in% names(results)))
  status_map <- c(agrees = "passed", disagrees = "failed", error = "failed",
                  `known difference` = "skipped", `no reference` = "skipped")
  unknown <- setdiff(unique(results$status), names(status_map))
  if (length(unknown)) {
    stop(sprintf("unknown sweep status: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  unname(lapply(split(results, factor(results$module, levels = unique(results$module))),
                function(rows) {
    module <- rows$module[[1L]]
    titles <- .cvs_unique(sprintf("%s: %s", rows$case, rows$quantity))
    assertions <- lapply(seq_len(nrow(rows)), function(i) {
      status <- status_map[[rows$status[[i]]]]
      reason <- if (is.na(rows$note[[i]])) rows$status[[i]] else
        sprintf("%s: %s", rows$status[[i]], rows$note[[i]])
      .cvs_assertion(module, titles[[i]], status,
                     if (identical(status, "passed")) character() else reason,
                     list(delta = abs(rows$difference[[i]]), tolerance = rows$tolerance[[i]],
                          module = module))
    })
    list(name = .cvs_sweep_entry(module),
         `_cvs` = list(tier = "equiv"), assertionResults = assertions)
  }))
}

#' Where the report goes: the CVS inbox, or NULL when reporting is off
#' @noRd
.cvs_inbox <- function() {
  if (identical(Sys.getenv("CVS_DISABLED"), "1")) return(NULL)
  candidates <- c(Sys.getenv("CVS_INBOX_DIR"), file.path("..", "validation", "data", "inbox"))
  candidates <- path.expand(candidates[nzchar(candidates)])
  usable <- candidates[dir.exists(candidates) & file.access(candidates, 2L) == 0L]
  if (length(usable)) usable[[1L]] else NULL
}

#' Write a CVS report
#'
#' Every assertion is tagged first (`cvs_tag_suites()`).
#' @param suites A list of suites from one of the `.cvs_from_*()` builders.
#' @param label Short name of the run (e.g. "harness"); part of the file name.
#' @param start_time POSIXct the run started.
#' @param output Directory for the local copy.
#' @param version Package version the run tested; the loaded version by default
#'   (a report built later from saved results passes the version they record).
#' @param name File name; by default latents-<label>-<start time>.json.
#' @return The paths written, invisibly.
#' @noRd
cvs_write_report <- function(suites, label, start_time, output = file.path("tmp", "cvs"),
                             version = as.character(utils::packageVersion("latents")),
                             name = sprintf("latents-%s-%s.json", label,
                                            format(start_time, "%Y%m%dT%H%M%S"))) {
  stopifnot(
    "`suites` must be a non-empty list" = is.list(suites) && length(suites) > 0L,
    "`label` must be one lowercase word" = length(label) == 1L && grepl("^[a-z0-9-]+$", label),
    "`start_time` must be a time" = inherits(start_time, "POSIXct") && length(start_time) == 1L,
    "`version` must be one string" = is.character(version) && length(version) == 1L,
    "`name` must be a latents-*.json file name" = length(name) == 1L &&
      grepl("^latents-[^/]+[.]json$", name))
  suites <- cvs_tag_suites(suites)
  counts <- cvs_axis_counts(suites)
  axes <- sprintf("; axes: %s", paste(sprintf("%s %d", counts$axis, counts$assertions),
                                      collapse = ", "))
  statuses <- unlist(lapply(suites, function(suite) {
    vapply(suite$assertionResults, `[[`, character(1), "status")
  }))
  suite_failed <- vapply(suites, function(suite) {
    any(vapply(suite$assertionResults, `[[`, character(1), "status") == "failed")
  }, logical(1))
  report <- list(
    numTotalTestSuites = length(suites), numPassedTestSuites = sum(!suite_failed),
    numFailedTestSuites = sum(suite_failed), numTotalTests = length(statuses),
    numPassedTests = sum(statuses == "passed"), numFailedTests = sum(statuses == "failed"),
    numPendingTests = sum(statuses == "skipped"),
    startTime = round(as.numeric(start_time) * 1000), success = !any(suite_failed),
    testResults = unname(suites),
    `_cvs` = list(target = "latents", run = label, version = version))
  # digits = NA: full precision; a fixed number of significant digits rounds startTime.
  json <- jsonlite::toJSON(report, auto_unbox = TRUE, null = "null", digits = NA)
  # Named by the run's start (the default `name`), so a report built later from saved
  # results is named by its run.
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(output, name)
  writeLines(json, paths[[1L]])
  inbox <- .cvs_inbox()
  if (is.null(inbox)) {
    message(sprintf("CVS report kept at %s (reporting off%s)", paths[[1L]], axes))
  } else {
    # Write under a temporary name and rename, so the watcher never reads half a file.
    partial <- file.path(inbox, paste0(".", name, ".part"))
    writeLines(json, partial)
    stopifnot("could not place the report in the CVS inbox" =
                file.rename(partial, file.path(inbox, name)))
    paths <- c(paths, file.path(inbox, name))
    message(sprintf("CVS report: %s (%d tests, %d failed)%s", file.path(inbox, name),
                    length(statuses), sum(statuses == "failed"), axes))
  }
  invisible(paths)
}

#' Tag a report written before tagging, from its saved file
#'
#' Reads a saved report, tags it (`cvs_tag_suites()`) and writes it again with
#' its own run label, start time, version and file name, so a run need not be
#' repeated.
#' @param path A saved report (latents-<run>-<time>.json).
#' @param output Directory for the local copy.
#' @return The paths written, invisibly (as `cvs_write_report()`).
#' @noRd
cvs_retag_report <- function(path, output = file.path("tmp", "cvs")) {
  stopifnot("`path` must be one existing file" = length(path) == 1L && file.exists(path))
  report <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  stopifnot("`path` is a latents CVS report" = identical(report[["_cvs"]]$target, "latents"))
  started <- as.POSIXct(report$startTime / 1000, origin = "1970-01-01", tz = "UTC")
  cvs_write_report(report$testResults, report[["_cvs"]]$run, started, output = output,
                   version = report[["_cvs"]]$version, name = basename(path))
}
