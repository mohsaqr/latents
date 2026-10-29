# Verification matrix: every major function against every reference it is
# checked by, built from saved evidence only. Nothing is fitted here.
#
#   Rscript equivalence/verification-matrix/build.R      (from the project root)
#
# Inputs are the harness report and the saved testthat results named in
# `.matrix_evidence_files()`, the attribution rules beside this script, and
# `candidates.csv` (references found but not yet wired in). Outputs, written
# beside this script: evidence.csv (one row per check), matrix.csv (one row per
# function and reference), summary.csv (one row per function), and
# VERIFICATION-MATRIX.html / .pdf. The PDF is printed by headless Chrome; set
# LATENTS_CHROME to its executable when it is not at the macOS default.

.matrix_dir <- file.path("equivalence", "verification-matrix")
source(file.path(.matrix_dir, "provenance.R"), local = TRUE)

.matrix_functions <- c(
  lpa = "lpa()", multilpa = "multilpa()", lca = "lca()", multilca = "multilca()",
  lta = "lta()", mixture_regression = "mixture_regression()",
  growth = "mixture_regression(random =)", three_step = "three_step() / r3step()",
  bootstrap_lrt = "bootstrap_lrt() / VLMR", enumeration = "enumerate_*() / compare_models()",
  pool_imputations = "pool_imputations()")

.matrix_external_kinds <- c("commercial program", "R package", "Python package")

# A check on a one-class reduction (one class equals lm(), glm.nb(), lme4, a
# linear mixed model...) confirms plumbing, not the mixture. Such checks are
# shown but do not make their reference a source. "One class per group" is
# group-level classification, a genuine mixture.
.matrix_one_class_pattern <- "(?i)(^|[^a-z])(one[- ]class(?! per)|1 class|one profile)([^a-z]|$)"
.matrix_kind_order <- c("commercial program", "R package", "Python package", "published",
                        "component check", "independent code", "internal")

#' Saved evidence the matrix reads, with the run it records
#' @return A data.frame: one row per input file.
#' @noRd
.matrix_evidence_files <- function() {
  files <- data.frame(
    role = c("harness", "tests", "tests", "tests", "tests"),
    dir = c(NA, "tests/equivalence", "equivalence", "equivalence", "equivalence"),
    path = file.path("validation", c(
      "release-fixes-2026-10-28/equivalence/report.csv",
      "release-fixes-2026-10-28/release/external-tests-consolidated.csv",
      "release-fixes-2026-10-28/release/retained-references-results.csv",
      "matrix-refresh-2026-10-28/equivalence-results.csv",
      "release-0.96.1-2026-10-20/outer-final/cases.csv")),
    run = c("2026-10-28", "2026-10-28", "2026-10-28", "2026-10-29", "2026-10-20"),
    code = c("0.98.0", "0.98.0", "0.98.0", "0.98.1", "0.96.1"),
    current = FALSE)
  files$current <- vapply(files$path, matrix_evidence_current, logical(1))
  files
}

#' Read one rule table beside this script
#' @noRd
.matrix_rules <- function(name) {
  utils::read.csv(file.path(.matrix_dir, name), stringsAsFactors = FALSE)
}

#' Index of the first pattern each text matches (NA when none does)
#' @noRd
.matrix_first_match <- function(text, patterns) {
  vapply(text, function(value) {
    hits <- which(vapply(patterns, function(pattern) grepl(pattern, value, perl = TRUE),
                         logical(1)))
    if (length(hits)) hits[[1L]] else NA_integer_
  }, integer(1), USE.NAMES = FALSE)
}

#' One evidence row per harness quantity
#' @noRd
.matrix_harness <- function(files) {
  source <- files[files$role == "harness", ]
  report <- utils::read.csv(source$path, stringsAsFactors = FALSE)
  fn_rules <- .matrix_rules("harness-functions.csv")
  ref_rules <- .matrix_rules("harness-references.csv")
  fn_index <- vapply(seq_len(nrow(report)), function(i) {
    candidates <- which(fn_rules$suite == report$suite[[i]])
    hit <- .matrix_first_match(report$quantity[[i]], fn_rules$pattern[candidates])
    if (is.na(hit)) NA_integer_ else candidates[[hit]]
  }, integer(1))
  ref_index <- .matrix_first_match(report$source, ref_rules$pattern)
  stopifnot(
    "every harness quantity is attributed to a function" = !anyNA(fn_index),
    "every harness source is mapped to a reference" = !anyNA(ref_index))
  data.frame(
    evidence_id = paste0("harness:", report$suite, ":", seq_len(nrow(report))),
    origin = "harness", location = paste0("harness/suites/", report$suite, ".R"),
    item = report$quantity, fn = fn_rules$fn[fn_index],
    reference = ref_rules$reference[ref_index], kind = ref_rules$kind[ref_index],
    checks = as.integer(report$agrees), values = report$n_values,
    failed = as.integer(!report$agrees), incomplete = 0L,
    run = source$run, code = source$code, current = source$current,
    dataset = ifelse(is.na(report$dataset) | report$dataset == "NA", NA_character_, report$dataset),
    stringsAsFactors = FALSE)
}

#' Saved testthat results: the current run where one exists, else the older
#' @noRd
.matrix_test_results <- function(files) {
  tests <- files[files$role == "tests" & file.exists(files$path), ]
  pieces <- lapply(seq_len(nrow(tests)), function(i) {
    rows <- utils::read.csv(tests$path[[i]], stringsAsFactors = FALSE)
    outcomes <- c("passed", "failed", "skipped", "error")
    complete <- all(outcomes %in% names(rows))
    passed <- if (complete) as.integer(rows$passed) else rep(0L, nrow(rows))
    errors <- if (complete) as.integer(rows$error) else rep(0L, nrow(rows))
    incomplete <- if (complete) as.integer(rows$skipped | rows$error) else rep(1L, nrow(rows))
    invalid <- is.na(passed) | is.na(errors) | is.na(incomplete) | is.na(rows$failed)
    passed[invalid] <- 0L
    incomplete[invalid] <- 1L
    failures <- as.integer(rows$failed) + errors
    failures[is.na(failures)] <- 0L
    data.frame(dir = tests$dir[[i]], file = rows$file, test = rows$test,
               nb = passed, failed = failures, incomplete = incomplete,
               run = tests$run[[i]], code = tests$code[[i]], current = tests$current[[i]],
               stringsAsFactors = FALSE)
  })
  results <- do.call(rbind, pieces)
  key <- paste(results$dir, results$file)
  # Prefer current evidence, then the latest run. Never combine an old and
  # replacement result for the same file just because both match the source.
  rank <- ifelse(results$current, "1", "0")
  priority <- paste(rank, results$run)
  newest <- tapply(priority, key, max)
  keep <- priority == newest[key]
  out <- results[keep, ]
  row.names(out) <- NULL
  out
}

#' Expand semicolon-separated functions and references into rows
#' @noRd
.matrix_expand <- function(base, fn, reference, kind) {
  fns <- strsplit(fn, ";", fixed = TRUE)[[1L]]
  refs <- strsplit(reference, ";", fixed = TRUE)[[1L]]
  kinds <- strsplit(kind, ";", fixed = TRUE)[[1L]]
  stopifnot("one kind per reference" = length(kinds) %in% c(1L, length(refs)))
  kinds <- rep_len(kinds, length(refs))
  grid <- expand.grid(f = seq_along(fns), r = seq_along(refs))
  data.frame(base[rep(1L, nrow(grid)), , drop = FALSE], fn = fns[grid$f],
             reference = refs[grid$r], kind = kinds[grid$r], shared = nrow(grid) > 1L,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' One evidence row per saved test (and per registered file with no result)
#' @noRd
.matrix_tests <- function(files) {
  rules <- .matrix_rules("test-attribution.csv")
  on_disk <- data.frame(
    dir = rep(c("tests/equivalence", "equivalence"), times = c(
      length(list.files("tests/equivalence", "^test-.*\\.R$")),
      length(list.files("equivalence", "^test-.*\\.R$")))),
    file = c(list.files("tests/equivalence", "^test-.*\\.R$"),
             list.files("equivalence", "^test-.*\\.R$")), stringsAsFactors = FALSE)
  ruled <- unique(paste(rules$dir, rules$file))
  stopifnot(
    "every equivalence test file has an attribution rule" =
      all(paste(on_disk$dir, on_disk$file) %in% ruled),
    "every attribution rule names a file on disk" =
      all(ruled %in% paste(on_disk$dir, on_disk$file)))
  results <- .matrix_test_results(files)
  results <- results[paste(results$dir, results$file) %in% ruled, ]
  saved <- lapply(seq_len(nrow(results)), function(i) {
    own <- which(rules$dir == results$dir[[i]] & rules$file == results$file[[i]])
    hit <- own[[.matrix_first_match(results$test[[i]], rules$pattern[own])]]
    base <- data.frame(
      evidence_id = paste0(results$dir[[i]], "/", results$file[[i]], ":", results$test[[i]]),
      origin = "testthat", location = file.path(results$dir[[i]], results$file[[i]]),
      item = results$test[[i]], checks = results$nb[[i]], values = NA_real_,
      failed = results$failed[[i]], incomplete = results$incomplete[[i]],
      run = results$run[[i]], code = results$code[[i]],
      current = results$current[[i]], dataset = NA_character_, stringsAsFactors = FALSE)
    .matrix_expand(base, rules$fn[[hit]], rules$reference[[hit]], rules$kind[[hit]])
  })
  unrun <- on_disk[!paste(on_disk$dir, on_disk$file) %in% paste(results$dir, results$file), ]
  registered <- lapply(seq_len(nrow(unrun)), function(i) {
    own <- which(rules$dir == unrun$dir[[i]] & rules$file == unrun$file[[i]])
    base <- data.frame(
      evidence_id = paste0(unrun$dir[[i]], "/", unrun$file[[i]]), origin = "testthat",
      location = file.path(unrun$dir[[i]], unrun$file[[i]]), item = "(no saved result)",
      checks = NA_integer_, values = NA_real_, failed = NA_integer_, incomplete = 1L,
      run = "no saved result", code = NA_character_, current = FALSE,
      dataset = NA_character_, stringsAsFactors = FALSE)
    do.call(rbind, lapply(own, function(r) {
      .matrix_expand(base, rules$fn[[r]], rules$reference[[r]], rules$kind[[r]])
    }))
  })
  out <- do.call(rbind, c(saved, registered))
  out[out$fn %in% names(.matrix_functions), ]
}

#' One evidence row per quantity of a comparison kit (Latent GOLD kits, the
#' lca() package comparison). Agreeing quantities are checks; disagreements a
#' kit's register explains are known differences; "not comparable" rows are
#' left out.
#' @noRd
.matrix_kits <- function() {
  kits <- .matrix_rules("kit-evidence.csv")
  missing <- !file.exists(kits$path)
  if (any(missing)) message("kit evidence not yet produced, left out: ", paste(kits$path[missing], collapse = ", "))
  kits <- kits[!missing, ]
  functions <- .matrix_rules("kit-functions.csv")
  do.call(rbind, lapply(seq_len(nrow(kits)), function(i) {
    kit <- kits[i, ]
    rows <- utils::read.csv(kit$path, stringsAsFactors = FALSE)
    standard <- switch(kit$format,
      status = data.frame(case = rows$case, quantity = rows$quantity,
                          reference = rows$reference %||% rep(kit$reference, nrow(rows)),
                          fn = rows$fn, agrees = rows$status == "agrees",
                          known = rows$status == "known difference", values = rows$n_values),
      agrees = data.frame(case = rows$case, quantity = rows$quantity, reference = kit$reference,
                          fn = NA_character_, agrees = rows$agrees, known = FALSE, values = 1),
      latentgold = {
        rows <- rows[rows$status != "not comparable", ]
        data.frame(case = rows$case, quantity = rows$quantity, reference = kit$reference,
                   fn = NA_character_, agrees = rows$status == "agree", known = FALSE, values = 1)
      },
      ordinal = data.frame(case = rep(rows$case, 2L),
                           quantity = rep(c("log likelihood", "every parameter"), each = nrow(rows)),
                           reference = kit$reference, fn = NA_character_,
                           agrees = c(rows$ll_difference <= 6e-5,
                                      rows$max_parameter_difference <= 1e-4),
                           known = FALSE, values = c(rep(1, nrow(rows)), rows$n_parameters)))
    if (anyNA(standard$fn)) {
      rules <- functions[functions$path == kit$path, ]
      hit <- .matrix_first_match(standard$case, rules$case_pattern)
      stopifnot("every kit case is attributed to a function" = !anyNA(hit))
      standard$fn <- rules$fn[hit]
    }
    data.frame(
      evidence_id = paste0("kit:", kit$path, ":", seq_len(nrow(standard))), origin = "kit",
      location = kit$path, item = paste(standard$case, standard$quantity), fn = standard$fn,
      reference = ifelse(nzchar(standard$reference), standard$reference, kit$reference),
      kind = rows$kind %||% rep(kit$kind, nrow(standard)), checks = as.integer(standard$agrees),
      values = ifelse(standard$agrees, standard$values, 0), failed = as.integer(!standard$agrees & !standard$known),
      incomplete = 0L,
      run = kit$run, code = "0.98.1", current = matrix_evidence_current(kit$path),
      shared = FALSE, known = standard$known,
      dataset = .matrix_kit_datasets(kit$path, standard$case), stringsAsFactors = FALSE)
  }))
}

#' The dataset of each case of a kit, by the first rule of kit-datasets.csv
#' whose pattern matches
#' @noRd
.matrix_kit_datasets <- function(path, cases) {
  rules <- .matrix_rules("kit-datasets.csv")
  rules <- rules[rules$path == path, ]
  hit <- .matrix_first_match(cases, rules$pattern)
  stopifnot("every kit case maps to a dataset (kit-datasets.csv)" = !anyNA(hit))
  vapply(seq_along(cases), function(i) {
    sub(rules$pattern[[hit[[i]]]], rules$replacement[[hit[[i]]]], cases[[i]], perl = TRUE)
  }, character(1))
}

#' Every dataset compared with another program, package or paper, one row
#' per dataset, function and reference
#'
#' Current checks on mixtures (not one-class reductions) against external or
#' published references only. Harness and kit rows name their dataset; a
#' testthat file's datasets come from test-datasets.csv (its simulation kit),
#' linked to every reference the file's current mixture checks compare with.
#' Names are matched case-insensitively and without a trailing description
#' in parentheses ("srl (10% MCAR, mice pmm m = 5)" is srl), so a dataset
#' counts once however many kits use it.
#' @noRd
.matrix_dataset_links <- function(evidence) {
  counted <- evidence[evidence$current & !evidence$one_class & evidence$checks > 0 &
                        evidence$failed == 0L & evidence$incomplete == 0L &
                        evidence$kind %in% c(.matrix_external_kinds, "published"), ]
  named <- counted[!is.na(counted$dataset), c("fn", "reference", "kind", "dataset")]
  tests <- unique(counted[counted$origin == "testthat", c("location", "fn", "reference", "kind")])
  listed <- .matrix_rules("test-datasets.csv")
  linked <- merge(tests, listed, by = c("location", "fn"))
  out <- unique(rbind(named, linked[c("fn", "reference", "kind", "dataset")]))
  out$key <- tolower(sub(" \\([^()]*\\)$", "", out$dataset))
  row.names(out) <- NULL
  out
}

#' Every check, one row per function and reference it bears on
#' @noRd
.matrix_evidence <- function(files = .matrix_evidence_files()) {
  harness <- transform(.matrix_harness(files), shared = FALSE, known = FALSE)
  tests <- transform(.matrix_tests(files), known = FALSE)
  wrappers <- .matrix_rules("kit-evidence.csv")$test_file
  tests <- tests[!tests$location %in% wrappers, ]
  out <- rbind(harness, tests, .matrix_kits())
  stopifnot("every evidence row names a known function" =
              all(out$fn %in% names(.matrix_functions)))
  out$one_class <- grepl(.matrix_one_class_pattern, out$item, perl = TRUE)
  row.names(out) <- NULL
  out
}

#' Sum checks once per evidence item
#' @noRd
.matrix_unique_checks <- function(rows) {
  unique_rows <- rows[!duplicated(rows$evidence_id), ]
  c(checks = sum(unique_rows$checks, na.rm = TRUE),
    values = sum(unique_rows$values, na.rm = TRUE),
    failed = sum(unique_rows$failed, na.rm = TRUE),
    incomplete = sum(unique_rows$incomplete, na.rm = TRUE))
}

#' One row per function and reference: what checks it, and how current
#' @noRd
.matrix_cells <- function(evidence, candidates, links = .matrix_dataset_links(evidence)) {
  tested <- evidence[!evidence$kind %in% c("provenance", "infrastructure"), ]
  keys <- unique(tested[c("fn", "reference", "kind")])
  cells <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    rows <- tested[tested$fn == keys$fn[[i]] & tested$reference == keys$reference[[i]] &
                     tested$kind == keys$kind[[i]], ]
    now <- .matrix_unique_checks(rows[rows$current, ])
    mixture <- .matrix_unique_checks(rows[rows$current & !rows$one_class, ])
    old <- .matrix_unique_checks(rows[!rows$current & rows$run != "no saved result", ])
    data.frame(keys[i, ], current_checks = now[["checks"]], current_values = now[["values"]],
               current_failed = now[["failed"]], current_incomplete = now[["incomplete"]],
               mixture_checks = mixture[["checks"]],
               older_checks = old[["checks"]],
               older_run = paste(unique(rows$code[!rows$current & !is.na(rows$code)]),
                                 collapse = ", "),
               unrun_files = sum(rows$run == "no saved result"),
               n_datasets = length(unique(links$key[links$fn == keys$fn[[i]] &
                                                      links$reference == keys$reference[[i]]])),
               shared = any(rows$shared), candidate = FALSE, stringsAsFactors = FALSE)
  }))
  untested <- candidates[!paste(candidates$fn, candidates$reference) %in%
                           paste(cells$fn, cells$reference), ]
  candidate_cells <- data.frame(
    fn = untested$fn, reference = untested$reference, kind = untested$kind,
    current_checks = rep(0L, nrow(untested)), current_values = rep(0, nrow(untested)),
    current_failed = rep(0L, nrow(untested)), current_incomplete = rep(0L, nrow(untested)),
    mixture_checks = rep(0L, nrow(untested)), older_checks = rep(0L, nrow(untested)),
    older_run = rep("", nrow(untested)), unrun_files = rep(0L, nrow(untested)),
    n_datasets = rep(0L, nrow(untested)), shared = rep(FALSE, nrow(untested)),
    candidate = rep(TRUE, nrow(untested)),
    stringsAsFactors = FALSE)
  out <- rbind(cells, candidate_cells)
  out$status <- ifelse(out$candidate, "candidate",
    ifelse(out$mixture_checks > 0L, "current",
      ifelse(out$current_checks > 0L, "one class only",
        ifelse(out$older_checks > 0L, "older run", "no saved result"))))
  out$status[!out$candidate & out$current_incomplete > 0L] <- "incomplete"
  out$status[!out$candidate & out$current_failed > 0L] <- "failed"
  out <- out[order(match(out$fn, names(.matrix_functions)),
                   match(out$kind, .matrix_kind_order), -out$current_checks,
                   out$reference), ]
  row.names(out) <- NULL
  out
}

#' Equivalence strength from the number of sources and published papers
#'
#' A source is an independent external program or package that agrees on a
#' model with two or more classes, with current results. One + per source and
#' one P per published paper; the colour (`tier`) marks none, one, two, or
#' three or more sources.
#' @param sources,papers Counts.
#' @return A list: the `tier` (a CSS class), the `plus` and `papers` marks,
#'   and the joined `label` ("0" without a source).
#' @noRd
.matrix_strength <- function(sources, papers) {
  tier <- if (sources >= 3L) "very_strong" else if (sources == 2L) "strong" else
    if (sources == 1L) "single" else "none"
  plus <- strrep("+", sources)
  marks <- strrep("P", papers)
  label <- paste(c(if (sources) plus else "0", if (papers) marks), collapse = " ")
  list(tier = tier, plus = plus, papers = marks, label = label)
}

#' One row per function: its equivalence strength and how to raise it
#' @noRd
.matrix_summary <- function(evidence, cells, links = .matrix_dataset_links(evidence)) {
  do.call(rbind, lapply(names(.matrix_functions), function(fn) {
    own <- cells[cells$fn == fn, ]
    rows <- evidence[evidence$fn == fn & !evidence$kind %in% c("provenance", "infrastructure"), ]
    external <- own$reference[own$kind %in% .matrix_external_kinds & own$status == "current"]
    published <- own$reference[own$kind == "published" & own$status == "current"]
    older_only <- own$reference[own$status == "older run" &
                                  own$kind %in% c(.matrix_external_kinds, "published")]
    one_class_only <- own$reference[own$status == "one class only" &
                                      own$kind %in% c(.matrix_external_kinds, "published")]
    candidate_papers <- own$reference[own$candidate & own$kind == "published"]
    candidate_sources <- own$reference[own$candidate & own$kind %in% .matrix_external_kinds]
    unrun <- length(unique(rows$location[rows$run == "no saved result"]))
    now <- .matrix_unique_checks(rows[rows$current & rows$kind != "internal", ])
    strength <- .matrix_strength(length(external), length(published))
    raise <- c(
      if (now[["failed"]] > 0) sprintf("%d failing checks", now[["failed"]]),
      if (now[["incomplete"]] > 0) sprintf("%d incomplete tests", now[["incomplete"]]),
      if (length(one_class_only)) sprintf("one-class checks only, not counted: %s",
                                          paste(one_class_only, collapse = ", ")),
      if (length(older_only)) sprintf("rerun stale evidence: adds %s",
                                      paste(older_only, collapse = ", ")),
      if (unrun) sprintf("%d test file(s) never saved a result", unrun),
      if (!length(published) && length(candidate_papers))
        sprintf("first paper from: %s", paste(candidate_papers, collapse = ", ")),
      if (!length(published) && !length(candidate_papers)) "no published paper found",
      if (length(external) < 3L && length(candidate_sources))
        sprintf("more sources from: %s", paste(candidate_sources, collapse = ", ")))
    data.frame(
      fn = fn, label = .matrix_functions[[fn]], strength = strength$label,
      tier = strength$tier, strength_plus = strength$plus,
      strength_papers = strength$papers, n_sources = length(external),
      n_datasets = length(unique(links$key[links$fn == fn])),
      sources = paste(external, collapse = ", "), n_papers = length(published),
      papers = paste(published, collapse = ", "),
      checks_current = now[["checks"]], values_current = now[["values"]],
      failed_current = now[["failed"]], incomplete_current = now[["incomplete"]],
      internal_checks = .matrix_unique_checks(rows[rows$current & rows$kind == "internal", ])[["checks"]],
      candidates = sum(own$candidate),
      to_strengthen = if (length(raise)) paste(raise, collapse = "; ") else "",
      stringsAsFactors = FALSE)
  }))
}

#' A strength badge
#' @noRd
.matrix_badge <- function(tier, label) {
  sprintf('<span class="badge %s">%s</span>', tier, .matrix_escape(label))
}

#' Escape text for HTML
#' @noRd
.matrix_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' A data.frame as an HTML table; `raw` names columns already holding HTML
#' @noRd
.matrix_html_table <- function(df, class = "", raw = character()) {
  header <- paste0("<th>", .matrix_escape(names(df)), "</th>", collapse = "")
  body <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(names(df), function(column) {
      value <- as.character(df[[column]][[i]])
      value <- if (is.na(value)) "" else value
      paste0("<td>", if (column %in% raw) value else .matrix_escape(value), "</td>")
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  }, character(1))
  sprintf('<table class="%s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>',
          class, header, paste(body, collapse = "\n"))
}

#' The text of one matrix cell
#' @noRd
.matrix_cell_label <- function(cell) {
  if (!nrow(cell)) return("")
  if (isTRUE(cell$candidate)) return('<span class="cand">candidate</span>')
  star <- if (isTRUE(cell$shared)) "*" else ""
  if (cell$status %in% c("failed", "incomplete")) return(cell$status)
  if (cell$mixture_checks > 0L) {
    return(sprintf('<span class="cur">%s%s</span>', format(cell$current_checks, big.mark = ","),
                   star))
  }
  if (cell$current_checks > 0L) {
    return(sprintf('<span class="one">%s%s (1 class)</span>',
                   format(cell$current_checks, big.mark = ","), star))
  }
  if (cell$older_checks > 0L) {
    return(sprintf('<span class="old">%s (%s)</span>',
                   format(cell$older_checks, big.mark = ","), cell$older_run))
  }
  '<span class="old">no result</span>'
}

#' Column headings for the grid: one per function, wrapped to fit
#' @noRd
.matrix_short_labels <- c(
  lpa = "lpa()", multilpa = "multilpa()", lca = "lca()", multilca = "multilca()",
  lta = "lta()", mixture_regression = "mixture_<br>regression()",
  growth = "growth<br>(random =)", three_step = "three_step()<br>r3step()",
  bootstrap_lrt = "bootstrap_<br>lrt()<br>VLMR", enumeration = "enumerate_*()<br>compare_<br>models()",
  pool_imputations = "pool_<br>imputations()")

#' The reference by function grid: one row per reference, grouped by kind
#' @noRd
.matrix_grid <- function(cells, summary) {
  shown <- cells[cells$kind %in% c(.matrix_external_kinds, "published", "component check", "independent code"), ]
  weight <- tapply(shown$current_checks + shown$older_checks, shown$reference, sum)
  references <- unique(shown[c("reference", "kind")])
  references <- references[order(match(references$kind, .matrix_kind_order),
                                 -weight[references$reference], references$reference), ]
  row_html <- function(reference, kind) {
    labels <- vapply(names(.matrix_functions), function(fn) {
      .matrix_cell_label(shown[shown$fn == fn & shown$reference == reference & shown$kind == kind, ])
    }, character(1))
    paste0("<tr><td class=\"ref\">", .matrix_escape(reference), "</td>",
           paste0("<td>", labels, "</td>", collapse = ""), "</tr>")
  }
  groups <- vapply(unique(references$kind), function(kind) {
    own <- references$reference[references$kind == kind]
    paste0(sprintf('<tr class="group"><td colspan="%d">%s</td></tr>',
                   length(.matrix_functions) + 1L, .matrix_escape(kind)),
           paste(vapply(own, row_html, character(1), kind = kind), collapse = ""))
  }, character(1))
  # Narrow grid cells: the + marks and the P marks on separate lines.
  stacked <- vapply(seq_len(nrow(summary)), function(i) {
    lines <- c(if (summary$n_sources[[i]]) summary$strength_plus[[i]] else "0",
               summary$strength_papers[[i]])
    sprintf('<span class="badge %s">%s</span>', summary$tier[[i]],
            paste(.matrix_escape(lines[nzchar(lines)]), collapse = "<br>"))
  }, character(1))
  status_row <- paste0("<tr class=\"status\"><td class=\"ref\">strength (+ source, P paper)</td>",
                       paste0("<td>", stacked, "</td>", collapse = ""), "</tr>")
  dataset_row <- paste0("<tr class=\"total\"><td class=\"ref\">datasets tested</td>",
                        paste0("<td>", format(summary$n_datasets, big.mark = ","), "</td>",
                               collapse = ""), "</tr>")
  total_row <- paste0("<tr class=\"total\"><td class=\"ref\">current external/component/oracle checks</td>",
                      paste0("<td>", format(summary$checks_current, big.mark = ","), "</td>",
                             collapse = ""), "</tr>")
  header <- paste0("<th>reference</th>",
                   paste0("<th>", .matrix_short_labels[names(.matrix_functions)], "</th>",
                          collapse = ""))
  sprintf('<table class="grid"><thead><tr>%s</tr></thead><tbody>%s%s%s%s</tbody></table>',
          header, status_row, dataset_row, paste(groups, collapse = ""), total_row)
}

#' The whole report as one HTML page
#' @noRd
.matrix_html <- function(evidence, cells, summary, candidates) {
  version <- read.dcf("DESCRIPTION", fields = "Version")[[1L]]
  files <- .matrix_evidence_files()
  totals <- .matrix_unique_checks(evidence[evidence$current &
                                             !evidence$kind %in% c("provenance", "infrastructure", "internal"), ])
  overview <- data.frame(
    `function` = summary$label,
    strength = .matrix_badge(summary$tier, summary$strength),
    sources = summary$n_sources,
    papers = summary$n_papers,
    datasets = format(summary$n_datasets, big.mark = ","),
    `checks passed` = format(summary$checks_current, big.mark = ","),
    failed = summary$failed_current, incomplete = summary$incomplete_current,
    `values (harness)` = ifelse(summary$values_current > 0,
                                format(summary$values_current, big.mark = ","), ""),
    `independent programs and packages` = summary$sources,
    `published papers` = summary$papers,
    `to strengthen` = summary$to_strengthen,
    check.names = FALSE, stringsAsFactors = FALSE)
  detail <- do.call(paste0, lapply(names(.matrix_functions), function(fn) {
    own <- cells[cells$fn == fn & !cells$candidate, ]
    table <- data.frame(
      reference = own$reference, kind = own$kind, status = own$status,
      `current checks` = format(own$current_checks, big.mark = ","),
      datasets = ifelse(own$n_datasets > 0L, own$n_datasets, ""),
      `values` = ifelse(own$current_values > 0, format(own$current_values, big.mark = ","), ""),
      `stale checks` = ifelse(own$older_checks > 0L, format(own$older_checks, big.mark = ","), ""),
      `files without result` = ifelse(own$unrun_files > 0L, own$unrun_files, ""),
      shared = ifelse(own$shared, "*", ""), check.names = FALSE, stringsAsFactors = FALSE)
    wanted <- candidates[candidates$fn == fn, ]
    candidate_table <- if (nrow(wanted)) paste0(
      "<h3>Candidates not yet wired in</h3>",
      .matrix_html_table(data.frame(reference = wanted$reference, kind = wanted$kind,
                                    data = wanted$data, `published numbers` = wanted$published,
                                    checked = wanted$mark, check.names = FALSE),
                         "detail")) else ""
    sprintf('<section class="fn"><h2>%s</h2>%s%s</section>',
            .matrix_escape(.matrix_functions[[fn]]),
            if (nrow(table)) .matrix_html_table(table, "detail") else "<p>No evidence.</p>",
            candidate_table)
  }))
  sources <- data.frame(input = files$path, run = files$run, `package version` = files$code,
                        used = ifelse(files$current, "current results",
                                      "only for files without a current result"),
                        check.names = FALSE)
  css <- paste(
    "@page { size: A4 landscape; margin: 12mm; }",
    ":root { --vstrong: #009E73; --strong: #56B4E9; --single: #E69F00; --none: #999999; --ink: #1d1d1f; --line: #d0d0d0; }",
    "body { font-family: -apple-system, 'Helvetica Neue', Arial, sans-serif; color: var(--ink); background: #ffffff; font-size: 9pt; }",
    "h1 { font-size: 17pt; margin: 0 0 4px; } h2 { font-size: 12pt; margin: 14px 0 4px; } h3 { font-size: 10pt; margin: 8px 0 3px; }",
    "p { margin: 3px 0; line-height: 1.35; } .meta { color: #555; }",
    "table { border-collapse: collapse; width: 100%; margin: 4px 0 8px; }",
    "th, td { border: 1px solid var(--line); padding: 2px 4px; vertical-align: top; text-align: left; }",
    "th { background: #f2f2f2; font-weight: 600; }",
    "table.grid { font-size: 7.5pt; table-layout: fixed; } table.grid th { text-align: center; }",
    "table.grid th:first-child, table.grid td.ref { width: 24%; text-align: left; }",
    "table.grid td { text-align: center; padding: 1px 3px; } table.grid tr { break-inside: avoid; }",
    "table.grid tr.group td { background: #f7f7f7; font-weight: 600; text-align: left; text-transform: uppercase; font-size: 6.5pt; letter-spacing: .04em; }",
    "table.grid tr.total td { border-top: 2px solid #888; font-weight: 600; }",
    ".badge { display: inline-block; padding: 0 5px; border-radius: 3px; font-weight: 600; white-space: nowrap; }",
    ".badge.very_strong { background: var(--vstrong); color: #fff; } .badge.strong { background: var(--strong); color: #000; }",
    ".badge.single { background: var(--single); color: #000; } .badge.none { background: var(--none); color: #fff; }",
    "table.grid tr.status td { font-size: 6.5pt; }",
    ".cur { font-weight: 600; } .old { color: #777; font-style: italic; } .cand { color: #0072B2; }",
    ".one { color: #D55E00; }",
    "section.fn { break-inside: avoid; } .pagebreak { break-before: page; }",
    "table.detail { font-size: 8pt; } table.overview { font-size: 7.5pt; }",
    "table.overview td, table.overview th { padding: 1px 4px; line-height: 1.25; }",
    "table.overview td:nth-child(3), table.overview td:nth-child(4), table.overview td:nth-child(5), table.overview td:nth-child(6), table.overview td:nth-child(7), table.overview td:nth-child(8) { text-align: right; }",
    ".legend { border: 1px solid var(--line); border-radius: 4px; padding: 5px 8px; margin: 6px 0 8px; background: #fafafa; }",
    ".legend h3 { margin: 0 0 3px; } .legend dl { display: grid; grid-template-columns: max-content 1fr; gap: 2px 10px; margin: 0; }",
    ".legend dt { font-weight: 600; } .legend dd { margin: 0; }", sep = "\n")
  paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>latents verification matrix</title>",
    "<style>", css, "</style></head><body>",
    "<h1>latents verification matrix</h1>",
    sprintf('<p class="meta">Package %s. Built from saved evidence on %s; nothing was refitted. %s current checks (%s values) against references other than latents itself.</p>',
            version, format(Sys.Date()), format(totals[["checks"]], big.mark = ","),
            format(totals[["values"]], big.mark = ",")),
    "<p>A <b>check</b> is one harness quantity (a whole output column counts once; <i>values</i> counts every number inside it) or one testthat expectation. ",
    "<b>Current</b> requires matching runtime-source and result-file hashes in an evidence sidecar. Historical 1.0.0 evidence carries an explicit retrospective source-equivalence attestation, not a claimed rerun. Missing or mismatched provenance is stale and does not count toward strength. ",
    "</p>",
    "<div class=\"legend\"><h3>Legend</h3><dl>",
    "<dt>Strength</dt><dd>One <b>+</b> per <i>source</i>: an independent program or package (another R package, Mplus, Latent GOLD, Stata, Python) whose results agree on the current code for a model with two or more classes. ",
    "One <b>P</b> per published paper whose printed numbers are reproduced. <span class=\"badge very_strong\">++++ PP</span> is four sources and two papers. ",
    "Colour: <span class=\"badge none\">0</span> no source, <span class=\"badge single\">+</span> one, <span class=\"badge strong\">++</span> two, <span class=\"badge very_strong\">+++</span> three or more.</dd>",
    "<dt>Not counted</dt><dd>Agreement on a one-class reduction only (one class equals lm(), glm.nb(), polr(), lme4, a linear mixed model), shown as \"(1 class)\"; component checks (formulas and helper routines); our own independent code (numDeriv, dense likelihoods, the kits' reference implementations); internal checks; stale results; sources with failed or incomplete tests.</dd>",
    "<dt>Sources, papers</dt><dd>The number of + and of P.</dd>",
    "<dt>Datasets</dt><dd>Distinct datasets on which the function was compared with a source or a paper on current code, for models with two or more classes: the harness's datasets, the kits' datasets (kit-datasets.csv) and the simulated datasets of the testthat kits that the other program ran (test-datasets.csv). A dataset used by several kits counts once. Tests that build their own small data inline are not counted, so this is a lower bound.</dd>",
    "<dt>Checks passed</dt><dd>Current checks against every reference other than latents itself: one harness quantity (a whole output column), one kit quantity or one testthat expectation. <i>Values</i> counts the numbers inside the harness checks.</dd>",
    "<dt>Failed</dt><dd>Current disagreements without a registered explanation, plus test errors. Skipped or incomplete tests provide no source credit.</dd>",
    "<dt>Grid cells</dt><dd><span class=\"cur\">bold</span>: current checks on mixtures; <span class=\"one\">(1 class)</span>: one-class reductions only; <span class=\"old\">grey italic</span>: stale provenance or no saved result; <span class=\"cand\">candidate</span>: found, not yet wired in; *: a test that checks several functions or references at once (counted once in each function's total).</dd>",
    "</dl></div>",
    "<h2>Strength by function</h2>", .matrix_html_table(overview, "overview", raw = "strength"),
    "<div class=\"pagebreak\"></div><h2>Reference by function</h2>",
    "<p>Each cell: the current checks of a function against a reference (styles in the legend). Internal checks (invariance, seed stability, refusal parity, identities) are in the evidence pages, not in this grid.</p>",
    .matrix_grid(cells, summary),
    "<div class=\"pagebreak\"></div><h2>Evidence by function</h2>", detail,
    "<div class=\"pagebreak\"></div><h2>Inputs</h2>", .matrix_html_table(sources, "detail"),
    "<p>Attribution rules: harness-functions.csv, harness-references.csv and test-attribution.csv beside build.R. ",
    "lca(): every indicator categorical, one level; multilca(): the same with groups; lpa(): every indicator Gaussian, one level; multilpa(): everything else, including count and ordinal indicators and tests that call the engine directly.</p>",
    "</body></html>")
}

#' A data.frame as a Markdown table (pipes in values escaped)
#' @noRd
.matrix_md_table <- function(df) {
  cell <- function(x) gsub("|", "\\|", ifelse(is.na(x), "", as.character(x)), fixed = TRUE)
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  rule <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- vapply(seq_len(nrow(df)), function(i) {
    paste0("| ", paste(vapply(df, function(column) cell(column[[i]]), character(1)),
                       collapse = " | "), " |")
  }, character(1))
  c(header, rule, body)
}

#' The matrix as a Markdown document that GitHub renders: the summary, the
#' legend, and per function the references with current checks
#' @noRd
.matrix_markdown <- function(evidence, cells, summary) {
  version <- read.dcf("DESCRIPTION", fields = "Version")[[1L]]
  totals <- .matrix_unique_checks(evidence[evidence$current &
                                             !evidence$kind %in% c("provenance", "infrastructure", "internal"), ])
  overview <- data.frame(
    `function` = paste0("`", summary$label, "`"),
    strength = paste0("`", summary$strength, "`"),
    sources = summary$n_sources, papers = summary$n_papers,
    datasets = trimws(format(summary$n_datasets, big.mark = ",")),
    checks = trimws(format(summary$checks_current, big.mark = ",")),
    failed = summary$failed_current, incomplete = summary$incomplete_current, check.names = FALSE, stringsAsFactors = FALSE)
  per_function <- unlist(lapply(seq_len(nrow(summary)), function(i) {
    fn <- summary$fn[[i]]
    own <- cells[cells$fn == fn & cells$status %in% c("current", "one class only", "failed", "incomplete"), ]
    table <- data.frame(
      reference = own$reference, kind = own$kind, status = own$status,
      checks = paste0(trimws(format(own$current_checks, big.mark = ",")),
                      ifelse(own$status == "one class only", " (1 class)", "")),
      datasets = ifelse(own$n_datasets > 0L, own$n_datasets, ""),
      stringsAsFactors = FALSE)
    c(sprintf("### `%s` — `%s`", summary$label[[i]], summary$strength[[i]]), "",
      sprintf("- Sources: %s.", if (nzchar(summary$sources[[i]])) summary$sources[[i]] else "none"),
      sprintf("- Published papers: %s.", if (nzchar(summary$papers[[i]])) summary$papers[[i]] else "none"),
      if (nzchar(summary$to_strengthen[[i]])) sprintf("- To strengthen: %s.", summary$to_strengthen[[i]]),
      "", .matrix_md_table(table), "")
  }))
  c("# latents verification matrix", "",
    sprintf(paste("Package %s, built %s from saved evidence (nothing refitted): %s current checks",
                  "(%s values) against references other than latents itself. Regenerated by",
                  "`Rscript equivalence/verification-matrix/build.R`; do not edit by hand.",
                  "Full detail, with the reference-by-function grid and older results:",
                  "[VERIFICATION-MATRIX.html](VERIFICATION-MATRIX.html) /",
                  "[.pdf](VERIFICATION-MATRIX.pdf)."),
            version, format(Sys.Date()), format(totals[["checks"]], big.mark = ","),
            format(totals[["values"]], big.mark = ",")), "",
    "## Strength by function", "", .matrix_md_table(overview), "",
    "## Legend", "",
    paste("- **Strength**: one `+` per *source*, an independent program or package (another R",
          "package, Mplus, Latent GOLD, Stata, Python) whose results agree on the current code",
          "for a model with two or more classes; one `P` per published paper whose printed",
          "numbers are reproduced. `++++ PP` is four sources and two papers."),
    paste("- **Not counted**: agreement on one-class reductions only (marked \"(1 class)\"), our",
          "own independent code, component checks, internal checks, stale results, and sources with failed or incomplete tests."),
    paste("- **Datasets**: distinct datasets on which the function was compared with a source or",
          "a paper on current code, for models with two or more classes (the harness, the kits",
          "and the testthat simulation kits that the other program ran). A dataset used by",
          "several kits counts once; tests with small inline data are not counted, so this is a",
          "lower bound."),
    paste("- **Current**: both the complete runtime-source manifest and the result-file hash match their evidence sidecar. Historical 1.0.0 evidence uses an explicitly labelled retrospective source-equivalence attestation; it is not a new run."),
    paste("- **Checks**: current harness quantities, kit quantities and passed testthat expectations",
          "against every reference other than latents itself. **Failed**: checks that disagree",
          "without a registered, explained difference, plus test errors. Skips are incomplete and earn no source credit."),
    "", "## By function", "", per_function)
}

#' Print an HTML file to PDF with headless Chrome
#' @noRd
.matrix_pdf <- function(html, pdf) {
  chrome <- Sys.getenv("LATENTS_CHROME", file.path(
    .Platform$file.sep, "Applications", "Google Chrome.app", "Contents", "MacOS", "Google Chrome"))
  if (!file.exists(chrome)) {
    warning(warningCondition(sprintf("Chrome not found at %s; only the HTML was written.", chrome),
                             class = "latents_matrix_no_chrome"))
    return(invisible(FALSE))
  }
  status <- system2(chrome, c("--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                              shQuote(paste0("--print-to-pdf=", pdf)),
                              shQuote(paste0("file://", normalizePath(html)))),
                    stdout = FALSE, stderr = FALSE)
  stopifnot("Chrome printed the PDF" = identical(status, 0L) && file.exists(pdf))
  invisible(TRUE)
}

#' Build every output
#' @noRd
build_verification_matrix <- function() {
  evidence <- .matrix_evidence()
  candidates <- .matrix_rules("candidates.csv")
  cells <- .matrix_cells(evidence, candidates)
  summary <- .matrix_summary(evidence, cells)
  utils::write.csv(evidence, file.path(.matrix_dir, "evidence.csv"), row.names = FALSE)
  utils::write.csv(cells, file.path(.matrix_dir, "matrix.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(.matrix_dir, "summary.csv"), row.names = FALSE)
  html <- file.path(.matrix_dir, "VERIFICATION-MATRIX.html")
  writeLines(.matrix_html(evidence, cells, summary, candidates), html, useBytes = TRUE)
  writeLines(.matrix_markdown(evidence, cells, summary),
             file.path(.matrix_dir, "VERIFICATION-MATRIX.md"), useBytes = TRUE)
  .matrix_pdf(html, file.path(.matrix_dir, "VERIFICATION-MATRIX.pdf"))
  invisible(summary[c("label", "strength", "n_sources", "n_papers", "n_datasets", "checks_current",
                      "failed_current")])
}

if (sys.nframe() == 0L) print(build_verification_matrix(), right = FALSE)
