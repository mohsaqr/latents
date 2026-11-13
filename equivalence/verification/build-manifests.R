# CVS verification manifests: one YAML file per verified function in
# modules/<fn>.yaml and the papers they cite in references.bib, for the
# per-target Verification tab of CVS (../validation/docs/TARGET-VERIFICATION.md).
#
#   Rscript equivalence/verification/build-manifests.R      (from the project root)
#
# Generated from the verification matrix (../verification-matrix/: summary.csv,
# evidence.csv, matrix.csv and its rule tables) and the help pages; never edit
# the output by hand. Rebuild after the matrix is rebuilt. The output is sorted
# and has no time stamp, so an unchanged matrix rebuilds byte for byte.
#
# Which report entries a manifest claims (`cvs_files`, `cvs_test_filter`) is
# decided with the CVS reporter's own tagging (../cvs/cvs-report.R): a block
# whose references fall on several CVS axes is never tagged, and CVS would
# count it as R, so no manifest claims it.

source(file.path("equivalence", "verification-matrix", "build.R"), local = TRUE)
source(file.path("equivalence", "cvs", "cvs-report.R"), local = TRUE)

.manifest_dir <- file.path("equivalence", "verification")

# CVS family of each function (FAMILY_ORDER in ../validation's verification tab).
.manifest_family <- c(
  lpa = "cluster", lca = "cluster", multilpa = "cluster", multilca = "cluster",
  enumeration = "cluster", bootstrap_lrt = "cluster", pool_imputations = "cluster",
  lta = "dynamics", mixture_regression = "regression", growth = "regression",
  three_step = "regression")

# Kinds of reference whose tolerances a manifest lists (internal checks excluded).
.manifest_reference_kinds <- c("R package", "commercial program", "published", "independent code")

#' Read one output table of the matrix
#' @noRd
.manifest_matrix_table <- function(name) {
  utils::read.csv(file.path(.matrix_dir, name), stringsAsFactors = FALSE)
}

#' Exported functions a matrix label names: "three_step() / r3step()" gives
#' both; "enumerate_*()" every export it matches
#' @noRd
.manifest_public_api <- function(label, exports) {
  names <- regmatches(label, gregexpr("[A-Za-z_.0-9*]+(?=\\()", label, perl = TRUE))[[1L]]
  api <- unlist(lapply(names, function(name) {
    if (grepl("*", name, fixed = TRUE)) {
      sort(grep(glob2rx(name), exports, value = TRUE), method = "radix")
    } else name
  }))
  stopifnot("every named function is exported" = length(api) > 0L && all(api %in% exports))
  unique(api)
}

#' The R/ file defining each function
#' @noRd
.manifest_module_files <- function(api) {
  files <- sort(list.files("R", "[.][Rr]$", full.names = TRUE), method = "radix")
  definitions <- lapply(files, function(path) {
    sub("^([A-Za-z_.0-9]+) *<- *function.*$", "\\1",
        grep("^[A-Za-z_.0-9]+ *<- *function", readLines(path, warn = FALSE), value = TRUE))
  })
  defined_in <- vapply(api, function(name) {
    hit <- which(vapply(definitions, function(defined) name %in% defined, logical(1)))
    stopifnot("each public function is defined in exactly one R/ file" = length(hit) == 1L)
    files[[hit]]
  }, character(1))
  sort(unique(unname(defined_in)), method = "radix")
}

#' The \title of a help page
#' @noRd
.manifest_rd_title <- function(name) {
  lines <- readLines(file.path("man", paste0(name, ".Rd")), warn = FALSE)
  title <- sub("^\\\\title\\{(.*)\\}$", "\\1", grep("^\\\\title\\{.*\\}$", lines, value = TRUE))
  stopifnot("the help page has a one-line title" = length(title) == 1L)
  title
}

#' BibTeX key, author and year of a published reference
#'
#' "Masyn (2013), Oxford Handbook ..." is masyn2013, author Masyn, year 2013;
#' a reference with no "(year)" after its authors keys on its first two words.
#' @return A data.frame with `reference`, `key`, `author`, `year`.
#' @noRd
.manifest_bib_keys <- function(reference) {
  dated <- regmatches(reference, regexec("^(.*?)\\s*\\(((?:19|20)[0-9]{2})\\)", reference, perl = TRUE))
  ascii <- function(x) {
    tolower(gsub("[^A-Za-z0-9]", "", iconv(x, "UTF-8", "ASCII//TRANSLIT", sub = "")))
  }
  rows <- lapply(seq_along(reference), function(i) {
    match <- dated[[i]]
    if (length(match)) {
      author <- match[[2L]]
      first <- strsplit(author, "[ ,&]+")[[1L]][[1L]]
      data.frame(reference = reference[[i]], key = paste0(ascii(first), match[[3L]]),
                 author = author, year = match[[3L]], stringsAsFactors = FALSE)
    } else {
      words <- strsplit(reference[[i]], "\\s+")[[1L]]
      data.frame(reference = reference[[i]], key = ascii(paste(utils::head(words, 2L), collapse = "")),
                 author = NA_character_, year = NA_character_, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}

#' Tolerances each function is held to, per reference
#'
#' From the harness comparisons the matrix reads and the per-quantity kits
#' that name their function. A known difference (registered, explained) is
#' left out; any other quantity outside its tolerance makes the entry "fail".
#' @return A data.frame: `fn`, `function`, `reference`, `tolerance`, `status`.
#' @noRd
.manifest_tolerances <- function() {
  files <- .matrix_evidence_files()
  report <- utils::read.csv(files$path[files$role == "harness"], stringsAsFactors = FALSE)
  fn_rules <- .matrix_rules("harness-functions.csv")
  ref_rules <- .matrix_rules("harness-references.csv")
  ref_index <- matrix_harness_reference(report$source, ref_rules)
  harness <- data.frame(
    fn = fn_rules$fn[matrix_harness_function(report$suite, report$quantity, fn_rules)],
    reference = ref_rules$reference[ref_index], kind = ref_rules$kind[ref_index],
    source = "harness", tolerance = report$tolerance, agrees = as.logical(report$agrees),
    stringsAsFactors = FALSE)
  kits <- .matrix_rules("kit-evidence.csv")
  kits <- kits[kits$format == "status" & file.exists(kits$path), ]
  kit_rows <- lapply(seq_len(nrow(kits)), function(i) {
    rows <- utils::read.csv(kits$path[[i]], stringsAsFactors = FALSE)
    rows <- rows[rows$status %in% c("agrees", "disagrees", "error"), ]
    reference <- rows$reference %||% rep(kits$reference[[i]], nrow(rows))
    data.frame(fn = rows$fn,
               reference = ifelse(nzchar(reference), reference, kits$reference[[i]]),
               kind = rows$kind %||% rep(kits$kind[[i]], nrow(rows)),
               source = basename(dirname(kits$path[[i]])), tolerance = rows$tolerance,
               agrees = rows$status == "agrees", stringsAsFactors = FALSE)
  })
  all_rows <- do.call(rbind, c(list(harness), kit_rows))
  all_rows <- all_rows[all_rows$kind %in% .manifest_reference_kinds & !is.na(all_rows$tolerance) &
                         !is.na(all_rows$fn), ]
  groups <- unique(all_rows[c("fn", "reference", "source")])
  entries <- lapply(seq_len(nrow(groups)), function(i) {
    rows <- all_rows[all_rows$fn == groups$fn[[i]] & all_rows$reference == groups$reference[[i]] &
                       all_rows$source == groups$source[[i]], ]
    bounds <- signif(range(rows$tolerance), 2L)
    tolerance <- if (isTRUE(all.equal(bounds[[1L]], bounds[[2L]]))) format(bounds[[1L]]) else
      sprintf("%s to %s", format(bounds[[1L]]), format(bounds[[2L]]))
    data.frame(fn = groups$fn[[i]],
               `function` = sprintf("%s: %d quantities (%s)", .matrix_functions[[groups$fn[[i]]]],
                                    nrow(rows), groups$source[[i]]),
               reference = groups$reference[[i]],
               tolerance = sprintf("%s (%d of %d within)", tolerance, sum(rows$agrees), nrow(rows)),
               status = if (all(rows$agrees)) "pass" else "fail",
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, entries)
  out[order(out$fn, out$reference, out$`function`, method = "radix"), ]
}

#' Does CVS's `cvs_test_filter` needle claim a test? (matchesNeedle in
#' ../validation/src/dump-verification.ts: a trailing `$` anchors to the end of
#' the title; otherwise a case-insensitive substring of title and ancestors)
#' @noRd
.manifest_needle_claims <- function(title, ancestor, needle) {
  needle <- tolower(needle)
  if (endsWith(needle, "$")) {
    return(endsWith(tolower(title), substr(needle, 1L, nchar(needle) - 1L)))
  }
  grepl(needle, tolower(paste0(title, "\n", ancestor)), fixed = TRUE)
}

#' Needles that claim exactly the wanted titles of one file
#'
#' Tries each candidate needle (a shared prefix, the whole title, the whole
#' title anchored with `$`) and keeps the first that claims none of the
#' `unwanted` titles.
#' @param wanted,unwanted Titles in the file that must and must not be claimed.
#' @param prefixes Shorter candidate needles, one per wanted title (or NULL).
#' @param ancestor The ancestor title CVS also matches against.
#' @noRd
.manifest_needles <- function(wanted, unwanted, ancestor, prefixes = NULL) {
  candidates <- lapply(seq_along(wanted), function(i) {
    c(if (!is.null(prefixes)) prefixes[[i]], wanted[[i]], paste0(wanted[[i]], "$"))
  })
  claims_unwanted <- function(needle) {
    any(vapply(unwanted, .manifest_needle_claims, logical(1), ancestor = ancestor, needle = needle))
  }
  needles <- vapply(candidates, function(options) {
    usable <- options[!vapply(options, claims_unwanted, logical(1))]
    if (!length(usable)) {
      stop(errorCondition(sprintf("no needle claims \"%s\" alone", options[[2L]]),
                          class = "latents_manifest_ambiguous_needle", call = NULL))
    }
    usable[[1L]]
  }, character(1))
  sort(unique(needles), method = "radix")
}

#' The report entries (and the tests within them) each function claims
#'
#' Locations of evidence.csv become report entry names as the CVS reporter
#' writes them: harness/suites/<suite>.R is equivalence/harness/suites/<suite>.R,
#' a kit is the testthat file that wraps it (kit-evidence.csv), a test file
#' keeps its path. tests/equivalence/ is not run by any reporting runner.
#' A file is claimed whole when every test in it is this function's evidence
#' or an internal check; otherwise only this function's tests, by needle. A
#' test whose references fall on several CVS axes is never claimed.
#' @return A data.frame: `fn`, `file`, `needle` (NA for the whole file).
#' @noRd
.manifest_claims <- function(evidence) {
  kits <- .matrix_rules("kit-evidence.csv")
  wrapper <- kits$test_file[match(evidence$location, kits$path)]
  evidence$file <- ifelse(
    evidence$origin == "harness", .cvs_harness_entry(sub("[.]R$", "", basename(evidence$location))),
    ifelse(evidence$origin == "kit", wrapper, evidence$location))
  evidence <- evidence[!is.na(evidence$file) & nzchar(evidence$file) &
                         !startsWith(evidence$file, file.path("tests", "equivalence", "")), ]
  rules <- .matrix_rules("test-attribution.csv")
  # Each rule: the functions it credits and its single CVS axis (NA when several or none).
  rule_axis <- vapply(seq_len(nrow(rules)), function(r) {
    axes <- unique(.cvs_axis(matrix_rule_references(rules$reference[[r]], rules$kind[[r]])$reference,
                             matrix_rule_references(rules$reference[[r]], rules$kind[[r]])$kind)$axis)
    if (length(axes) == 1L) axes[[1L]] else NA_character_
  }, character(1))
  rule_fns <- strsplit(rules$fn, ";", fixed = TRUE)
  pairs <- unique(evidence[c("fn", "file", "origin")])
  pairs <- pairs[order(pairs$fn, pairs$file, method = "radix"), ]
  claims <- lapply(seq_len(nrow(pairs)), function(i) {
    fn <- pairs$fn[[i]]
    file <- pairs$file[[i]]
    in_file <- unique(evidence[evidence$file == file, c("fn", "item", "origin")])
    if (identical(pairs$origin[[i]], "harness")) {
      if (all(in_file$fn == fn)) return(data.frame(fn = fn, file = file, needle = NA_character_))
      wanted <- unique(in_file$item[in_file$fn == fn])
      unwanted <- unique(in_file$item[in_file$fn != fn])
      needles <- .manifest_needles(wanted, unwanted, sub("[.]R$", "", basename(file)),
                                   prefixes = sub(": .*$", ": ", wanted))
      return(data.frame(fn = fn, file = file, needle = needles))
    }
    own <- which(rules$dir == dirname(file) & rules$file == basename(file))
    credits <- vapply(own, function(r) fn %in% rule_fns[[r]] && !is.na(rule_axis[[r]]), logical(1))
    harmless <- vapply(own, function(r) identical(rule_axis[[r]], "internal"), logical(1))
    if (length(own) && all(credits | harmless)) {
      return(data.frame(fn = fn, file = file, needle = NA_character_))
    }
    titles <- unique(in_file$item[in_file$origin == "testthat" & in_file$item != "(no saved result)"])
    hit <- matrix_test_rule(dirname(file), basename(file), titles, rules)
    wanted_rule <- vapply(hit, function(r) {
      !is.na(r) && fn %in% rule_fns[[r]] && !is.na(rule_axis[[r]])
    }, logical(1))
    if (!any(wanted_rule)) return(NULL)
    needles <- .manifest_needles(titles[wanted_rule], titles[!wanted_rule], basename(file))
    data.frame(fn = fn, file = file, needle = needles)
  })
  out <- unique(do.call(rbind, claims))
  out[order(out$fn, out$file, out$needle, method = "radix", na.last = FALSE), ]
}

#' A YAML double-quoted scalar
#' @noRd
.yaml_quote <- function(x) {
  x <- gsub("\\", "\\\\", enc2utf8(as.character(x)), fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

#' YAML lines for a block list: of scalars (a character vector) or of flow
#' maps (a data.frame, one map per row, in column order)
#' @noRd
.yaml_list <- function(key, values) {
  if (is.data.frame(values)) {
    if (!nrow(values)) return(sprintf("%s: []", key))
    items <- vapply(seq_len(nrow(values)), function(i) {
      fields <- vapply(names(values), function(name) {
        sprintf("%s: %s", name, .yaml_quote(values[[name]][[i]]))
      }, character(1))
      sprintf("  - { %s }", paste(fields, collapse = ", "))
    }, character(1))
    return(c(sprintf("%s:", key), items))
  }
  if (!length(values)) return(sprintf("%s: []", key))
  c(sprintf("%s:", key), sprintf("  - %s", .yaml_quote(values)))
}

#' Every manifest, as a list per function
#' @return A named list (by function) of manifests: scalars, character vectors
#'   and data.frames, in the order they are written.
#' @noRd
manifest_build <- function() {
  summary <- .manifest_matrix_table("summary.csv")
  evidence <- .manifest_matrix_table("evidence.csv")
  cells <- .manifest_matrix_table("matrix.csv")
  stopifnot(
    "summary.csv has one row per matrix function" = setequal(summary$fn, names(.matrix_functions)),
    "every function has a family" = all(summary$fn %in% names(.manifest_family)))
  exports <- sub("^export\\((.*)\\)$", "\\1",
                 grep("^export\\(", readLines("NAMESPACE", warn = FALSE), value = TRUE))
  claims <- .manifest_claims(evidence)
  tolerances <- .manifest_tolerances()
  wired <- cells[!cells$candidate, ]
  sources <- .cvs_axis(wired$reference, wired$kind)
  wired$software <- sources$software
  papers <- wired[wired$kind == "published", ]
  papers <- merge(papers, unique(.manifest_bib_keys(papers$reference)), by = "reference")
  fns <- sort(summary$fn, method = "radix")
  manifests <- lapply(fns, function(fn) {
    row <- summary[summary$fn == fn, ]
    api <- .manifest_public_api(row$label, exports)
    own_claims <- claims[claims$fn == fn, ]
    filtered <- own_claims[!is.na(own_claims$needle), ]
    expected <- unique(wired[wired$fn == fn & !is.na(wired$software), c("software", "reference")])
    expected <- expected[order(expected$software, expected$reference, method = "radix"), ]
    own_tolerances <- tolerances[tolerances$fn == fn, c("function", "reference", "tolerance", "status")]
    list(
      runner_id = fn, family = .manifest_family[[fn]], label = row$label,
      slug = gsub("_", "-", fn), name = .manifest_rd_title(api[[1L]]),
      summary = sprintf(paste0(
        "%s. Verification strength %s in latents' verification matrix: %d source(s)",
        "%s, %d published paper(s), %d dataset(s), %d current checks of which %d failed."),
        .manifest_rd_title(api[[1L]]), row$strength, row$n_sources,
        if (nzchar(row$sources)) sprintf(" (%s)", row$sources) else "", row$n_papers,
        row$n_datasets, row$checks_current, row$failed_current),
      module_files = .manifest_module_files(api), public_api = api,
      papers = sort(unique(papers$key[papers$fn == fn]), method = "radix"),
      tolerances = own_tolerances,
      cvs_files = sort(unique(own_claims$file), method = "radix"),
      cvs_test_filter = data.frame(file = filtered$file, test_name_includes = filtered$needle,
                                   stringsAsFactors = FALSE),
      cross_software_expected = data.frame(software = expected$software, package = expected$reference,
                                           stringsAsFactors = FALSE))
  })
  names(manifests) <- fns
  attr(manifests, "papers") <- unique(papers[c("key", "reference", "author", "year")])
  manifests
}

#' YAML text of one manifest
#' @noRd
manifest_yaml <- function(manifest) {
  lines <- unlist(lapply(names(manifest), function(key) {
    value <- manifest[[key]]
    if (is.data.frame(value) || key %in% c("module_files", "public_api", "papers", "cvs_files")) {
      .yaml_list(key, value)
    } else {
      sprintf("%s: %s", key, .yaml_quote(value))
    }
  }))
  c("# Generated by equivalence/verification/build-manifests.R from the verification",
    "# matrix; do not edit by hand.", lines)
}

#' BibTeX for every cited paper, one entry per key (the longest reference
#' string of a key is its title)
#' @noRd
manifest_bib <- function(papers) {
  papers <- papers[order(papers$key, -nchar(papers$reference), papers$reference, method = "radix"), ]
  papers <- papers[!duplicated(papers$key), ]
  unlist(lapply(seq_len(nrow(papers)), function(i) {
    fields <- c(
      title = sprintf("{%s}", papers$reference[[i]]),
      author = if (!is.na(papers$author[[i]])) sprintf("{%s}", sub(" et al[.]?$", " and others", gsub(" & ", " and ", papers$author[[i]]))),
      year = if (!is.na(papers$year[[i]])) sprintf("{%s}", papers$year[[i]]),
      note = "{As recorded in latents' verification matrix (equivalence/verification-matrix/)}")
    c(sprintf("@misc{%s,", papers$key[[i]]),
      sprintf("  %s = %s,", names(fields), fields), "}", "")
  }))
}

#' Write every manifest and references.bib
#' @return The files written, invisibly.
#' @noRd
build_manifests <- function(output = .manifest_dir) {
  manifests <- manifest_build()
  modules <- file.path(output, "modules")
  dir.create(modules, recursive = TRUE, showWarnings = FALSE)
  stale <- setdiff(list.files(modules, "[.]yaml$"), paste0(names(manifests), ".yaml"))
  unlink(file.path(modules, stale))
  paths <- vapply(names(manifests), function(fn) {
    path <- file.path(modules, paste0(fn, ".yaml"))
    writeLines(enc2utf8(manifest_yaml(manifests[[fn]])), path, useBytes = TRUE)
    path
  }, character(1))
  bib <- file.path(output, "references.bib")
  writeLines(enc2utf8(manifest_bib(attr(manifests, "papers"))), bib, useBytes = TRUE)
  invisible(c(unname(paths), bib))
}

if (sys.nframe() == 0L) print(build_manifests())
