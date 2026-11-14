# The CVS verification manifests (verification/build-manifests.R): current with
# the matrix, valid YAML with the fields CVS requires, and claiming only report
# entries the CVS reporter writes (../validation/docs/TARGET-VERIFICATION.md).

manifest_env <- function() {
  environment <- new.env(parent = globalenv())
  # The generator, like build.R, works from the project root.
  withr::with_dir("..", sys.source(file.path("equivalence", "verification", "build-manifests.R"),
                                   envir = environment))
  environment
}

# The parts of a manifest line the generator writes: a comment, `key: "scalar"`,
# `key:` / `key: []`, `  - "item"` and `  - { key: "value", ... }`.
manifest_line_ok <- function(line) {
  quoted <- "\"([^\"\\\\]|\\\\.)*\""
  patterns <- c("^# ", sprintf("^[a-z_]+: %s$", quoted), "^[a-z_]+:( \\[\\])?$", sprintf("^  - %s$", quoted),
                sprintf("^  - \\{ [a-z_]+: %s(, [a-z_]+: %s)* \\}$", quoted, quoted))
  any(vapply(patterns, grepl, logical(1), x = line, perl = TRUE))
}

test_that("the committed manifests are the ones the matrix generates", {
  generator <- manifest_env()
  output <- withr::local_tempdir()
  written <- withr::with_dir("..", generator$build_manifests(output = output))
  expect_identical(sort(basename(written)),
                   sort(c(paste0(names(generator$.matrix_functions), ".yaml"), "references.bib")))
  committed <- file.path("verification", c(file.path("modules", basename(written[-length(written)])),
                                            "references.bib"))
  expect_true(all(file.exists(committed)))
  # Deterministic: an unchanged matrix rebuilds byte for byte.
  expect_identical(lapply(written, readLines), lapply(committed, readLines))
})

test_that("every manifest parses and carries what CVS requires", {
  paths <- sort(list.files(file.path("verification", "modules"), "[.]yaml$", full.names = TRUE))
  expect_length(paths, 11L)
  lines <- unlist(lapply(paths, readLines, encoding = "UTF-8"))
  expect_true(all(vapply(lines, manifest_line_ok, logical(1))))
  skip_if_not_installed("yaml")
  bib_keys <- sub("^@misc\\{(.*),$", "\\1",
                  grep("^@misc\\{", readLines(file.path("verification", "references.bib")), value = TRUE))
  expect_false(anyDuplicated(bib_keys) > 0L)
  for_each <- lapply(paths, function(path) {
    manifest <- yaml::read_yaml(path)
    required <- c("runner_id", "family", "label", "slug", "summary", "module_files", "public_api",
                  "papers", "tolerances", "cvs_files", "cross_software_expected")
    expect_true(all(required %in% names(manifest)), info = path)
    expect_identical(paste0(manifest$runner_id, ".yaml"), basename(path))
    expect_true(manifest$family %in% c("cluster", "regression", "dynamics"))
    expect_true(all(file.exists(file.path("..", unlist(manifest$module_files)))), info = path)
    expect_true(all(unlist(manifest$papers) %in% bib_keys), info = path)
    statuses <- vapply(manifest$tolerances, `[[`, character(1), "status")
    expect_true(any(statuses == "pass"), info = path)
    software <- vapply(manifest$cross_software_expected, `[[`, character(1), "software")
    expect_true(all(software %in% c("latentgold", "Mplus", "stata", "sas", "paper", "independent")), info = path)
    filtered <- vapply(manifest$cvs_test_filter, `[[`, character(1), "file")
    expect_true(all(filtered %in% unlist(manifest$cvs_files)), info = path)
    unlist(manifest$cvs_files)
  })
  expect_gt(length(unlist(for_each)), 0L)
})

test_that("every claimed file is an entry name the CVS reporter writes", {
  source(file.path("cvs", "cvs-report.R"), local = TRUE)
  paths <- list.files(file.path("verification", "modules"), "[.]yaml$", full.names = TRUE)
  claimed <- unique(unlist(lapply(paths, function(path) {
    lines <- readLines(path)
    start <- which(lines == "cvs_files:")
    if (!length(start)) return(character())
    block <- lines[seq.int(start + 1L, length(lines))]
    block <- block[cumsum(!startsWith(block, "  - ")) == 0L]
    sub("^  - \"(.*)\"$", "\\1", block)
  })))
  # The names the runners report under: equivalence/run.R names each test file
  # file.path(<its folder>, file) (equivalence/, or tests/equivalence/ with
  # --dir=); the harness names each suite with .cvs_harness_entry(), the sweep
  # each module with .cvs_sweep_entry().
  suites <- sub("[.]R$", "", list.files(file.path("harness", "suites"), "[.]R$"))
  harness <- .cvs_from_harness(
    data.frame(suite = suites, source = "poLCA", quantity = "q", difference = 0, tolerance = 1,
               agrees = TRUE, output = NA_character_, n_values = 1L, known_open = NA_character_,
               stringsAsFactors = FALSE),
    data.frame(suite = suites, status = "ran", reason = NA_character_, stringsAsFactors = FALSE))
  modules <- sub("[.]R$", "", list.files(file.path("argument-sweep", "modules"), "[.]R$"))
  reported <- c(file.path("equivalence", list.files(".", "^test-.*[.]R$")),
                file.path("tests", "equivalence",
                          list.files(file.path("..", "tests", "equivalence"), "^test-.*[.]R$")),
                vapply(harness, `[[`, character(1), "name"), .cvs_sweep_entry(modules))
  expect_gt(length(claimed), 20L)
  expect_identical(setdiff(claimed, reported), character())
})

test_that("needles claim exactly as CVS matches them", {
  generator <- manifest_env()
  claims <- generator$.manifest_needle_claims
  expect_true(claims("Faithful VVV 2: loglik", "invariance", "faithful vvv 2: "))
  expect_true(claims("loglik", "invariance", "INVARIANCE"))          # ancestors are searched too
  expect_false(claims("loglik [2]", "x", "loglik$"))                 # `$` anchors to the title's end
  expect_true(claims("a: loglik", "x", "loglik$"))
  needles <- generator$.manifest_needles
  # A prefix that would also claim an unwanted title falls back to the whole title.
  expect_identical(needles(c("fit A: x", "fit A: y"), "fit A: z", "suite", prefixes = c("fit A: ", "fit A: ")),
                   c("fit A: x", "fit A: y"))
  # A title contained in an unwanted one is anchored.
  expect_identical(needles("agrees", "agrees with Mplus", "file"), "agrees$")
  expect_error(needles("agrees", c("agrees", "it agrees"), "file"),
               class = "latents_manifest_ambiguous_needle")
})

test_that("references get stable BibTeX keys and YAML strings are escaped", {
  generator <- manifest_env()
  keys <- generator$.manifest_bib_keys(c(
    "Masyn (2013), Oxford Handbook of Quantitative Methods, Table 25.12", "Masyn (2013)",
    "Muthén & Asparouhov (2022)", "Proust-Lima, Philipps and Liquet (2017), JSS 78(2), section 6",
    "SAS/STAT 14.2 PROC FMM Example 40.3 (Wang et al. 1996 model)"))
  expect_identical(keys$key, c("masyn2013", "masyn2013", "muthen2022", "proustlima2017", "sasstat142"))
  expect_identical(keys$year, c("2013", "2013", "2022", "2017", NA))
  bib <- generator$manifest_bib(keys)
  expect_identical(sum(startsWith(bib, "@misc{")), 4L)              # one entry per key
  expect_true("  title = {Masyn (2013), Oxford Handbook of Quantitative Methods, Table 25.12}," %in% bib)
  quoted <- generator$.yaml_quote("a \"b\" \\ c")
  expect_identical(quoted, "\"a \\\"b\\\" \\\\ c\"")
  skip_if_not_installed("yaml")
  expect_identical(yaml::yaml.load(paste("k:", quoted))$k, "a \"b\" \\ c")
})
