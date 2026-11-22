# Writes tests/equivalence/fixtures/mplus-avepp.json from two Mplus 9 runs kept in
# equivalence/features/vlmr/returned/ (a three-class latent profile model with
# equal variances on 500 rows, and a three-class latent class model of six
# binary items on 1000 rows), and copies their data files to
# tests/equivalence/data/. For each case: the H0 log likelihood, the relative
# entropy, the estimated class counts (from the estimated model) and the AvePP
# matrix (average latent class probabilities for most likely class membership,
# row = most likely class, column = class), as Mplus printed them.
#
# Run from the package root: Rscript tests/equivalence/fixtures/make-mplus-avepp.R

source_dir <- file.path("equivalence", "features", "vlmr", "returned")
cases <- c(eei_k3poor = "eei_k3poor_01_k3", lca_k3 = "lca_k3_01_k3")
stopifnot("run from the package root" = file.exists("DESCRIPTION"), dir.exists(source_dir))

number_after <- function(lines, pattern) {
  at <- grep(pattern, lines)[1L]
  stopifnot("pattern found" = !is.na(at))
  as.numeric(sub(sprintf(".*%s\\s+", pattern), "", lines[at]))
}

parse_case <- function(stem) {
  lines <- readLines(file.path(source_dir, paste0(stem, ".out")), warn = FALSE)
  version <- trimws(sub("^Mplus VERSION", "", lines[grep("^Mplus VERSION", lines)[1L]]))
  counts_at <- grep("^BASED ON THE ESTIMATED MODEL", lines)[1L]
  count_rows <- lines[counts_at + 5:7]
  counts <- t(vapply(strsplit(trimws(count_rows), "\\s+"), as.numeric, numeric(3)))
  avepp_at <- grep("^Average Latent Class Probabilities for Most Likely", lines)[1L]
  avepp_rows <- lines[avepp_at + 5:7]
  avepp <- t(vapply(strsplit(trimws(avepp_rows), "\\s+"), as.numeric, numeric(4)))[, -1L]
  cells <- expand.grid(class = 1:3, assigned_class = 1:3)
  list(version = version,
       log_likelihood = number_after(lines, "H0 Value"),
       relative_entropy = number_after(lines, "Entropy"),
       class_counts = data.frame(class = counts[, 1L], estimated_n = counts[, 2L], proportion = counts[, 3L]),
       avepp = data.frame(assigned_class = cells$assigned_class, class = cells$class,
                          average_posterior = avepp[cbind(cells$assigned_class, cells$class)]))
}

data_dir <- file.path("tests", "equivalence", "data")
parsed <- lapply(cases, parse_case)
inputs <- lapply(names(cases), \(case) {
  target <- file.path(data_dir, sprintf("mplus-%s.dat", case))
  stopifnot(file.copy(file.path(source_dir, paste0(cases[[case]], ".dat")), target, overwrite = TRUE))
  list(path = target, sha256 = unname(tools::sha256sum(target)))
})
versions <- unique(vapply(parsed, `[[`, character(1), "version"))
stopifnot("one Mplus version" = length(versions) == 1L)
fixture <- c(lapply(parsed, \(case) case[setdiff(names(case), "version")]),
             list(provenance = list(
               program = "Mplus", version = versions,
               generated_by = "tests/equivalence/fixtures/make-mplus-avepp.R",
               generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
               data = inputs,
               source = sprintf("%s/%s.out (and .inp)", source_dir, unname(cases)))))
jsonlite::write_json(fixture, file.path("tests", "equivalence", "fixtures", "mplus-avepp.json"),
                     auto_unbox = TRUE, digits = NA, pretty = TRUE)
