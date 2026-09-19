# Published-information-criterion arithmetic; this is not an estimation replication.
# Run from the repository root:
#   Rscript validation/houle-2026/check-fit-indices.R [path/to/supplement.docx]
# No external R packages are required. Optional DOCX input verifies the fixture
# against the original first table without storing the copyrighted document.

stopifnot(file.exists("DESCRIPTION"), dir.exists("R"))
output_dir <- file.path("validation", "houle-2026")
fixture <- utils::read.csv(file.path(output_dir, "table-s1.csv"),
                           stringsAsFactors = FALSE)
stopifnot(nrow(fixture) == 24L)

read_docx_fit_table <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, file.exists(path))
  scratch <- tempfile("houle-docx-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  utils::unzip(path, files = "word/document.xml", exdir = scratch)
  xml <- paste(readLines(file.path(scratch, "word", "document.xml"),
                        warn = FALSE, encoding = "UTF-8"), collapse = "")
  extract <- function(pattern, text) {
    regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
  }
  table <- extract("<w:tbl\\b[^>]*>.*?</w:tbl>", xml)[1L]
  rows <- extract("<w:tr\\b[^>]*>.*?</w:tr>", table)
  numeric_rows <- lapply(rows, function(row) {
    cells <- trimws(gsub("<[^>]+>", "", extract("<w:tc\\b[^>]*>.*?</w:tc>", row)))
    if (length(cells) != 11L || !grepl("^-?[0-9]+\\.[0-9]+$", cells[2L])) {
      return(NULL)
    }
    as.numeric(cells[c(2L, 3L, 5L, 6L, 7L, 8L)])
  })
  do.call(rbind, Filter(Negate(is.null), numeric_rows))
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
source_checked <- length(args) == 1L
if (source_checked) {
  original <- read_docx_fit_table(args[[1L]])
  source_columns <- c("log_likelihood", "n_parameters", "aic", "caic", "bic", "sabic")
  stopifnot(identical(dim(original), c(24L, 6L)),
            max(abs(original - as.matrix(fixture[, source_columns]))) == 0)
}

# Derive N from BIC - AIC. Each printed criterion can differ from its true
# value by 0.0005, so their difference can differ by 0.001. N=10,000 is the
# unique integer allowed by the intersection over all 24 rows.
half_unit <- 0.0005
q <- fixture$n_parameters
sample_lower <- max(exp(2 + (fixture$bic - fixture$aic - 2 * half_unit) / q))
sample_upper <- min(exp(2 + (fixture$bic - fixture$aic + 2 * half_unit) / q))
stopifnot(ceiling(sample_lower) == floor(sample_upper))
n_individuals <- ceiling(sample_lower)
stopifnot(n_individuals == 10000)

# Independent definitions, calculated without the package implementation.
penalties <- cbind(aic = 2 * q, caic = q * (log(n_individuals) + 1),
                   bic = q * log(n_individuals),
                   sabic = q * log((n_individuals + 2) / 24))
reported <- as.matrix(fixture[, colnames(penalties)])
deviance <- -2 * fixture$log_likelihood
expected <- penalties + deviance

# Check joint consistency, not just four loose tolerances. An unrounded
# deviance must lie in the LL-derived interval AND each criterion-derived
# interval simultaneously. Multiplying LL by -2 doubles its rounding radius.
lower <- pmax(deviance - 2 * half_unit,
              apply(reported - penalties - half_unit, 1L, max))
upper <- pmin(deviance + 2 * half_unit,
              apply(reported - penalties + half_unit, 1L, min))
stopifnot(all(lower <= upper),
          all(abs(fixture$caic - fixture$bic - q) < 1e-9))

# Load the current source tree, ensuring a stale installed package cannot
# determine these results. Functions are isolated in a separate environment.
package_env <- new.env(parent = globalenv())
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE),
                 sys.source, envir = package_env))

package_values <- t(vapply(seq_len(nrow(fixture)), function(i) {
  # Arithmetic-only adapter: LL and q come from the published table, and
  # posteriors are placeholders needed by information_criteria(). No model is
  # fitted, no entropy is validated, and no group count is inferred here.
  # A different dummy group count deliberately distinguishes the conventions.
  adapter <- structure(list(
    log_likelihood = fixture$log_likelihood[i], n_parameters = q[i],
    n_groups = 100L, n_observations = n_individuals,
    n_informative = n_individuals,
    group_posteriors = matrix(1, nrow = 100L, ncol = 1L),
    subject_posteriors = matrix(1, nrow = n_individuals, ncol = 1L)
  ), class = "multilpa")
  result <- package_env$information_criteria(adapter)
  selected <- result[result$convention %in% c("none", "individuals"), ]
  stats::setNames(selected$value[match(colnames(penalties), selected$criterion)],
                  colnames(penalties))
}, stats::setNames(numeric(ncol(penalties)), colnames(penalties))))
package_error <- abs(package_values - expected)
stopifnot(all(is.finite(package_values)), all(package_error < 1e-9))

results <- data.frame(
  fixture[, c("block", "classes_enumerated")],
  rounding_interval_lower = lower, rounding_interval_upper = upper,
  rounding_interval_width = upper - lower, rounding_consistent = lower <= upper,
  max_printed_difference = apply(abs(reported - expected), 1L, max),
  max_package_difference = apply(package_error, 1L, max)
)
utils::write.csv(results, file.path(output_dir, "fit-index-results.csv"), row.names = FALSE)
summary <- c(
  "Houle et al. (2026), supplemental Table S1: arithmetic validation",
  sprintf("Original DOCX fixture verification: %s", if (source_checked) "24/24 rows pass" else "not requested"),
  sprintf("Joint LL/AIC/BIC/CAIC/ABIC rounding consistency: %d/24 rows pass", sum(lower <= upper)),
  sprintf("N interval inferred from rounded BIC - AIC: [%.9f, %.9f]", sample_lower, sample_upper),
  sprintf("Unique integer N: %.0f", n_individuals),
  sprintf("Package information_criteria arithmetic: %d/96 values pass", sum(package_error < 1e-9)),
  sprintf("Maximum package/formula absolute difference: %.12g", max(package_error)),
  sprintf("Maximum printed/formula absolute difference: %.12g (rounding-compatible)", max(abs(reported - expected))),
  "Package source: current repository R/; individual sample-size convention.",
  "Scope: arithmetic only; no published-data estimation, entropy, SE, or test replication."
)
writeLines(summary, file.path(output_dir, "fit-index-results.txt"))
cat(paste(summary, collapse = "\n"), "\n")
