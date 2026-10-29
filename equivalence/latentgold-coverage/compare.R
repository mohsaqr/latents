# Compares latents with the retained Latent GOLD 6.1 listings in returned/.
# Run from the project root: Rscript equivalence/latentgold-coverage/compare.R
# Writes comparison.csv (one row per compared quantity, with the function it
# checks and its status), search.csv (latents' own search against Latent GOLD's
# maximum for the enumeration steps) and RUN.txt beside this script. Exits
# non-zero when a quantity disagrees without a registered known difference to
# explain it, or when a registered difference no longer matches anything.
# Latent GOLD is not needed: its output is kept in returned/.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "latentgold-coverage")
source(file.path(directory, "compare-functions.R"))

#' Mark disagreements explained by the register; record stale entries
#' @noRd
apply_known_differences <- function(comparison, known) {
  explained <- vapply(seq_len(nrow(comparison)), function(i) {
    hits <- which(known$case == comparison$case[[i]] &
                    vapply(known$quantity, grepl, logical(1), x = comparison$quantity[[i]]))
    if (length(hits) && !comparison$agrees[[i]]) hits[[1L]] else NA_integer_
  }, integer(1))
  comparison$status <- ifelse(comparison$agrees, "agrees",
                              ifelse(is.na(explained), "disagrees", "known difference"))
  comparison$known_reason <- ifelse(is.na(explained), "", known$reason[explained])
  list(comparison = comparison, stale = known[setdiff(seq_len(nrow(known)), explained), ])
}

#' One row per case: how many quantities agree, and the worst agreeing gap
#' @noRd
case_summary <- function(comparison) {
  do.call(rbind, lapply(split(comparison, comparison$case), function(rows) {
    agreeing <- rows$difference[rows$status == "agrees" & rows$tolerance > 0]
    data.frame(case = rows$case[[1L]], fn = rows$fn[[1L]], compared = nrow(rows),
               agreed = sum(rows$status == "agrees"),
               known = sum(rows$status == "known difference"),
               worst_difference = if (length(agreeing)) max(agreeing) else NA_real_)
  }))
}

#' Run the comparison, write its files, and report
#' @return `TRUE` when nothing disagrees unexplained and no entry is stale.
run_comparison <- function(directory) {
  cases <- readRDS(file.path(directory, "kit", "cases.rds"))
  known <- utils::read.csv(file.path(directory, "known-differences.csv"), stringsAsFactors = FALSE)
  raw <- compare_all(cases, file.path(directory, "returned"))
  result <- apply_known_differences(raw, known)
  comparison <- result$comparison
  utils::write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
  utils::write.csv(attr(raw, "search"), file.path(directory, "search.csv"), row.names = FALSE)
  writeLines(c(sprintf("run: %s", format(Sys.Date())),
               sprintf("latents: %s", utils::packageVersion("latents")),
               sprintf("R: %s", R.version.string)),
             file.path(directory, "RUN.txt"))
  print(case_summary(comparison), row.names = FALSE)
  print(attr(raw, "search"), row.names = FALSE)
  cat(sprintf("\n%d quantities compared: %d agree, %d known differences, %d disagree\n",
              nrow(comparison), sum(comparison$status == "agrees"),
              sum(comparison$status == "known difference"), sum(comparison$status == "disagrees")))
  if (nrow(result$stale)) {
    cat("Registered differences that no longer match a disagreement:\n")
    print(result$stale, row.names = FALSE)
  }
  disagreeing <- subset(comparison, status == "disagrees")
  if (nrow(disagreeing)) print(disagreeing, row.names = FALSE)
  !nrow(disagreeing) && !nrow(result$stale)
}

if (!run_comparison(directory)) quit(status = 1L, save = "no")
