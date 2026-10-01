# Compares multilpa(family = ) group-class fits with Latent GOLD 6.1 on the
# kit cases, using the tolerances declared in README.md. Run from the project
# root: Rscript equivalence/latentgold-families/compare.R
# Writes comparison.csv and REPORT.md; exits with status 1 on any disagreement.

devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-families")
cases <- readRDS(file.path(base, "kit", "cases.rds"))

listing_value <- function(lines, label) {
  hit <- grep(label, lines, fixed = TRUE, value = TRUE)
  stopifnot("listing label must appear once" = length(hit) >= 1L)
  as.numeric(strsplit(hit[1L], "\t", fixed = TRUE)[[1L]][2L])
}

permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  smaller <- permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))
  }))
}

compare_case <- function(name) {
  case <- cases[[name]]
  data <- utils::read.delim(file.path(base, "kit", paste0(name, ".dat")))
  vars <- paste0("y", seq_len(case$d))
  fit <- multilpa(data, vars, "g", n_group_classes = case$H, family = case$family,
                  between_variance = case$between, tol = 1e-14, max_iter = 20000,
                  n_starts = 20, seed = 1)
  lines <- iconv(readLines(file.path(base, "returned", paste0(name, ".lst")),
                           warn = FALSE), "latin1", "UTF-8")
  lg_ll <- listing_value(lines, "Log-likelihood (LL)")
  lg_npar <- listing_value(lines, "Number of parameters (Npar)")
  posteriors <- utils::read.delim(
    file.path(base, "returned", paste0(name, "_posteriors.txt")),
    check.names = FALSE, fileEncoding = "latin1")
  class_columns <- grep("^GClass#[0-9]+$", names(posteriors), value = TRUE)
  lg <- unique(cbind(g = posteriors[[1L]], posteriors[class_columns]))
  stopifnot("one Latent GOLD posterior row per group" =
              nrow(lg) == fit$n_groups && !anyDuplicated(lg$g))
  mine <- fit$group_posteriors[as.character(lg$g), , drop = FALSE]
  theirs <- as.matrix(lg[class_columns])
  options <- permutations(case$H)
  gaps <- vapply(seq_len(nrow(options)), function(k) {
    max(abs(mine - theirs[, options[k, ], drop = FALSE]))
  }, numeric(1))
  data.frame(
    case = name, family = case$family, between_variance = case$between,
    quantity = c("parameter count", "log likelihood", "group posteriors"),
    latents = c(fit$n_parameters, fit$log_likelihood, NA),
    latent_gold = c(lg_npar, lg_ll, NA),
    difference = c(abs(fit$n_parameters - lg_npar), abs(fit$log_likelihood - lg_ll),
                   min(gaps)),
    tolerance = c(0, 1e-4, 1e-3),
    stringsAsFactors = FALSE)
}

comparison <- do.call(rbind, lapply(names(cases), compare_case))
comparison$agree <- comparison$difference <= comparison$tolerance
utils::write.csv(comparison, file.path(base, "comparison.csv"), row.names = FALSE)

report <- c(
  "# Latent GOLD comparison: group-class families — results", "",
  sprintf("Run %s; %s; latents %s. Tolerances declared in README.md before running.",
          format(Sys.Date()), R.version.string, utils::packageVersion("latents")), "",
  sprintf("**%d of %d quantities agree.**", sum(comparison$agree), nrow(comparison)), "",
  "| Case | Family | Between | Quantity | latents | Latent GOLD | Difference | Tolerance | Agree |",
  "|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %s | %s | %s | %.2e | %.0e | %s |", comparison$case,
          comparison$family, comparison$between_variance, comparison$quantity,
          ifelse(is.na(comparison$latents), "", format(comparison$latents, digits = 10)),
          ifelse(is.na(comparison$latent_gold), "", format(comparison$latent_gold, digits = 10)),
          comparison$difference, comparison$tolerance,
          ifelse(comparison$agree, "yes", "**no**")))
writeLines(report, file.path(base, "REPORT.md"))
print(comparison, row.names = FALSE, digits = 8)
if (!all(comparison$agree)) quit(status = 1L)
