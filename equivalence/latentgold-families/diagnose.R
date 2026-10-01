# Separates quadrature error from optimizer differences: evaluates the exact
# latents likelihood at Latent GOLD's own (printed, 4-decimal) estimates.
# Run from the project root after compare.R's inputs exist.
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-families")
cases <- readRDS(file.path(base, "kit", "cases.rds"))

numbers <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  suppressWarnings(as.numeric(fields))
}

lg_parameters <- function(name, case) {
  lines <- iconv(readLines(file.path(base, "returned", paste0(name, ".lst")),
                           warn = FALSE), "latin1", "UTF-8")
  vars <- paste0("y", seq_len(case$d))
  profile <- which(startsWith(lines, "Profile"))[1L]
  size <- numbers(lines[profile + grep("^Size", lines[profile:length(lines)])[1L] - 1L])
  weights <- size[is.finite(size)][seq(1L, by = 2L, length.out = case$H)]
  means <- t(vapply(vars, function(v) {
    at <- profile + which(trimws(sapply(strsplit(lines[profile:length(lines)], "\t"),
                                       `[`, 1L)) == v)[1L] - 1L
    row <- numbers(lines[at + 1L])
    row[is.finite(row)][seq(1L, by = 2L, length.out = case$H)]
  }, numeric(case$H)))
  variance_rows <- function(start_label, prefix) {
    start <- which(lines == start_label)[1L]
    block <- lines[(start + 2L):(start + 60L)]
    block <- block[seq_len(which(!grepl("^[A-Za-z]", block))[1L] - 1L)]
    lapply(vars, function(v) {
      label <- if (prefix == "F") sub("^y", "F", v) else v
      rows <- block[startsWith(block, paste0(label, "\t")) & !grepl("(chol)", block, fixed = TRUE) &
                      !grepl("<->", block, fixed = TRUE)]
      vapply(rows, function(r) {
        values <- numbers(r)
        values[is.finite(values)][1L]
      }, numeric(1), USE.NAMES = FALSE)
    })
  }
  within <- variance_rows("Variances\t", "y")
  between <- variance_rows("Variances / Covariances continuous latent\t", "F")
  expand <- function(per_var) {
    vapply(per_var, function(v) rep_len(v, case$H), numeric(case$H))
  }
  list(means = matrix(t(means), case$H, case$d), between = matrix(expand(between), case$H, case$d),
       within = matrix(expand(within), case$H, case$d), weights = weights / sum(weights))
}

rows <- lapply(names(cases), function(name) {
  case <- cases[[name]]
  data <- utils::read.delim(file.path(base, "kit", paste0(name, ".dat")))
  vars <- paste0("y", seq_len(case$d))
  stats <- .additive_prepare(data, vars, "g")
  parameters <- lg_parameters(name, case)
  exact_at_lg <- .additive_expectation(stats, parameters)$log_likelihood
  comparison <- utils::read.csv(file.path(base, "comparison.csv"))
  ll <- subset(comparison, case == name & quantity == "log likelihood")
  data.frame(case = name, latents_max = ll$latents, lg_reported = ll$latent_gold,
             exact_at_lg_estimates = exact_at_lg,
             quadrature_error = ll$latent_gold - exact_at_lg,
             latents_minus_exact_at_lg = ll$latents - exact_at_lg)
})
result <- do.call(rbind, rows)
print(result, digits = 10, row.names = FALSE)
utils::write.csv(result, file.path(base, "diagnose.csv"), row.names = FALSE)
