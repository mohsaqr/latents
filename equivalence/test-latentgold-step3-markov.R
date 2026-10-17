# The three-step latent Markov model (r3step(initial_covariates =,
# transition_covariates =)) against Latent GOLD 6.1's `step3 modal ml` on a
# dynamic latent variable. Latent GOLD read each occasion's step-two
# posteriors from the kit, assigned modally, estimated the classification
# errors itself and fitted the chain with dummy coding against the last
# profile (initial) and against staying (`~tra`, transitions), the coding
# lta() and r3step() use. Its listings are kept in returned/, so this runs
# offline; rebuild with latentgold-step3-markov/make-kit.R and rerun Latent
# GOLD if the kit changes. Estimates and standard errors are printed to four
# decimals, the log likelihood to four.

base <- file.path("latentgold-step3-markov")
cases <- readRDS(file.path(base, "kit", "cases.rds"))

lg_listing <- function(name) {
  lines <- iconv(readLines(file.path(base, "returned", paste0(name, ".lst")), warn = FALSE),
                 "latin1", "UTF-8")
  sub("\r$", "", lines)
}

lg_log_likelihood <- function(lines) {
  row <- grep("^Log-likelihood \\(LL\\)", lines, value = TRUE)[1L]
  as.numeric(strsplit(row, "\t")[[1L]][2L])
}

# Rows of the Regression Parameters table with an estimated coefficient.
lg_parameters <- function(lines) {
  start <- which(trimws(lines) == "Regression Parameters")[1L] + 2L
  end <- which(trimws(lines) == "Paired Comparisons")[1L] - 1L
  rows <- strsplit(lines[start:end], "\t")
  rows <- rows[vapply(rows, function(r) length(r) >= 7L && r[2L] == "<-", logical(1))]
  frame <- do.call(rbind, lapply(rows, function(r) {
    profile <- as.integer(sub(".*\\(([0-9]+)\\)$", "\\1", r[1L]))
    initial <- startsWith(r[1L], "State[=0]")
    data.frame(block = if (initial) "initial" else "transition",
               from = if (initial) NA_character_ else
                 paste0("profile_", sub(".*\\(([0-9]+)\\)$", "\\1", r[5L])),
               to = paste0("profile_", profile),
               term = if (identical(trimws(r[3L]), "1")) "(Intercept)" else trimws(r[3L]),
               estimate = as.numeric(r[6L]), standard_error = suppressWarnings(as.numeric(r[7L])),
               stringsAsFactors = FALSE)
  }))
  # Reference cells are printed as 0 with no standard error.
  frame[is.finite(frame$standard_error), , drop = FALSE]
}

invisible(lapply(names(cases), function(name) {
  test_that(sprintf("three-step latent Markov model matches Latent GOLD: %s", name), {
    case <- cases[[name]]
    fit <- readRDS(file.path(base, "kit", paste0(name, "-fit.rds")))
    data <- readRDS(file.path(base, "kit", paste0(name, "-data.rds")))
    ours <- r3step(fit, data, initial_covariates = case$initial_covariates,
                   transition_covariates = case$transition_covariates,
                   vcov_type = "observed")
    lines <- lg_listing(name)
    # Absolute gaps against Latent GOLD's printed rounding.
    expect_lt(abs(attr(ours, "log_likelihood") - lg_log_likelihood(lines)), 5.1e-5)
    theirs <- lg_parameters(lines)
    expect_identical(nrow(theirs), nrow(ours))
    key <- function(frame) paste(frame$block, frame$from, frame$to, frame$term)
    matched <- ours[match(key(theirs), key(ours)), ]
    expect_false(anyNA(matched$estimate))
    expect_lt(max(abs(matched$estimate - theirs$estimate)), 5.1e-5)
    expect_lt(max(abs(matched$standard_error - theirs$standard_error)), 5.1e-5)
  })
}))
