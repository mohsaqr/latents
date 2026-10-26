# Global-maximum audit (Equivalence 2.0, E2.4).
#
# Every Monte Carlo case is fitted from five seeds with the package's default
# number of starts (10) and with 20, and each fit's log likelihood is set
# against the best of Mplus and Latent GOLD. A fit below the reference by more
# than 5e-4 stopped at a local maximum. The result is written to
# audit-maxima.csv beside this script; test-montecarlo.R checks it.
#
# Run from the project root:
#   Rscript equivalence/montecarlo/audit-maxima.R

source(file.path("validation", "fork-runtime.R"))

#' Fit every case from several seeds and starts, against the references
#' @return One row per case, seed and number of starts.
audit_maxima <- function(seeds = 1:5, starts = c(10L, 20L)) {
  cases <- montecarlo_cases()
  cores <- max(1L, min(8L, parallel::detectCores() - 1L))
  grid <- expand.grid(case = names(cases), seed = seeds, n_starts = starts,
                      stringsAsFactors = FALSE)
  grid$log_likelihood <- unlist(parallel::mclapply(seq_len(nrow(grid)), function(i) {
    fit <- suppressWarnings(montecarlo_fit(cases[[grid$case[i]]], seed = grid$seed[i],
                                           n_starts = grid$n_starts[i]))
    fit$log_likelihood
  }, mc.cores = .validation_fork_cores(cores), mc.preschedule = FALSE))
  references <- vapply(names(cases), function(name)
    max(montecarlo_mplus(name, cases[[name]])$log_likelihood,
        montecarlo_latentgold(name, cases[[name]])$log_likelihood), numeric(1))
  grid$config <- vapply(cases[grid$case], `[[`, character(1), "config")
  grid$reference <- references[grid$case]
  grid$gap <- grid$log_likelihood - grid$reference
  grid$reached <- grid$gap > -5e-4
  grid
}

#' Share of fits reaching the reference maximum, by configuration and starts
audit_summary <- function(audit) {
  reached <- aggregate(reached ~ config + n_starts, data = audit, FUN = mean)
  worst <- aggregate(gap ~ config + n_starts, data = audit, FUN = min)
  merge(reached, worst, by = c("config", "n_starts"))
}

suppressMessages(pkgload::load_all(".", quiet = TRUE))
setwd("equivalence")
invisible(lapply(c("helper-mplus-classification.R", "helper-montecarlo.R"), source))
audit <- audit_maxima()
utils::write.csv(audit, file.path("montecarlo", "audit-maxima.csv"), row.names = FALSE)
print(audit_summary(audit))
