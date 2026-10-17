# Latent GOLD 6.1 kit for the three-step latent Markov model of
# r3step(initial_covariates =, transition_covariates =) (Vermunt 2010; Di
# Mari, Oberski and Vermunt 2016). Each case is simulated, fitted with lta(),
# and each occasion's step-two posteriors (from its own indicators, see
# .lta_occasion_posteriors()) are written with the covariates. Latent GOLD
# reads those posteriors, assigns modally, estimates the classification
# errors itself and fits the step-three chain (`step3 modal ml`).
# Run from the package root:
#   Rscript equivalence/latentgold-step3-markov/make-kit.R
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-step3-markov")
kit_dir <- file.path(base, "kit")
wine_dir <- file.path(Sys.getenv("HOME"), ".wine-latentgold", "drive_c", "lgstep3markov")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wine_dir, showWarnings = FALSE, recursive = TRUE)

# Per case: persons, occasions, profile means of two indicators, initial
# logits (intercept, slope on w; the last profile is the reference) and, per
# origin, transition logits (intercept, slope on z) of each other profile.
cases <- list(
  m01 = list(seed = 301L, J = 400L, T = 5L, means = rbind(c(-1, -0.8), c(1, 0.8)),
             initial = rbind(c(0.3, 0.8)),
             moves = list(rbind(c(-1.2, 0.8)), rbind(c(-1.0, -0.6))),
             initial_covariates = "w", transition_covariates = "z", drop = 0),
  m02 = list(seed = 302L, J = 500L, T = 4L,
             means = rbind(c(-1.5, -1), c(0, 0.2), c(1.5, 1)),
             initial = rbind(c(0.4, 0), c(0.2, 0)),
             moves = list(rbind(c(-1.5, 0.5), c(-2.5, 0)), rbind(c(-1.8, 0), c(-1.6, 0.6)),
                          rbind(c(-2.4, 0), c(-1.4, -0.5))),
             initial_covariates = character(), transition_covariates = "z", drop = 0),
  # Unbalanced: a quarter of the later occasions are missing, so sequences
  # end at different lengths.
  m03 = list(seed = 303L, J = 450L, T = 6L, means = rbind(c(-1.1, -0.9), c(0.9, 1.1)),
             initial = rbind(c(-0.2, 0.6)),
             moves = list(rbind(c(-1.4, 0.7)), rbind(c(-1.1, -0.4))),
             initial_covariates = "w", transition_covariates = "z", drop = 0.25))

simulate <- function(case) {
  set.seed(case$seed)
  k <- nrow(case$means)
  w <- stats::rnorm(case$J)
  z <- matrix(stats::rnorm(case$J * case$T), case$J)
  softmax_draw <- function(logits) {
    p <- exp(logits - apply(logits, 1L, max))
    p <- p / rowSums(p)
    1L + rowSums(stats::runif(nrow(p)) > t(apply(p, 1L, cumsum)))
  }
  first <- cbind(vapply(seq_len(k - 1L), function(l) {
    case$initial[l, 1L] + case$initial[l, 2L] * w
  }, numeric(case$J)), 0)
  state <- matrix(NA_integer_, case$J, case$T)
  state[, 1L] <- softmax_draw(first)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(case$T)[-1L]) {  # simulation only
    logits <- t(vapply(seq_len(case$J), function(i) {
      from <- state[i, t - 1L]
      out <- numeric(k)
      out[setdiff(seq_len(k), from)] <- case$moves[[from]][, 1L] +
        case$moves[[from]][, 2L] * z[i, t]
      out
    }, numeric(k)))
    state[, t] <- softmax_draw(matrix(logits, case$J))
  }
  frame <- data.frame(id = rep(seq_len(case$J), each = case$T),
                      time = rep(seq_len(case$T), case$J),
                      w = rep(w, each = case$T), z = as.vector(t(z)))
  s <- as.vector(t(state))
  frame$y1 <- stats::rnorm(length(s), case$means[s, 1L], 1)
  frame$y2 <- stats::rnorm(length(s), case$means[s, 2L], 1)
  keep <- frame$time <= 2L | stats::runif(nrow(frame)) >= case$drop
  # A person's sequence ends at the first dropped occasion.
  ended <- stats::ave(!keep, frame$id, FUN = cumsum) > 0
  frame[!ended, ]
}

syntax <- function(name, case) {
  k <- nrow(case$means)
  posterior <- paste(sprintf("post%d", seq_len(k)), collapse = " ")
  independent <- c(case$initial_covariates, case$transition_covariates)
  initial_terms <- paste(c("1", case$initial_covariates), collapse = " + ")
  transition_terms <- paste(sprintf("(~tra) %s | State[-1]",
                                    c("1", case$transition_covariates)), collapse = " + ")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261017 sets=20 tolerance=1e-005 iterations=100;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   step3 modal ml;",
    "   output", "      parameters=last standarderrors estimatedvalues;",
    "variables",
    "   caseid id;",
    sprintf("   independent %s;", paste(unique(independent), collapse = ", ")),
    sprintf("   latent State nominal dynamic coding=last posterior=(%s);", posterior),
    "equations",
    sprintf("   State[=0] <- %s;", initial_terms),
    sprintf("   State <- %s;", transition_terms),
    "end model")
}

invisible(lapply(names(cases), function(name) {
  case <- cases[[name]]
  data <- simulate(case)
  fit <- lta(data, c("y1", "y2"), "id", n_profiles = nrow(case$means), time = "time",
             n_starts = 5, seed = 1)
  posteriors <- .lta_occasion_posteriors(fit)
  colnames(posteriors) <- sprintf("post%d", seq_len(ncol(posteriors)))
  out <- cbind(data, posteriors)
  utils::write.table(out, file.path(kit_dir, paste0(name, ".dat")), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  writeLines(syntax(name, case), file.path(kit_dir, paste0(name, ".lgs")))
  saveRDS(fit, file.path(kit_dir, paste0(name, "-fit.rds")))
  saveRDS(data, file.path(kit_dir, paste0(name, "-data.rds")))
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
file.copy(list.files(kit_dir, pattern = "[.](dat|lgs)$", full.names = TRUE), wine_dir,
          overwrite = TRUE)
cat("kit written:", length(cases), "cases\n")
