# Latent GOLD 6.1 kit for the three-step latent Markov model with group
# classes: r3step(initial_covariates =, transition_covariates =) on
# lta(n_group_classes = 2). Each case is simulated from a mixture of two
# Markov chains (a stable and a mobile group class), fitted with lta(), and
# each occasion's step-two posteriors (.lta_occasion_posteriors()) are
# written with the covariates. Latent GOLD reads the posteriors of the
# dynamic State, assigns modally, estimates the classification errors and
# fits the step-three model (`step3 modal ml`) with a static latent GClass
# whose classes have their own initial and transition logits: the model
# r3step() fits.
#
# STATUS: the syntax has NOT been run. Latent GOLD is not installed on the
# machine where this kit was built (2026-11-15); returned/ is empty and
# equivalence/test-latentgold-step3-markov-groups.R skips until the listings
# are there. Whether `step3` accepts a latent variable without `posterior=`
# (GClass) is to be confirmed on the first run.
# Run from the package root:
#   Rscript equivalence/latentgold-step3-markov-groups/make-kit.R
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-step3-markov-groups")
kit_dir <- file.path(base, "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)

# Per case: persons, occasions, profile means of two indicators, group-class
# shares, and per group class the initial logit (intercept, slope on w; the
# last profile is the reference) and per origin the log odds of moving
# (intercept, slope on z).
cases <- list(
  g01 = list(seed = 401L, J = 500L, T = 5L, means = rbind(c(-1.3, -1.1), c(1.3, 1.1)),
             shares = c(0.5, 0.5), initial = rbind(c(1, 0.8), c(-1, 0.8)),
             moves = list(rbind(c(-2.2, 0.8), c(-2.0, -0.6)),
                          rbind(c(-0.5, 0.8), c(-0.3, -0.6))),
             initial_covariates = "w", transition_covariates = "z"),
  g02 = list(seed = 402L, J = 500L, T = 6L, means = rbind(c(-1.2, -1), c(1.2, 1)),
             shares = c(0.6, 0.4), initial = rbind(c(0.5, 0), c(-0.5, 0)),
             moves = list(rbind(c(-2.5, 0.6), c(-2.3, 0)),
                          rbind(c(-0.7, 0.6), c(-0.6, 0))),
             initial_covariates = character(), transition_covariates = "z"))

simulate <- function(case) {
  set.seed(case$seed)
  g <- 1L + (stats::runif(case$J) > case$shares[1L])
  w <- stats::rnorm(case$J)
  z <- matrix(stats::rnorm(case$J * case$T), case$J)
  state <- matrix(NA_integer_, case$J, case$T)
  first <- case$initial[g, 1L] + case$initial[g, 2L] * w
  state[, 1L] <- ifelse(stats::runif(case$J) < stats::plogis(first), 1L, 2L)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(case$T)[-1L]) {  # simulation only
    from <- state[, t - 1L]
    coefficients <- t(vapply(seq_len(case$J), function(i) case$moves[[g[i]]][from[i], ],
                             numeric(2)))
    moved <- stats::runif(case$J) < stats::plogis(coefficients[, 1L] +
                                                    coefficients[, 2L] * z[, t])
    state[, t] <- ifelse(moved, 3L - from, from)
  }
  frame <- data.frame(id = rep(seq_len(case$J), each = case$T),
                      time = rep(seq_len(case$T), case$J),
                      w = rep(w, each = case$T), z = as.vector(t(z)))
  s <- as.vector(t(state))
  frame$y1 <- stats::rnorm(length(s), case$means[s, 1L], 1)
  frame$y2 <- stats::rnorm(length(s), case$means[s, 2L], 1)
  frame
}

syntax <- function(name, case) {
  posterior <- paste(sprintf("post%d", seq_len(nrow(case$means))), collapse = " ")
  independent <- c(case$initial_covariates, case$transition_covariates)
  initial_terms <- paste(sprintf("%s | GClass", c("1", case$initial_covariates)),
                         collapse = " + ")
  transition_terms <- paste(sprintf("(~tra) %s | State[-1] GClass",
                                    c("1", case$transition_covariates)), collapse = " + ")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261115 sets=50 tolerance=1e-005 iterations=100;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   step3 modal ml;",
    "   output", "      parameters=last standarderrors estimatedvalues;",
    "variables",
    "   caseid id;",
    sprintf("   independent %s;", paste(unique(independent), collapse = ", ")),
    "   latent GClass nominal 2 coding=last,",
    sprintf("      State nominal dynamic coding=last posterior=(%s);", posterior),
    "equations",
    "   GClass <- 1;",
    sprintf("   State[=0] <- %s;", initial_terms),
    sprintf("   State <- %s;", transition_terms),
    "end model")
}

invisible(lapply(names(cases), function(name) {
  case <- cases[[name]]
  data <- simulate(case)
  fit <- lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time",
             n_group_classes = 2, n_starts = 5, seed = 1)
  posteriors <- .lta_occasion_posteriors(fit)
  colnames(posteriors) <- sprintf("post%d", seq_len(ncol(posteriors)))
  utils::write.table(cbind(data, posteriors), file.path(kit_dir, paste0(name, ".dat")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(syntax(name, case), file.path(kit_dir, paste0(name, ".lgs")))
  saveRDS(fit, file.path(kit_dir, paste0(name, "-fit.rds")))
  saveRDS(data, file.path(kit_dir, paste0(name, "-data.rds")))
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat("kit written:", length(cases), "cases\n")
