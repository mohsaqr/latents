# Latent GOLD 6.1 kit for three_step(occasion =, method = "ml"): a distal
# outcome of an lta() chain, related to the profile at one occasion. Each
# case is simulated, fitted with lta(), and each occasion's step-two
# posteriors (from its own indicators, .lta_occasion_posteriors(), as in the
# step3-markov kit) are written with the outcome, which sits on the record of
# the chosen occasion and is missing on the others. Latent GOLD reads the
# posteriors, assigns modally, estimates the classification errors itself and
# fits the step-three chain with the outcome depending on the state
# (`step3 modal ml` on a dynamic latent variable).
#
# Two syntaxes per case:
# - <case>.lgs holds the chain's initial and transition logits fixed at the
#   lta() estimates (`~fix`), the convention of three_step(): this is the
#   comparison of record, to Latent GOLD's printed precision;
# - <case>-free.lgs lets Latent GOLD re-estimate the chain with the outcome.
#   Its outcome means are expected to differ from three_step()'s in the
#   second or third decimal; it is kept to measure that gap, not to pass.
#
# STATUS: prepared 2026-11-15 on a machine without Latent GOLD; not yet run.
# The `~fix` labels and the order of the values in their braces follow the
# Latent GOLD 6.1 syntax manual as remembered, not as executed: on the first
# run, check that the listing prints the fixed chain values (and no standard
# errors for them) before trusting the comparison.
#
# Run from the package root:
#   Rscript equivalence/latentgold-step3-distal/make-kit.R
# then run every .lgs in kit/ with Latent GOLD 6.1 (batch: lg61 <file>.lgs
# /b /o <file>.lst) and copy the listings to returned/.
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-step3-distal")
kit_dir <- file.path(base, "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base, "returned"), showWarnings = FALSE)

# Per case: persons, occasions, profile means of two indicators, initial
# probabilities, a transition matrix (rows: from), the chosen occasion and
# the outcome's mean and standard deviation in each profile.
cases <- list(
  d01 = list(seed = 401L, J = 500L, T = 4L, means = rbind(c(-1, -0.8), c(1, 0.8)),
             initial = c(0.55, 0.45), transition = rbind(c(0.8, 0.2), c(0.15, 0.85)),
             occasion = 3L, outcome_mean = c(0, 1.5), outcome_sd = c(1, 0.8)),
  d02 = list(seed = 402L, J = 600L, T = 4L,
             means = rbind(c(-1.5, -1), c(0, 0.2), c(1.5, 1)),
             initial = c(0.4, 0.35, 0.25),
             transition = rbind(c(0.8, 0.15, 0.05), c(0.1, 0.8, 0.1), c(0.05, 0.15, 0.8)),
             occasion = 2L, outcome_mean = c(-1, 0, 1.2), outcome_sd = c(1, 1.2, 0.9)))

simulate <- function(case) {
  set.seed(case$seed)
  k <- nrow(case$means)
  draw <- function(p) 1L + rowSums(stats::runif(nrow(p)) > t(apply(p, 1L, cumsum)))
  state <- matrix(NA_integer_, case$J, case$T)
  state[, 1L] <- draw(matrix(case$initial, case$J, k, byrow = TRUE))
  for (t in seq_len(case$T)[-1L]) {  # simulation only
    state[, t] <- draw(case$transition[state[, t - 1L], , drop = FALSE])
  }
  frame <- data.frame(id = rep(seq_len(case$J), each = case$T),
                      time = rep(seq_len(case$T), case$J))
  s <- as.vector(t(state))
  frame$y1 <- stats::rnorm(length(s), case$means[s, 1L], 1)
  frame$y2 <- stats::rnorm(length(s), case$means[s, 2L], 1)
  at <- state[, case$occasion]
  frame$distal <- rep(stats::rnorm(case$J, case$outcome_mean[at], case$outcome_sd[at]),
                      each = case$T)
  frame
}

# Latent GOLD's codings: initial logits against the last state
# (coding=last), transition logits of moving against staying (~tra).
chain_logits <- function(fit) {
  initial <- as.vector(fit$initial_probabilities[1L, ])
  transition <- matrix(fit$transition_probabilities[, , 1L], length(initial))
  k <- length(initial)
  moves <- unlist(lapply(seq_len(k), function(from) {
    log(transition[from, -from] / transition[from, from])
  }))
  list(initial = log(initial[-k] / initial[k]), transition = moves)
}

syntax <- function(name, case, fit, fixed) {
  k <- nrow(case$means)
  posterior <- paste(sprintf("post%d", seq_len(k)), collapse = " ")
  logits <- chain_logits(fit)
  format_values <- function(v) paste(sprintf("%.12f", v), collapse = " ")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s%s';", name, if (fixed) "" else "-free"), "options",
    "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261115 sets=20 tolerance=1e-005 iterations=100;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing includeall;",
    "   step3 modal ml;",
    "   output", "      parameters=last standarderrors estimatedvalues;",
    "variables",
    "   caseid id;",
    "   dependent distal continuous;",
    sprintf("   latent State nominal dynamic coding=last posterior=(%s);", posterior),
    "equations",
    if (fixed) "   State[=0] <- (a~fix) 1;" else "   State[=0] <- 1;",
    if (fixed) "   State <- (b~fix) (~tra) 1 | State[-1];" else
      "   State <- (~tra) 1 | State[-1];",
    "   distal <- 1 | State;",
    "   distal | State;",
    if (fixed) sprintf("   a = {%s};", format_values(logits$initial)),
    if (fixed) sprintf("   b = {%s};", format_values(logits$transition)),
    "end model")
}

results <- do.call(rbind, lapply(names(cases), function(name) {
  case <- cases[[name]]
  data <- simulate(case)
  fit <- lta(data, c("y1", "y2"), "id", n_profiles = nrow(case$means), time = "time",
             n_starts = 5, seed = 1)
  posteriors <- .lta_occasion_posteriors(fit)
  colnames(posteriors) <- sprintf("post%d", seq_len(ncol(posteriors)))
  out <- cbind(data[c("id", "time", "y1", "y2")], posteriors)
  # The outcome on the chosen occasion's record only ("." is missing).
  out$distal <- ifelse(data$time == case$occasion, sprintf("%.15g", data$distal), ".")
  utils::write.table(out, file.path(kit_dir, paste0(name, ".dat")), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  writeLines(syntax(name, case, fit, TRUE), file.path(kit_dir, paste0(name, ".lgs")))
  writeLines(syntax(name, case, fit, FALSE), file.path(kit_dir, paste0(name, "-free.lgs")))
  saveRDS(fit, file.path(kit_dir, paste0(name, "-fit.rds")))
  saveRDS(data, file.path(kit_dir, paste0(name, "-data.rds")))
  # latents' own answer at kit time, observed-information and sandwich
  # standard errors, for the person who runs Latent GOLD to compare by eye.
  ours <- three_step(fit, data, "distal", occasion = case$occasion, method = "ml")
  data.frame(case = name, occasion = case$occasion, profile = ours$class,
             mean = ours$estimate, variance = ours$variance,
             standard_error_sandwich = ours$standard_error)
}))
utils::write.csv(results, file.path(kit_dir, "latents-results.csv"), row.names = FALSE)
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat("kit written:", length(cases), "cases\n")
