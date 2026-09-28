# Interval coverage for latent transition models (`parameter_inference()` on
# an `lta()` fit). Run from the package root:
#   Rscript validation/transition-inference-coverage.R
#
# Data: 150 sequences of five occasions, two profiles, one group class; the
# initial distribution is (0.6, 0.4) and the transition matrix has rows
# (0.85, 0.15) and (0.25, 0.75); two Gaussian indicators with profile means
# (0, 0) and (2, 1.6), unit variances. Profiles are aligned to the generating
# labels by the first indicator's mean. Coverage is of the nominal 95% Wald
# interval under the observed information and under the sandwich clustered on
# sequences.
pkgload::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260927)
replications <- 300L
truth <- c(initial_1 = 0.6, stay_1 = 0.85, stay_2 = 0.75, mean_a_1 = 0,
           mean_a_2 = 2, variance_a_1 = 1)
initial <- c(0.6, 0.4)
transition <- matrix(c(0.85, 0.15, 0.25, 0.75), 2L, byrow = TRUE)

replicate_once <- function(replication) {
  frame <- do.call(rbind, lapply(seq_len(150L), function(person) {
    path <- Reduce(function(previous, wave) {
      sample.int(2L, 1L, prob = transition[previous, ])
    }, seq_len(4L), init = sample.int(2L, 1L, prob = initial), accumulate = TRUE)
    data.frame(person = person, wave = seq_len(5L), state = path)
  }))
  frame$a <- stats::rnorm(nrow(frame), c(0, 2)[frame$state])
  frame$b <- stats::rnorm(nrow(frame), c(0, 1.6)[frame$state])
  fit <- withCallingHandlers(
    lta(frame, c("a", "b"), "person", 2L, time = "wave", n_starts = 3,
        seed = 1, tol = 1e-10, max_iter = 5000),
    latents_unconverged = function(condition) invokeRestart("muffleWarning"))
  if (!fit$converged) return(NULL)
  low <- if (fit$means[1L, "a"] < fit$means[2L, "a"]) 1L else 2L
  high <- 3L - low
  pick <- function(table) {
    row <- function(parameter, outcome, term) {
      table[table$parameter == parameter & table$outcome == outcome &
              table$term == term, ]
    }
    profile <- function(k) sprintf("profile_%d", k)
    rbind(row("initial_probability", profile(low), "group_class_1"),
          row("transition_probability", profile(low),
              sprintf("group_class_1:%s", profile(low))),
          row("transition_probability", profile(high),
              sprintf("group_class_1:%s", profile(high))),
          row("mean", profile(low), "a"), row("mean", profile(high), "a"),
          row("variance", profile(low), "a"))
  }
  observed <- pick(parameter_inference(fit))
  robust <- pick(parameter_inference(fit, vcov_type = "robust"))
  data.frame(replication = replication, parameter = names(truth),
             estimate = observed$estimate,
             se_observed = observed$standard_error,
             se_robust = robust$standard_error,
             covered_observed = observed$conf_low <= truth & truth <= observed$conf_high,
             covered_robust = robust$conf_low <= truth & truth <= robust$conf_high)
}

# One row per parameter.
summarise_coverage <- function(results) {
  table <- do.call(rbind, lapply(split(results, results$parameter), function(part) {
    data.frame(parameter = part$parameter[1L], usable = nrow(part),
               bias = mean(part$estimate) - truth[[part$parameter[1L]]],
               empirical_sd = stats::sd(part$estimate),
               mean_se_observed = mean(part$se_observed),
               mean_se_robust = mean(part$se_robust),
               coverage_observed = mean(part$covered_observed),
               coverage_robust = mean(part$covered_robust))
  }))
  table[match(names(truth), table$parameter), ]
}

results <- do.call(rbind, parallel::mclapply(
  seq_len(replications), replicate_once, mc.set.seed = TRUE,
  mc.cores = max(1L, parallel::detectCores() - 2L)))
summary_table <- summarise_coverage(results)
print(summary_table, digits = 3, row.names = FALSE)
utils::write.csv(summary_table, "validation/transition-inference-coverage.csv",
                 row.names = FALSE)
