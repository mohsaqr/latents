pkgload::load_all(".", quiet = TRUE)
set.seed(110506)
options(width = 110)
dd <- data.frame(id = rep(1:220, each = 2), time = rep(1:2, 220))
first <- sample(0:1, 220, TRUE)
second <- ifelse(runif(220) < .8, first, 1 - first)
state <- as.vector(rbind(first, second))
dd$a <- rnorm(nrow(dd), 4 * state, .7)
dd$b <- rnorm(nrow(dd), -3 * state, .8)
dd$c <- rnorm(nrow(dd), 2 * state, .9)
warnings <- character()
fit <- withCallingHandlers(lta(dd, c("a", "b", "c"), "id", 2, "time",
  mover_stayer = TRUE, n_starts = 3, seed = 39, tol = 1e-9, cores = 1),
  warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1], conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
print(warnings)
pars <- .lta_parameters(fit)
softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
initial_mover <- softmax(c(pars$initial[1, 1, 1], 0))
initial_stayer <- softmax(c(pars$initial[1, 1, 2], 0))
transition <- rbind(softmax(c(0, pars$transition[1, 1, 1, 1])),
                    softmax(c(pars$transition[1, 1, 2, 1], 0)))
joint <- (1 - pars$group_probabilities[2]) * diag(initial_mover) %*% transition +
  pars$group_probabilities[2] * diag(initial_stayer)
print(joint)
# The two-wave observed mixture depends on the four cells of joint only.
# Decompose it into a stayer diagonal and mover matrix in several ways.
new_shares <- c(.01, .10, .20)
alternatives <- do.call(rbind, lapply(new_shares, function(share) {
  stationary <- c(.5, .5)
  mover <- (joint - share * diag(stationary)) / (1 - share)
  stopifnot(all(mover > 0))
  p <- pars
  initial <- rowSums(mover)
  moves <- mover / initial
  p$group_probabilities <- c(1 - share, share)
  p$initial[1, 1, 1] <- log(initial[1] / initial[2])
  p$initial[1, 1, 2] <- 0
  p$transition[1, 1, 1, 1] <- log(moves[1, 2] / moves[1, 1])
  p$transition[1, 1, 2, 1] <- log(moves[2, 1] / moves[2, 2])
  data.frame(stayer_share = share, LL = .lta_expectation(fit$x, fit$codes, fit$layout,
    fit$designs, p, fit$occasion_of_row)$log_likelihood)
}))
cat("Reported share:", pars$group_probabilities[2], "LL:", fit$log_likelihood,
    "converged:", fit$converged, "replicated:", fit$n_best_replicated, "\n")
print(alternatives, digits = 12)
stopifnot(max(abs(alternatives$LL - fit$log_likelihood)) < 1e-8)
saveRDS(list(fit = fit, alternatives = alternatives, warnings = warnings),
  "validation/review-2026-11-05/engines/transition-ridge.rds")
