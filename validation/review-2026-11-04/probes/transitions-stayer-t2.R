# T = 2 with a stayer class: show the stayer share is arbitrary. Build other
# parameter points from the fit's implied latent joint distribution and evaluate
# them with the package's own likelihood.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 3)[s]), y2 = rnorm(T, c(1, -2)[s], .8), y3 = rnorm(T, c(0, 2)[s]))
  }))
}
vars <- c("y1","y2","y3")
d2 <- sim(300, 2, 1)
fit <- suppressWarnings(lta(d2, vars, "id", 2, "time", mover_stayer = TRUE, n_starts = 6, seed = 1))
pars <- .lta_parameters(fit)
softmax <- function(e) { e <- e - max(e); exp(e)/sum(exp(e)) }
mover_init <- softmax(c(pars$initial[1,1,1], 0)); stayer_init <- softmax(c(pars$initial[1,1,2], 0))
tr <- rbind(softmax(c(0, pars$transition[1,1,1,1])), softmax(c(pars$transition[1,1,2,1], 0)))
s <- pars$group_probabilities[2]
joint <- (1 - s) * diag(mover_init) %*% tr + s * diag(stayer_init)
cat("latent joint implied by the fit:\n"); print(round(joint, 5))
cat("reported stayer share", round(s, 4), " logLik", fit$log_likelihood, "\n")
ll_at <- function(s_new, r_new) {
  mover <- joint - s_new * diag(r_new)
  stopifnot(all(mover > 0))
  mover <- mover / (1 - s_new)
  init <- rowSums(mover); trn <- mover / init
  p <- pars
  p$group_probabilities <- c(1 - s_new, s_new)
  p$initial[1,1,1] <- log(init[1]/init[2]); p$initial[1,1,2] <- log(r_new[1]/r_new[2])
  p$transition[1,1,1,1] <- log(trn[1,2]/trn[1,1]); p$transition[1,1,2,1] <- log(trn[2,1]/trn[2,2])
  .lta_expectation(fit$x, fit$codes, fit$layout, fit$designs, p, fit$occasion_of_row)$log_likelihood
}
for (s_new in c(.05, .2, .4)) for (r1 in c(.3, .6))
  cat(sprintf("stayer share %.2f, stayer initial %.1f: logLik %.8f (diff from fit %.2e)\n", s_new, r1, ll_at(s_new, c(r1, 1-r1)), ll_at(s_new, c(r1, 1-r1)) - fit$log_likelihood))
pi <- parameter_inference(fit)
print(pi[, c("block","group_class","from","to","term","estimate","standard_error","estimated")])
cat("eigen range of the information:\n"); inf <- .lta_inference(fit)
ev <- eigen(inf$information[inf$identified, inf$identified], only.values = TRUE)$values
print(range(ev))
