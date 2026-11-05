suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# 1. Brute-force path sum of the three-step chain on a tiny case, plain base R.
K <- 3L
D <- rbind(c(.8,.15,.05), c(.1,.7,.2), c(.05,.25,.7))
# persons: person 1 full T=4, person 2 T=3 with a missing emission at t=2, person 3 T=2, person 4 T=4
modal <- list(c(1,1,2,3), c(2,NA,3), c(3,1), c(2,2,2,1))
z0 <- c(0.5,-1,0.2,1.3)              # initial covariate, one per person
zt <- list(c(.3,-.2,.7,.1), c(0,.5,-1), c(1.1,.4), c(-.6,.2,.9,1.5)) # transition covariate by occasion
# Build layout the package's way
ids <- rep(seq_along(modal), lengths(modal))
times <- unlist(lapply(lengths(modal), seq_len))
layout <- .multilpa_sequence_layout(ids, times, length(modal), "observed")
log_density <- log(t(D))[ifelse(is.na(unlist(modal)),1,unlist(modal)), , drop=FALSE]
log_density[is.na(unlist(modal)), ] <- 0
emission <- .multilpa_sequence_emission(log_density, layout)
designs <- .lta_designs(layout, matrix(unlist(zt), ncol=1, dimnames=list(NULL,"zt")),
                        matrix(z0[ids], ncol=1, dimnames=list(NULL,"z0")), "homogeneous")
inputs <- list(layout=layout, designs=designs, emission=emission, errors=D, modal=unlist(modal),
               n_profiles=K, informative_groups=rep(TRUE,4), n_observations=length(ids))
p0 <- 2L; p <- 2L
n_par <- p0*(K-1)+p*(K-1)*K
theta <- round(rnorm(n_par,0,.7),3)
ev <- .lta_three_step_evaluate(theta, inputs)

# independent coding: initial logits against last, transition logits against stay
init <- matrix(theta[1:(p0*(K-1))], p0, K-1)
tr <- array(theta[-(1:(p0*(K-1)))], c(p,K-1,K))
init_prob <- function(z) { e <- c(init[1,]+init[2,]*z, 0); exp(e)/sum(exp(e)) }
trans_prob <- function(k, z) {
  e <- numeric(K); others <- setdiff(1:K,k)
  e[others] <- tr[1,,k] + tr[2,,k]*z
  e[k] <- 0
  exp(e)/sum(exp(e))
}
brute <- function(i, theta_ok = TRUE) {
  m <- modal[[i]]; TT <- length(m)
  paths <- as.matrix(expand.grid(rep(list(1:K), TT)))
  sum(apply(paths, 1, function(s) {
    pr <- init_prob(z0[i])[s[1]]
    if (!is.na(m[1])) pr <- pr * D[s[1], m[1]]
    if (TT>1) for (t in 2:TT) {    # brute force by design, a tiny reference
      pr <- pr * trans_prob(s[t-1], zt[[i]][t])[s[t]]
      if (!is.na(m[t])) pr <- pr * D[s[t], m[t]]
    }
    pr
  }))
}
bf <- log(vapply(seq_along(modal), brute, 0))
cat("per-person LL package:", format(ev$log_likelihood, digits=12), "\n")
cat("per-person LL brute  :", format(bf, digits=12), "\n")
cat("max abs diff:", max(abs(ev$log_likelihood-bf)), "\n")
# staying reference: zero theta row for origin k gives equal mix; check transition prob with package coder
lt <- .lta_log_transition(cbind(1,0.7), tr)
print(round(exp(lt[1,,]),6)); print(round(rbind(trans_prob(1,.7),trans_prob(2,.7),trans_prob(3,.7)),6))
# analytic scores vs numDeriv (of the package's own LL) and vs brute force
if (requireNamespace("numDeriv", quietly=TRUE)) {
  brute_total <- function(th) {
    init <<- matrix(th[1:(p0*(K-1))], p0, K-1); tr <<- array(th[-(1:(p0*(K-1)))], c(p,K-1,K))
    sum(log(vapply(seq_along(modal), brute, 0)))
  }
  g_b <- numDeriv::grad(brute_total, theta)
  brute_total(theta)
  cat("max |score - FD of brute force|:", max(abs(colSums(ev$scores)-g_b)), "\n")
  per <- numDeriv::jacobian(function(th) .lta_three_step_evaluate(th, inputs)$log_likelihood, theta)
  cat("max |per-person score - FD|:", max(abs(ev$scores-per)), "\n")
}
