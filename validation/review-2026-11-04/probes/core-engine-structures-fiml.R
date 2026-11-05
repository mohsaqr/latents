suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Independent likelihood (plain base R), independent parameter count (numerical
# rank of the map from a hand parameterization) and stationarity under FIML
# (optim from the fit in the structure's own coordinates) for all 14 structures.
n <- 150; d <- 3; K <- 2
z <- rbinom(n, 1, .4)
S1 <- matrix(c(1,.5,.2,.5,1.5,.3,.2,.3,.8),3); S2 <- matrix(c(2,-.6,.1,-.6,1,.2,.1,.2,.5),3)
X <- t(vapply(z, \(g) if (g) c(3,0,1)+drop(rnorm(3)%*%chol(S1)) else c(0,2,-1)+drop(rnorm(3)%*%chol(S2)), numeric(3)))
colnames(X) <- c("a","b","c")
Xm <- X; Xm[sample(length(X), round(.12*length(X)))] <- NA
Xm <- Xm[rowSums(!is.na(Xm))>0,]
dat <- as.data.frame(Xm)

ldmvn <- function(row, mu, S) {
  o <- which(!is.na(row)); if (!length(o)) return(0)
  R <- chol(S[o,o,drop=FALSE]); e <- backsolve(R, row[o]-mu[o], transpose=TRUE)
  -.5*(length(o)*log(2*pi)+2*sum(log(diag(R)))+sum(e^2))
}
hand_ll <- function(pi_, M, Ss, data = Xm) {
  sum(apply(data, 1, function(r) {
    l <- vapply(seq_along(pi_), \(k) log(pi_[k])+ldmvn(r, M[k,], Ss[[k]]), 0)
    m <- max(l); m+log(sum(exp(l-m)))
  }))
}
hand_count <- function(s, k, d) { # by hand, volume + shape + orientation, written from the decomposition
  tab <- list(EII=1,VII=k,EEI=d,VEI=k+d-1,EVI=1+k*(d-1),VVI=k*d,
    EEE=d*(d+1)/2,VEE=k+d*(d+1)/2-1,EVE=1+k*(d-1)+d*(d-1)/2,VVE=k*d+d*(d-1)/2,
    EEV=1+(d-1)+k*d*(d-1)/2,VEV=k+(d-1)+k*d*(d-1)/2,EVV=1+k*(d*(d+1)/2-1),VVV=k*d*(d+1)/2)
  tab[[s]]+k*d+k-1
}
expm_skew <- function(u) { A <- matrix(0,3,3); A[lower.tri(A)] <- u; A <- A-t(A)
  # matrix exponential by scaling and squaring through series
  B <- A/1024; E <- diag(3); Tm <- diag(3); for (i in 1:12) { Tm <- Tm%*%B/i; E <- E+Tm }
  for (i in 1:10) E <- E%*%E; E }
# Parameterization around the fit, local coordinates theta = 0 at the fit
param <- function(s, fit_covs) {
  v <- substr(s,1,1); sh <- substr(s,2,2); o <- substr(s,3,3)
  dec <- lapply(fit_covs, function(S) { e <- eigen(S, symmetric=TRUE); list(l=prod(e$values)^(1/3), A=e$values/prod(e$values)^(1/3), D=e$vectors) })
  if (o=="I") { dec <- lapply(fit_covs, function(S) { a <- diag(S); l <- prod(a)^(1/3); list(l=l,A=a/l,D=diag(3)) }) }
  if (o=="E") { D0 <- dec[[1]]$D; dec <- lapply(fit_covs, function(S) { a <- diag(t(D0)%*%S%*%D0); l <- prod(a)^(1/3); list(l=l,A=a/l,D=D0) }) }
  nv <- if (v=="E") 1 else K; ns <- if (sh=="I") 0 else if (sh=="E") 2 else 2*K; nr <- if (o=="I") 0 else if (o=="E") 3 else 3*K
  list(n=nv+ns+nr, build=function(t) {
    tv <- t[seq_len(nv)]; ts <- t[nv+seq_len(ns)]; tr <- t[nv+ns+seq_len(nr)]
    lapply(seq_len(K), function(k) {
      l <- dec[[k]]$l*exp(if (v=="E") tv else tv[k])
      a <- dec[[k]]$A
      if (sh!="I") { da <- if (sh=="E") ts else ts[(k-1)*2+1:2]; a <- a*exp(c(da,-sum(da))) }
      Dk <- dec[[k]]$D
      if (nr>0) Dk <- Dk%*%expm_skew(if (o=="E") tr else tr[(k-1)*3+1:3])
      m <- l*Dk%*%diag(a)%*%t(Dk); (m+t(m))/2
    })
  })
}
res <- do.call(rbind, lapply(.multilpa_structures(), function(s) {
  fit <- withCallingHandlers(lpa(dat, c("a","b","c"), n_profiles=K, model=s, missing="fiml", n_starts=6, seed=3,
             tol=1e-11, max_iter=20000), warning=function(w) invokeRestart("muffleWarning"))
  Ss <- if (!is.null(fit$covariances)) lapply(1:K, \(k) matrix(fit$covariances[,,k],3,3)) else lapply(1:K, \(k) diag(fit$variances[k,]))
  pi_ <- drop(fit$profile_probabilities); M <- unname(fit$means)
  ll0 <- hand_ll(pi_, M, Ss)
  P <- param(s, Ss)
  obj <- function(t) { # local coordinates: logit pi, means, spread
    p1 <- plogis(qlogis(pi_[1])+t[1]); -hand_ll(c(p1,1-p1), M+matrix(t[2:7],2,3), P$build(t[-(1:7)]))
  }
  op <- optim(rep(0, 7+P$n), obj, method="BFGS", control=list(maxit=300, reltol=1e-12))
  # count: numerical rank of the Jacobian of (pi, means, distinct cov entries)
  implied <- function(t) { Sm <- P$build(t); c(unlist(lapply(Sm, \(S) S[lower.tri(S,TRUE)]))) }
  J <- vapply(seq_len(P$n), function(i) { h <- 1e-5; e <- numeric(P$n); e[i] <- h; (implied(e)-implied(-e))/(2*h) }, numeric(12))
  rk <- qr(J, tol=1e-6)$rank
  data.frame(s=s, ll_fit=fit$log_likelihood, ll_hand=ll0, diff=fit$log_likelihood-ll0,
             optim_gain=(-op$value)-ll0, converged=fit$converged, n_par=fit$n_parameters,
             n_hand=hand_count(s,K,d), n_rank=rk+K*d+K-1, iters=fit$iterations)
}))
print(res, digits=8)
# Sanity of the stationarity check itself: a fit stopped early must show a positive gain.
san <- do.call(rbind, lapply(c("VVI","EVE","VEV"), function(s) {
  fit <- withCallingHandlers(lpa(dat, c("a","b","c"), n_profiles=K, model=s, missing="fiml", n_starts=1, seed=3, max_iter=3L, acceleration="none"), warning=function(w) invokeRestart("muffleWarning"))
  Ss <- if (!is.null(fit$covariances)) lapply(1:K, \(k) matrix(fit$covariances[,,k],3,3)) else lapply(1:K, \(k) diag(fit$variances[k,]))
  pi_ <- drop(fit$profile_probabilities); M <- unname(fit$means); ll0 <- hand_ll(pi_, M, Ss); P <- param(s, Ss)
  obj <- function(t) { p1 <- plogis(qlogis(pi_[1])+t[1]); -hand_ll(c(p1,1-p1), M+matrix(t[2:7],2,3), P$build(t[-(1:7)])) }
  op <- optim(rep(0, 7+P$n), obj, method="BFGS", control=list(maxit=300, reltol=1e-12))
  data.frame(s=s, stopped_early_ll=ll0, optim_gain=(-op$value)-ll0)
}))
print(san, digits=6)
