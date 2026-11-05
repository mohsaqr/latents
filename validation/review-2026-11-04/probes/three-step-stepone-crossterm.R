suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
stopifnot(requireNamespace("numDeriv", quietly=TRUE))
# One data set per separation. The cross term of the stacked M-estimator, and a bootstrap reference.
run <- function(sep, n=800, seed=1, B=300) {
  set.seed(seed)
  tr <- 1L+(runif(n)<.45)
  d <- data.frame(a=rnorm(n,c(-sep,sep)[tr]), b=rnorm(n,c(-sep,sep)[tr]), y=rnorm(n,c(0,1)[tr],1))
  fit <- lpa(d,c("a","b"),n_profiles=2,model="VVI",n_starts=5,seed=1,tol=1e-11,max_iter=20000)
  X <- as.matrix(d[,c("a","b")])
  point <- c(as.vector(t(fit$means)), as.vector(t(log(sqrt(fit$variances)))), qlogis(fit$profile_probabilities[1,1]))
  state <- function(th, X.=X) {
    mu <- matrix(th[1:4],2,byrow=TRUE); s <- matrix(exp(th[5:8]),2,byrow=TRUE); p <- plogis(th[9])
    l <- cbind(log(p)+dnorm(X.[,1],mu[1,1],s[1,1],log=TRUE)+dnorm(X.[,2],mu[1,2],s[1,2],log=TRUE),
               log(1-p)+dnorm(X.[,1],mu[2,1],s[2,1],log=TRUE)+dnorm(X.[,2],mu[2,2],s[2,2],log=TRUE))
    m <- apply(l,1,max); list(lli=m+log(rowSums(exp(l-m))), post=exp(l-m)/rowSums(exp(l-m)))
  }
  S <- numDeriv::jacobian(function(th) state(th)$lli, point)
  H <- numDeriv::hessian(function(th) sum(state(th)$lli), point)
  V1 <- solve(-H)
  W <- max.col(state(point)$post,"first")
  D_of <- function(P,W.=W) t(vapply(1:2,function(t) vapply(1:2,function(s.) sum(P[W.==s.,t])/sum(P[,t]),0),c(0,0)))
  bch <- function(th){ P <- state(th)$post; w <- solve(D_of(P))[W,]; as.vector(crossprod(w,d$y))/colSums(w) }
  J <- numDeriv::jacobian(bch, point)
  P0 <- state(point)$post; w <- solve(D_of(P0))[W,]; Tk <- colSums(w); mu <- as.vector(crossprod(w,d$y))/Tk
  psi <- w*outer(d$y,mu,"-")/matrix(Tk,n,2,byrow=TRUE)
  V3 <- crossprod(psi)
  inc <- J%*%V1%*%t(J)
  cross <- J%*%V1%*%crossprod(S,psi)
  full <- V3+inc+cross+t(cross)
  f <- three_step(fit,d,"y"); e <- three_step(fit,d,"y",step_one="estimated")
  cat(sprintf("\nsep %.2f, mean entropy-R2 %.3f  n=%d\n",sep,1+sum(P0*log(pmax(P0,1e-300)))/(n*log(2)),n))
  cat("pkg SE fixed     ",format(f$standard_error,digits=5),"\n")
  cat("pkg SE estimated ",format(e$standard_error,digits=5),"\n")
  cat("hand sqrt(V3+JV1J')",format(sqrt(diag(V3+inc))*sqrt(n/n),digits=5),"\n")
  cat("with cross term   ",format(sqrt(diag(full)),digits=5),"\n")
  cat("cross contribution to variance", format(2*diag(cross),digits=4)," vs JV1J' ",format(diag(inc),digits=4),"\n")
  # Nonparametric bootstrap with warm-started EM (one-unit resampling, refitting step one; modal reassigned)
  em <- function(Xb, th) { 
    ll_old <- -Inf
    for (it in seq_len(2000)) { # EM on a small base-R reference; each iterate needs the last
      st <- state(th,Xb); r <- st$post; ll <- sum(st$lli)
      if (abs(ll-ll_old)<1e-8*(1+abs(ll))) break
      ll_old <- ll
      nk <- colSums(r); mu <- crossprod(r,Xb)/nk
      s2 <- rbind(colSums(r[,1]*sweep(Xb,2,mu[1,])^2)/nk[1], colSums(r[,2]*sweep(Xb,2,mu[2,])^2)/nk[2])
      th <- c(as.vector(t(mu)), as.vector(t(log(sqrt(s2)))), qlogis(nk[1]/nrow(Xb)))
    }
    th
  }
  draws <- t(vapply(seq_len(B), function(b) {
    idx <- sample.int(n,n,TRUE); Xb <- X[idx,]; yb <- d$y[idx]
    th <- em(Xb, point); P <- state(th,Xb)$post; Wb <- max.col(P,"first")
    w <- solve(D_of(P,Wb))[Wb,]; as.vector(crossprod(w,yb))/colSums(w)
  }, c(0,0)))
  # also: bootstrap with step one held at the point estimate (assignments and D recomputed from the fixed theta)
  fixed_draws <- t(vapply(seq_len(B), function(b) {
    idx <- sample.int(n,n,TRUE); w. <- w[idx,]; as.vector(crossprod(w.,d$y[idx]))/colSums(w.)
  }, c(0,0)))
  cat("bootstrap SD, step one refitted   ",format(apply(draws,2,sd),digits=5),"\n")
  cat("bootstrap SD, weights held fixed  ",format(apply(fixed_draws,2,sd),digits=5),"\n")
}
run(0.8); run(0.6); run(0.8,seed=2)
