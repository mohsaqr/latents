suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Synthetic: 2 profiles, 4 occasions, covariates on initial and transition
n <- 250; TT <- 4
z0 <- rnorm(n); zt <- matrix(rnorm(n*TT, 0, 1), n, TT)
state <- matrix(0L, n, TT)
state[,1] <- 1L + (runif(n) < plogis(.3 + .9*z0)) * 1L   # P(1)... profile labels handled below
invisible(lapply(2:TT, function(t) {
  move <- runif(n) < plogis(-1.2 + .8*zt[,t])
  state[,t] <<- ifelse(move, 3L-state[,t-1], state[,t-1])
}))
long <- data.frame(id = rep(seq_len(n), TT), time = rep(seq_len(TT), each=n), z0 = rep(z0,TT), zt = as.vector(zt))
s <- as.vector(state)
long$a <- rnorm(n*TT, c(-1.3,1.3)[s]); long$b <- rnorm(n*TT, c(-1.1,1.1)[s]); long$c <- rnorm(n*TT, c(1,-1)[s])
# rows in a shuffled order on purpose
long_s <- long[sample(nrow(long)), ]
long <- long[order(long$id, long$time), ]
fit <- lta(long, c("a","b","c"), id="id", time="time", n_profiles=2, n_starts=3, seed=1)
cat("class:", class(fit), "\n")
res <- r3step(fit, long, initial_covariates="z0", transition_covariates="zt")
print(as.data.frame(res)); cat("LL attr", attr(res,"log_likelihood"), "\n"); print(attr(res,"classification_errors"))

# Hand: occasion posteriors from fit parameters, pooled shares, modal, D
mu <- fit$means; var <- fit$variances
X <- as.matrix(long[,c("a","b","c")])
cat("means dims"); print(dim(mu)); print(mu); print(var)
ld <- vapply(1:2, function(k) rowSums(mapply(function(col) dnorm(X[,col], mu[k,col]+ (if(!is.null(fit$centers)) 0 else 0), sqrt(var[k,col]), log=TRUE), 1:3)), numeric(nrow(X)))
# if means are centred on data use fit's indicator data
ld2 <- vapply(1:2, function(k) rowSums(mapply(function(col) dnorm(fit$indicator_data[,col], mu[k,col], sqrt(var[k,col]), log=TRUE), 1:3)), numeric(nrow(X)))
share <- colMeans(fit$subject_posteriors)
post1 <- exp(sweep(ld,2,log(share),"+")); post1 <- post1/rowSums(post1)
post2 <- exp(sweep(ld2,2,log(share),"+")); post2 <- post2/rowSums(post2)
pkg <- .lta_occasion_posteriors(fit)
cat("max diff occasion posteriors, raw X:", max(abs(pkg-post1)), " indicator_data:", max(abs(pkg-post2)), "\n")
modal <- max.col(pkg, "first")
D <- rbind(colSums(pkg[,1]*outer(modal,1:2,"==")), colSums(pkg[,2]*outer(modal,1:2,"=="))) / colSums(pkg)
cat("D diff:", max(abs(D-attr(res,"classification_errors"))), "\n")

# Independent ML in base R (forward recursion, no package helpers)
ll_i <- function(par, i) {
  # par: init (int, z0), from1->2 (int, zt), from2->1 (int, zt)
  m <- modal[long$id==i][order(long$time[long$id==i])]; z <- zt[i,]
  p1 <- plogis(par[1]+par[2]*z0[i]); pi0 <- c(p1,1-p1)
  a <- pi0*D[,m[1]]
  for (t in 2:TT) {    # recursion in time; tiny reference
    a12 <- plogis(par[3]+par[4]*z[t]); a21 <- plogis(par[5]+par[6]*z[t])
    A <- rbind(c(1-a12,a12), c(a21,1-a21))
    a <- as.vector(a %*% A)*D[,m[t]]
  }
  log(sum(a))
}
nll <- function(par) -sum(vapply(seq_len(n), function(i) ll_i(par,i), 0))
op <- optim(c(0,0,-1,0,-1,0), nll, method="BFGS", control=list(reltol=1e-13, maxit=2000))
cat("independent par:", format(op$par,digits=8), "\n")
cat("independent LL :", -op$value, " vs package ", attr(res,"log_likelihood"), "\n")
est <- res$estimate
cat("max est diff:", max(abs(est-op$par)), "\n")
if (requireNamespace("numDeriv", quietly=TRUE)) {
  H <- numDeriv::hessian(nll, op$par)
  obs <- r3step(fit, long, initial_covariates="z0", transition_covariates="zt", vcov_type="observed")
  cat("observed SE: pkg", format(obs$standard_error,digits=7), "\n   indep", format(sqrt(diag(solve(H))),digits=7),"\n")
  S <- numDeriv::jacobian(function(par) vapply(seq_len(n), function(i) ll_i(par,i), 0), op$par)
  V <- solve(H) %*% (crossprod(S)*n/(n-1)) %*% solve(H)
  cat("robust SE: pkg", format(res$standard_error,digits=7), "\n   indep", format(sqrt(diag(V)),digits=7),"\n")
}
# Row-shuffled refit: same table up to labels
fit_s <- lta(long_s, c("a","b","c"), id="id", time="time", n_profiles=2, n_starts=3, seed=1)
res_s <- r3step(fit_s, long_s, initial_covariates="z0", transition_covariates="zt")
cat("shuffled rows LL:", attr(res_s,"log_likelihood"), "\n")
print(cbind(as.data.frame(res)[,c("block","from","to","term","estimate","standard_error")], s=res_s$estimate, s_se=res_s$standard_error))
# p_adjusted family
print(all.equal(res$p_adjusted[res$term!="(Intercept)"], p.adjust(res$p_value[res$term!="(Intercept)"],"BH")))
# Equivariance in units/origin
long2 <- transform(long, z0=z0*1000+5e4, zt=zt/100-7)
res2 <- r3step(fit, long2, initial_covariates="z0", transition_covariates="zt")
slope <- res$term!="(Intercept)"
cat("slope ratio z0, zt, se ratio:", res2$estimate[slope]/res$estimate[slope], " | ", res2$standard_error[slope]/res$standard_error[slope], "\n")
cat("LL diff after rescale", attr(res2,"log_likelihood")-attr(res,"log_likelihood"), "\n")
# A general-engine fit with the same first-order model (occasion-varying transitions) : classifications
g <- lta(long, c("a","b","c"), id="id", time="time", n_profiles=2, transitions="occasion", n_starts=3, seed=1)
cat("general class:", class(g), "\n")
gp <- .lta_occasion_posteriors(g)
parms <- .lta_parameters(g)
print(str(parms$measurement,1,give.attr=FALSE))
m <- parms$measurement[[1]]
ld3 <- vapply(1:2, function(k) rowSums(mapply(function(col) dnorm(X[,col], m$means[k,col]+g$centers[col], sqrt(m$variances[k,col]), log=TRUE), 1:3)), numeric(nrow(X)))
sh <- colMeans(g$expectation$subject_posteriors)
p3 <- exp(sweep(ld3,2,log(sh),"+")); p3 <- p3/rowSums(p3)
cat("general occasion posteriors vs hand:", max(abs(gp-p3)), "\n")
print(as.data.frame(r3step(g, long, transition_covariates="zt"))[,c("block","from","to","term","estimate","standard_error")])
