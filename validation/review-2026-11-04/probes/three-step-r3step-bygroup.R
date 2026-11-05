suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
G <- 80; m <- 12
gc <- rep(1:2, length.out=G)
g <- rep(seq_len(G), each=m)
x1 <- rnorm(G*m); x2 <- rnorm(G*m, 50, 8)
eta <- c(1.2,-1.2)[gc[g]] + .8*x1 - .03*(x2-50)
truth <- 1L + (runif(G*m) < plogis(eta))
d <- data.frame(g=g, x1=x1, x2=x2,
  a=rnorm(G*m, c(-1.4,1.4)[truth]), b=rnorm(G*m, c(-1.2,1.2)[truth]))
# Unsorted rows, non-numeric id order
d$g <- paste0("s", 100-d$g)
d <- d[sample(nrow(d)),]
fit <- multilpa(d, c("a","b"), "g", n_profiles=2, n_group_classes=2, n_starts=5, seed=1)
print(fit$profile_probabilities); print(fit$group_probabilities)
r <- r3step(fit, d, c("x1","x2"), by_group_class=TRUE)
print(as.data.frame(r)); cat("LL", attr(r,"log_likelihood"),"\n")
ob <- r3step(fit, d, c("x1","x2"), by_group_class=TRUE, vcov_type="observed")

# Independent likelihood in base R
pieces <- .multilpa_level_assignments(fit,"individuals")
post <- fit$subject_posteriors; W <- max.col(post,"first")
cat("modal equal to package:", all(W==fit$subject_profiles),"\n")
D <- t(vapply(1:2, function(t) vapply(1:2, function(s) sum(post[W==s,t])/sum(post[,t]), 0), c(0,0)))
gi <- match(d$g, unique(d$g))
cat("group_index same:", identical(as.integer(gi), as.integer(fit$group_index)),"\n")
nll_g <- function(par) {
  # alpha_1, alpha_2 (one free class), b1, b2, gamma
  a <- par[1:2]; b <- par[3:4]; gam <- par[5]
  etah <- c(plogis(gam), 1-plogis(gam))
  lg <- vapply(1:2, function(h) {
    p1 <- plogis(a[h]+b[1]*d$x1+b[2]*d$x2)
    lik <- p1*D[1,W]+(1-p1)*D[2,W]
    as.vector(rowsum(log(lik), gi, reorder=FALSE)) + log(etah[h])
  }, numeric(G))
  mx <- apply(lg,1,max)
  mx + log(rowSums(exp(lg-mx)))
}
nll <- function(par) -sum(nll_g(par))
start <- r$estimate
op <- optim(start+rnorm(5,0,.05), nll, method="BFGS", control=list(reltol=1e-14,maxit=2000,parscale=c(1,1,1,.05,1)))
cat("independent est:", format(op$par,digits=8), "\nLL",-op$value,"\n pkg est:", format(r$estimate,digits=8),"\n")
if (requireNamespace("numDeriv", quietly=TRUE)) {
  H <- numDeriv::hessian(nll, r$estimate)
  S <- numDeriv::jacobian(nll_g, r$estimate)
  Vo <- solve(H); Vr <- Vo %*% (crossprod(S)*G/(G-1)) %*% Vo
  cat("obs SE pkg", format(ob$standard_error,digits=7),"\n   indep  ", format(sqrt(diag(Vo)),digits=7),"\n")
  cat("rob SE pkg", format(r$standard_error,digits=7),"\n   indep  ", format(sqrt(diag(Vr)),digits=7),"\n")
}
# Equivariance
d2 <- transform(d, x2=x2*100+1e5)
r2 <- r3step(fit, d2, c("x1","x2"), by_group_class=TRUE)
cat("LL diff rescale", attr(r2,"log_likelihood")-attr(r,"log_likelihood"),"\n")
print(cbind(r[,c("outcome","term","estimate")], r2=r2$estimate, ratio=r2$estimate/r$estimate, seratio=r2$standard_error/r$standard_error))
# Multi-start: other starts reach the same max?  Start from random points
others <- vapply(1:6, function(i) { o <- optim(rnorm(5,0,1), nll, method="BFGS", control=list(reltol=1e-12,maxit=2000,parscale=c(1,1,1,.05,1))); -o$value },0)
cat("independent multi-start LLs:", format(others,digits=9),"\n")
# H=1 reduces to pooled
fit1 <- multilpa(d, c("a","b"), "g", n_profiles=2, n_group_classes=1, n_starts=5, seed=1)
rp <- r3step(fit1, d, c("x1","x2"))
rg <- r3step(fit1, d, c("x1","x2"), by_group_class=TRUE)
print(as.data.frame(rp)[,c("outcome","term","estimate","standard_error","p_adjusted")])
print(as.data.frame(rg)[,c("level","outcome","term","estimate","standard_error","p_adjusted")])
# three profiles and three group classes: gradient FD at random point
fit3 <- multilpa(d, c("a","b"), "g", n_profiles=3, n_group_classes=2, n_starts=5, seed=1)
piec <- .multilpa_level_assignments(fit3,"individuals")
model <- .multilpa_r3step_group_model(fit3, as.matrix(d[,c("x1","x2")]), piec, .multilpa_error_matrix(piec))
th <- model$start+rnorm(length(model$start),0,.2)
if (requireNamespace("numDeriv", quietly=TRUE)) cat("grad FD diff", max(abs(model$gradient(th)-numDeriv::grad(model$objective,th))),"\n")
print(model$start); print(fit3$profile_probabilities)
# Pooled r3step: 3 classes, hand likelihood
rr <- r3step(fit3, d, "x1")
print(as.data.frame(rr)); print(attr(rr,"reference_class"))
pi3 <- .multilpa_level_assignments(fit3,"individuals"); D3 <- .multilpa_error_matrix(pi3)
nl3 <- function(par) { e <- cbind(par[1]+par[2]*d$x1, par[3]+par[4]*d$x1, 0); pr <- exp(e)/rowSums(exp(e)); -sum(log(rowSums(pr*t(D3)[pi3$modal,]))) }
if (requireNamespace("numDeriv", quietly=TRUE)) {
  o3 <- optim(rr$estimate+.02, nl3, method="BFGS", control=list(reltol=1e-14,maxit=2000))
  cat("3-class independent", format(o3$par,digits=7), "pkg", format(rr$estimate,digits=7),"\n")
  Hr <- numDeriv::hessian(nl3, rr$estimate)
  Sr <- numDeriv::jacobian(function(par){ e <- cbind(par[1]+par[2]*d$x1, par[3]+par[4]*d$x1, 0); pr <- exp(e)/rowSums(exp(e)); log(rowSums(pr*t(D3)[pi3$modal,]))}, rr$estimate)
  Sc <- rowsum(Sr, gi, reorder=FALSE)
  V <- solve(Hr)%*%(crossprod(Sc)*G/(G-1))%*%solve(Hr)
  cat("rob SE pkg", format(rr$standard_error,digits=7),"\n   indep  ", format(sqrt(diag(V)),digits=7),"\n")
}
# group-level r3step alignment with shuffled rows
w <- rnorm(G); names(w) <- unique(d$g); d$w <- w[d$g]
gq <- r3step(fit, d, "w", level="groups")
print(as.data.frame(gq))
p <- .multilpa_level_assignments(fit,"groups"); Dg <- .multilpa_error_matrix(p)
cat("group_values order == unique(d$g):", identical(as.character(fit$group_values), unique(d$g)),"\n")
wv <- w[as.character(p$units)]
nlg <- function(par){ p1 <- plogis(par[1]+par[2]*wv); -sum(log(p1*Dg[1,p$modal]+(1-p1)*Dg[2,p$modal])) }
print(optim(c(0,0), nlg, method="BFGS", control=list(reltol=1e-14))$par)
