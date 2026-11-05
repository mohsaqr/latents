suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
cat("pkgs:", vapply(c("numDeriv","mice","nnet","sandwich","poLCA","flexmix","mclust"), requireNamespace, TRUE, quietly=TRUE), "\n")
# ---- independent base-R one-step covariate model: 2 group classes, 3 profiles, 2 Gaussians
sim <- function(G=60, n=8, seed=1){ set.seed(seed)
  do.call(rbind, lapply(seq_len(G), function(g){
    w <- rnorm(1); h <- rbinom(1,1,plogis(.5*w))+1L
    z <- rnorm(n)
    eta <- cbind(c(1,-.5)[h] + .8*z, c(-.5,1)[h] - .6*z, 0)
    p <- exp(eta)/rowSums(exp(eta))
    k <- apply(p,1,function(q) sample.int(3,1,prob=q))
    data.frame(g=g,w=w,z=z,k=k,
      y1=rnorm(n,c(-2,0,2)[k],.8), y2=rnorm(n,c(1,-1,1.5)[k],.6))
  }))}
d <- sim()
fit <- multilpa(d, c("y1","y2"), "g", n_profiles=3, n_group_classes=2,
  profile_covariates="z", group_covariates="w", n_starts=4, seed=1, tol=1e-10, max_iter=5000)
print(class(fit)); print(fit$converged); print(fit$log_likelihood); print(fit$n_parameters)
print(fit$profile_coefficients); print(fit$group_coefficients)
theta <- .multilpa_cov_encode(fit)
print(length(theta))
# independent likelihood written from scratch
indep_ll <- function(means, vars, B, Gm){
  # B: (2 intercepts + 1 slope) x 2 profiles; Gm: 2 x 1 (intercept,w)
  log_dens <- sapply(1:3, function(k) dnorm(d$y1, means[k,1], sqrt(vars[k,1]), log=TRUE)+dnorm(d$y2, means[k,2], sqrt(vars[k,2]), log=TRUE))
  ids <- unique(d$g); ll <- 0
  for (g in ids){ r <- which(d$g==g); w <- d$w[r[1]]
    pg2 <- 1/(1+exp(sum(c(1,w)*Gm))); pg1 <- 1-pg2
    lg <- sapply(1:2, function(h){
      X <- cbind(as.numeric(h==1), as.numeric(h==2), d$z[r])
      eta <- cbind(X%*%B, 0); pk <- exp(eta)/rowSums(exp(eta))
      sum(log(rowSums(pk*exp(log_dens[r,,drop=FALSE]))))})
    m <- max(lg); ll <- ll + m + log(pg1*exp(lg[1]-m)+pg2*exp(lg[2]-m)) }
  ll }
cat("hand LL:", sprintf("%.8f", indep_ll(fit$means, fit$variances, fit$profile_coefficients, fit$group_coefficients)),
    " package:", sprintf("%.8f", fit$log_likelihood), "\n")
# natural-scale reproduction via a packed vector, then numDeriv SEs of the slope etc.
unpack <- function(th){ means <- matrix(th[1:6],3,2,byrow=TRUE); vars <- exp(matrix(th[7:12],3,2,byrow=TRUE))
  B <- matrix(th[13:18],3,2); Gm <- matrix(th[19:20],2,1); list(means=means,vars=vars,B=B,Gm=Gm)}
th0 <- c(as.vector(t(fit$means)), as.vector(t(log(fit$variances))), as.vector(fit$profile_coefficients), as.vector(fit$group_coefficients))
f <- function(th){ u <- unpack(th); indep_ll(u$means,u$vars,u$B,u$Gm)}
H <- numDeriv::hessian(f, th0)
V <- solve(-H)
pi <- parameter_inference(fit)
natural <- c(sqrt(diag(V))[1:6], sqrt(diag(V))[7:12]*exp(th0[7:12]), sqrt(diag(V))[13:20])
cat("max rel diff SE (independent Hessian vs package):", max(abs(natural/pi$standard_error-1)),"\n")
# hand robust: group scores numerically
gscore <- function(th){ u <- unpack(th); ids <- unique(d$g)
  log_dens <- sapply(1:3, function(k) dnorm(d$y1,u$means[k,1],sqrt(u$vars[k,1]),log=TRUE)+dnorm(d$y2,u$means[k,2],sqrt(u$vars[k,2]),log=TRUE))
  vapply(ids,function(g){r<-which(d$g==g); w<-d$w[r[1]]; pg2<-1/(1+exp(sum(c(1,w)*u$Gm))); pg1<-1-pg2
    lg<-sapply(1:2,function(h){X<-cbind(as.numeric(h==1),as.numeric(h==2),d$z[r]);eta<-cbind(X%*%u$B,0);pk<-exp(eta)/rowSums(exp(eta));sum(log(rowSums(pk*exp(log_dens[r,,drop=FALSE]))))})
    m<-max(lg); m+log(pg1*exp(lg[1]-m)+pg2*exp(lg[2]-m))},0)}
J <- numDeriv::jacobian(gscore, th0)
G <- nrow(J); R <- V%*%crossprod(J)%*%V*G/(G-1)
pr <- parameter_inference(fit, vcov_type="robust")
natr <- c(sqrt(diag(R))[1:6], sqrt(diag(R))[7:12]*exp(th0[7:12]), sqrt(diag(R))[13:20])
cat("robust max rel diff:", max(abs(natr/pr$standard_error-1)),"\n")
po <- parameter_inference(fit, vcov_type="opg"); O <- solve(crossprod(J))
nato <- c(sqrt(diag(O))[1:6], sqrt(diag(O))[7:12]*exp(th0[7:12]), sqrt(diag(O))[13:20])
cat("opg max rel diff:", max(abs(nato/po$standard_error-1)),"\n")
# ---- relabel to a different reference: coef contrasts & SE vs linear contrast of V
al <- .multilpa_cov_permute(fit, c(3L,1L,2L), c(2L,1L))
pa <- parameter_inference(al); cat("permuted LL", al$log_likelihood, "\n")
# slope SE after rebase: new reference = old profile 2. new profile_1 = old 3 vs old 2: -b2 ; new profile_2 = old 1 vs old 2: b1-b2
Vb <- V[13:18,13:18]
se1 <- sqrt(Vb[6,6]); se2 <- sqrt(Vb[3,3]+Vb[6,6]-2*Vb[3,6])
print(pa[pa$level=="profile",c("outcome","term","estimate","standard_error")])
cat("by hand z slope new profile_1:", -fit$profile_coefficients[3,2], se1, " new profile_2:", fit$profile_coefficients[3,1]-fit$profile_coefficients[3,2], se2,"\n")
print(pa[pa$level=="group",c("outcome","term","estimate","standard_error")])
cat("hand group slope flip:", -fit$group_coefficients[2,1],"\n")
