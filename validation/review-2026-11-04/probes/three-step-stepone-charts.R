suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
stopifnot(requireNamespace("numDeriv", quietly=TRUE))
n <- 500
tr <- 1L+(runif(n)<.45)
d <- data.frame(a=rnorm(n,c(-1,1)[tr],c(1,1.2)[tr]), y=rnorm(n,c(0,1.5)[tr],1))
d$b <- .5*d$a+rnorm(n,c(-.6,.6)[tr],c(.9,.8)[tr])
d$cc <- factor(sapply(tr, function(t) sample(c("p","q","r"),1,prob=if(t==1) c(.6,.3,.1) else c(.15,.35,.5))))
d$u1 <- rbinom(n,1,c(.25,.75)[tr]); d$u2 <- rbinom(n,1,c(.7,.3)[tr]); d$u3 <- rbinom(n,1,c(.3,.65)[tr]); d$u4 <- rbinom(n,1,c(.65,.2)[tr])
d$cnt <- rpois(n,c(1,3)[tr])
d$yb <- factor(ifelse(runif(n)<c(.3,.7)[tr],"yes","no"))

check <- function(label, fit, post_fun, point) {
  ll <- function(th) post_fun(th)$ll
  V1 <- solve(-numDeriv::hessian(ll, point))
  P0 <- post_fun(point)$post
  cat("\n==",label,": base-point posterior diff", max(abs(P0-fit$subject_posteriors)), " LL ",fit$log_likelihood-ll(point), "\n")
  W <- max.col(P0,"first")
  stopifnot(all(W==fit$subject_profiles))
  errfun <- function(P) t(vapply(1:2,function(t) vapply(1:2,function(s) sum(P[W==s,t])/sum(P[,t]),0),c(0,0)))
  est_bch <- function(th){ P <- post_fun(th)$post; w <- solve(errfun(P))[W,]; as.vector(crossprod(w,d$y))/colSums(w) }
  est_prop <- function(th){ P <- post_fun(th)$post; as.vector(crossprod(P,d$y))/colSums(P) }
  est_ml <- function(th){ P <- post_fun(th)$post; wt <- t(errfun(P)[,W])*matrix(colMeans(P),n,2,byrow=TRUE)
    mu <- start_ml[1:2]; v <- exp(2*start_ml[3:4]); ch <- Inf; it <- 0L
    while (ch > 1e-13 && it < 20000L) { it <- it+1L
      j <- wt*cbind(dnorm(d$y,mu[1],sqrt(v[1])),dnorm(d$y,mu[2],sqrt(v[2]))); r <- j/rowSums(j)
      m2 <- as.vector(crossprod(r,d$y))/colSums(r); v2 <- colSums(r*outer(d$y,m2,"-")^2)/colSums(r)
      ch <- max(abs(c(m2-mu,v2-v))); mu <- m2; v <- v2 }
    mu }
  fx <- three_step(fit,d,"y",method="ml",variance="varying")
  start_ml <<- c(fx$estimate, log(sqrt(fx$variance)))
  for_method <- list(bch=est_bch, proportional=est_prop, ml=est_ml)
  invisible(lapply(names(for_method), function(m) {
    f <- three_step(fit,d,"y",method=m); fp <- three_step(fit,d,"y",method=m,contrast="pairs")
    e <- three_step(fit,d,"y",method=m,step_one="estimated"); ep <- three_step(fit,d,"y",method=m,step_one="estimated",contrast="pairs")
    J <- numDeriv::jacobian(for_method[[m]], point)
    inc <- J%*%V1%*%t(J)
    cat(sprintf("%-13s estimated SE %s | hand %s | pair SE %.6f hand %.6f\n", m,
       paste(format(e$standard_error,digits=6),collapse=" "),
       paste(format(sqrt(f$standard_error^2+diag(inc)),digits=6),collapse=" "),
       ep$standard_error,
       sqrt(fp$standard_error^2 + inc[1,1]+inc[2,2]-2*inc[1,2])))
  }))
}
# VVV: independent Cholesky coordinates
fit <- lpa(d,c("a","b"),n_profiles=2,model="VVV",n_starts=5,seed=1,tol=1e-12,max_iter=20000)
Sg <- lapply(1:2,function(k){ if(length(dim(fit$covariances))==3) fit$covariances[,,k] else fit$covariances[[k]] })
print(Sg)
ch <- lapply(Sg,chol)
point <- c(as.vector(t(fit$means)), log(ch[[1]][1,1]),ch[[1]][1,2],log(ch[[1]][2,2]), log(ch[[2]][1,1]),ch[[2]][1,2],log(ch[[2]][2,2]), qlogis(fit$profile_probabilities[1,1]))
X <- as.matrix(d[,c("a","b")])
cen <- if (!is.null(fit$centers)) fit$centers else c(0,0)
ld_mvn <- function(mu,R){ z <- backsolve(R, t(X)-mu, transpose=TRUE); -colSums(z^2)/2 - sum(log(diag(R))) - log(2*pi) }
vvv <- function(th){
  mu <- matrix(th[1:4],2,byrow=TRUE)
  R1 <- rbind(c(exp(th[5]),th[6]),c(0,exp(th[7]))); R2 <- rbind(c(exp(th[8]),th[9]),c(0,exp(th[10])))
  p <- plogis(th[11])
  l <- cbind(log(p)+ld_mvn(mu[1,],R1), log(1-p)+ld_mvn(mu[2,],R2))
  m <- apply(l,1,max); list(ll=sum(m+log(rowSums(exp(l-m)))), post=exp(l-m)/rowSums(exp(l-m)))
}
check("VVV", fit, vvv, point)
# LCA: 4 binary, class logits
fl <- lca(d,c("u1","u2","u3","u4"),n_classes=2,n_starts=5,seed=1,tol=1e-12,max_iter=20000)
print(names(fl))
str(fl$item_probabilities %||% fl$categorical_probabilities %||% fl$response_probabilities, 2)
pm <- fl$item_probabilities %||% fl$categorical_probabilities %||% fl$response_probabilities
pp <- t(vapply(1:2,function(k) vapply(names(pm),function(v){ pm[[v]][k,2] },0),numeric(4)))
pointl <- c(as.vector(qlogis(pp)), qlogis(fl$profile_probabilities[1,1]))
U <- as.matrix(d[,c("u1","u2","u3","u4")])
lcaf <- function(th){
  q <- matrix(plogis(th[1:8]),2); p <- plogis(th[9])
  l <- cbind(log(p)+rowSums(dbinom(U,1,matrix(q[1,],n,4,byrow=TRUE),log=TRUE)), log(1-p)+rowSums(dbinom(U,1,matrix(q[2,],n,4,byrow=TRUE),log=TRUE)))
  m <- apply(l,1,max); list(ll=sum(m+log(rowSums(exp(l-m)))), post=exp(l-m)/rowSums(exp(l-m)))
}
check("LCA binary", fl, lcaf, pointl)
# Mixed: continuous + 3-category + Poisson count
fm <- multilpa(transform(d,id=seq_len(n)),c("a","cc","cnt"),n_profiles=2,categorical="cc",count="cnt",n_starts=5,seed=1,tol=1e-12,max_iter=20000)
print(grep("count|rate|lambda|poisson",names(fm),value=TRUE))
str(fm$count_means %||% fm$count_rates %||% fm$rates %||% fm$extra_parameters,2)
str(fm$categorical_probabilities %||% fm$response_probabilities,2)
