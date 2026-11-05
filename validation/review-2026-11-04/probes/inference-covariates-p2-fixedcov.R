suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
G <- 60; n <- 6
d <- do.call(rbind, lapply(seq_len(G), function(g){ z <- rnorm(n); k <- rbinom(n,1,plogis(.3+1*z))+1L
  q <- rbind(c(.8,.15,.05),c(.1,.2,.7)); 
  data.frame(g=g,z=z,k=k,
   a=sapply(k,function(kk) sample(c("x","y","z"),1,prob=q[kk,])),
   b=sapply(k,function(kk) sample(c("x","y","z"),1,prob=q[kk,])),
   c=sapply(k,function(kk) sample(c("x","y"),1,prob=q[kk,1:2]/sum(q[kk,1:2]))),
   y=rnorm(n,c(-1,1)[k],1))}))
for (v in c("a","b","c")) d[[v]] <- factor(d[[v]])
base <- multilpa(d, c("a","b","c","y"), "g", n_profiles=2, n_group_classes=1, categorical=c("a","b","c"), n_starts=3, seed=1, tol=1e-10)
st <- starting_values(base, what="measurement")
str(st, max.level=1)
fit <- multilpa(d, c("a","b","c","y"), "g", n_profiles=2, n_group_classes=1, categorical=c("a","b","c"),
                profile_covariates="z", start=st, fixed="measurement", n_starts=1, seed=1, tol=1e-10, max_iter=5000)
cat("class:",class(fit),"fixed:",fit$fixed,"n_par:",fit$n_parameters,"\n")
print(coef(fit, "unconstrained"))
pi <- parameter_inference(fit); print(pi[,c("level","outcome","term","parameter","estimate","standard_error")])
vu <- vcov(fit, scale="unconstrained"); print(dim(vu)); print(rownames(vu))
# independent conditional SE: only (intercept, slope) free
ld <- with(list(), {
  rp <- fit$response_probabilities
  ll <- sapply(1:2, function(k) log(rp$a[k,as.character(d$a)]) + log(rp$b[k,as.character(d$b)]) + log(rp$c[k,as.character(d$c)]) + dnorm(d$y, fit$means[k,1], sqrt(fit$variances[k,1]), log=TRUE)); ll })
f <- function(b){ eta <- b[1]+b[2]*d$z; p <- plogis(eta); sum(log(p*exp(ld[,1])+(1-p)*exp(ld[,2]))) }
cf <- as.vector(fit$profile_coefficients); cat("hand LL",f(cf)," pkg ",fit$log_likelihood,"\n")
H <- numDeriv::hessian(f, cf); V <- solve(-H); cat("conditional SE by hand:", sqrt(diag(V)),"\n")
print(numDeriv::grad(f,cf))
# Same fit without categorical: gaussian only for contrast
cat("package SE:", pi$standard_error, "ratio pkg/hand:", pi$standard_error/sqrt(diag(V)), "\n")
print(all.equal(unname(fit$means[,1]), unname(st$means[,1])))
print(round(cbind(rp=unlist(fit$response_probabilities), st=unlist(st$response_probabilities)),6)[1:4,])
# joint Hessian SE (everything but gaussian held free) i.e. what the code effectively inverts
free_resp <- function(th){ # th: 10 response logits + 2 membership, gaussian held
  dec <- function(v,cats) { v <- matrix(v,2,cats-1,byrow=TRUE); m <- cbind(exp(v),1); m/rowSums(m) }
  pa <- dec(th[1:4],3); pb <- dec(th[5:8],3); pc <- dec(th[9:10],2)
  lv <- sapply(1:2, function(k) log(pa[k,as.integer(d$a)])+log(pb[k,as.integer(d$b)])+log(pc[k,as.integer(d$c)])+dnorm(d$y,fit$means[k,1],sqrt(fit$variances[k,1]),log=TRUE))
  p <- plogis(th[11]+th[12]*d$z); sum(log(p*exp(lv[,1])+(1-p)*exp(lv[,2]))) }
cu <- coef(fit,"unconstrained")[-(1:4)]
cat("LL at packed:", free_resp(cu), "\n")
Hj <- numDeriv::hessian(free_resp, cu); Vj <- solve(-Hj)
cat("SE if the response logits are treated as free:", sqrt(diag(Vj))[11:12], "\n")
cat("vcov(unconstrained) diag sqrt membership:", sqrt(diag(vu))[11:12], "\n")
cat("vcov natural:\n"); print(vcov(fit))
print(max(abs(sqrt(diag(vu))[1:10]/sqrt(diag(Vj))[1:10]-1)))
# robust too
pr <- parameter_inference(fit, vcov_type="robust"); print(pr[,c("term","standard_error")])
