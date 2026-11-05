suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
stopifnot(requireNamespace("numDeriv", quietly=TRUE))
n <- 500; tr <- sample.int(3,n,TRUE,c(.4,.35,.25)); mus <- rbind(c(-1.5,-1),c(1.5,-1),c(0,1.5))
d <- data.frame(a=rnorm(n,mus[tr,1]),b=rnorm(n,mus[tr,2]),y=rnorm(n,c(0,1,3)[tr],c(1,1.3,.8)[tr]),x=rnorm(n)+tr/2)
d$cat <- factor(vapply(tr,function(t) sample(c("p","q","r"),1,prob=rbind(c(.6,.3,.1),c(.2,.5,.3),c(.1,.2,.7))[t,]),""))
f <- lpa(d,c("a","b"),n_profiles=3,n_starts=5,seed=1,tol=1e-10)
pc <- .multilpa_level_assignments(f,"individuals")
rel <- function(a,b) max(abs(a-b)/(abs(b)+1e-8))
invisible(lapply(c(1,1e-4,1e6), function(u) {
  y <- d$y*u
  for_var <- lapply(c("varying","equal"), function(v) {
    fit <- .multilpa_step_ml_fit(pc,y,v)
    ch <- fit$outcome_coordinates; y0 <- (y-ch$center)/ch$scale
    t0 <- fit; t0$means <- ch$means; t0$variances <- fit$variances/ch$scale^2
    th <- .multilpa_step_ml_coordinates(t0,v); sc <- .multilpa_step_ml_scale(t0,v)
    score <- function(p) -colSums(.multilpa_step_ml_scores(t0,y0,v,p))
    obj <- function(p) -.multilpa_step_ml_log_likelihood(t0,y0,v,p)
    g <- numDeriv::grad(obj,th+.01); h_a <- .multilpa_step_ml_hessian(score,th,sc); h_n <- numDeriv::hessian(obj,th)
    cat(sprintf("ML unit %g %-7s |score at estimate| %.2e; grad vs FD rel %.1e; Hessian vs numDeriv rel %.1e\n",u,v,max(abs(score(th))),rel(score(th+.01),g),rel(h_a,h_n)))
  })
}))
# Categorical ML: information closed form vs numDeriv
ind <- outer(as.integer(d$cat),1:3,"==")*1
cf <- .multilpa_step_categorical_fit(pc,ind)
sizes <- cf$sizes; E <- t(.multilpa_error_matrix(pc)[,pc$modal])
ll <- function(th){ m <- matrix(th,3,2,byrow=TRUE); pr <- t(apply(cbind(m,0),1,function(v) exp(v)/sum(exp(v)))); sum(log(rowSums(E*tcrossprod(ind,pr)*matrix(sizes,n,3,byrow=TRUE)))) }
th0 <- as.vector(t(log(cf$probabilities[,1:2]/cf$probabilities[,3])))
cat("categorical info rel:",rel(.multilpa_step_categorical_information(cf,ind),-numDeriv::hessian(ll,th0)),"; score sum:",max(abs(colSums(.multilpa_step_categorical_scores(cf,ind)))),"; numeric grad of ll at estimate",max(abs(numDeriv::grad(ll,th0))),"\n")
# Categorical sandwich SE by hand against the table
Inv <- solve(-numDeriv::hessian(ll,th0)); S <- numDeriv::jacobian(function(th){ m <- matrix(th,3,2,byrow=TRUE); pr <- t(apply(cbind(m,0),1,function(v) exp(v)/sum(exp(v)))); log(rowSums(E*tcrossprod(ind,pr)*matrix(sizes,n,3,byrow=TRUE))) },th0)
Vt <- Inv%*%(crossprod(S)*n/(n-1))%*%Inv
Jg <- numDeriv::jacobian(function(th){ m <- matrix(th,3,2,byrow=TRUE); as.vector(t(t(apply(cbind(m,0),1,function(v) exp(v)/sum(exp(v)))))) },th0)
se_h <- sqrt(diag(Jg%*%Vt%*%t(Jg))); tab <- three_step(f,d,"cat",method="ml")
cat("categorical ML SE rel:",rel(tab$standard_error,se_h),"\n")
# Overall Wald on logits by hand, categorical ML
L <- rbind(c(-1,0,1,0,0,0),c(0,-1,0,1,0,0),c(-1,0,0,0,1,0),c(0,-1,0,0,0,1))
th_mat <- th0; d0 <- L%*%th_mat; W <- as.numeric(t(d0)%*%solve(L%*%Vt%*%t(L))%*%d0)
print(c(hand=W, pkg=three_step(f,d,"cat",method="ml",contrast="overall")$statistic, df=three_step(f,d,"cat",method="ml",contrast="overall")$df))
# r3step pooled: Hessian with step 1e-4 vs numDeriv in orthonormal coordinates; separated covariate near the boundary
des <- cbind(`(Intercept)`=1,x=d$x); basis <- .multilpa_design_basis(des); Z <- t(backsolve(basis,t(des),transpose=TRUE))
ba <- t(.multilpa_error_matrix(pc))[pc$modal,]
obj <- function(th){ pr <- .multilpa_logit_prior(matrix(th,2,2),Z); -sum(log(rowSums(pr*ba))) }
gr <- function(th){ pr <- .multilpa_logit_prior(matrix(th,2,2),Z); j <- pr*ba; po <- j/rowSums(j); r <- (po-pr)[,1:2]; -colSums(cbind(Z*r[,1],Z*r[,2])) }
th <- optim(rep(0,4),obj,gr,method="BFGS",control=list(reltol=1e-12,maxit=1000))
cat("r3step |gradient at BFGS stop|:",max(abs(gr(th$par))),"\n")
hh <- .multilpa_observed_hessian(function(s) obj(th$par+s),function(s) gr(th$par+s),rep(1,4),1e-4)
cat("r3step Hessian vs numDeriv rel:",rel(hh$scaled,numDeriv::hessian(obj,th$par)),"\n")
# Step-one J via .multilpa_step_one_correction's step vs numDeriv at unit changes
d2 <- transform(d,a=a*1e-4,b=b*1e4)
f2 <- lpa(d2,c("a","b"),n_profiles=3,n_starts=5,seed=1,tol=1e-10)
a1 <- three_step(f,d,"y",step_one="estimated"); a2 <- three_step(f2,d2,"y",step_one="estimated")
cat("step-one SE after rescaling indicators (matched by sorted estimate):",format(sort(a1$standard_error),digits=6)," | ",format(sort(a2$standard_error),digits=6),"\n")
# Zero V1: the correction must vanish exactly
ns <- asNamespace("latents"); orig <- get(".multilpa_step_one_correction",ns)
src <- paste(deparse(orig,width.cutoff=500L),collapse="\n")
zero <- eval(parse(text=sub("root <- .multilpa_covariance_root(covariance)","root <- .multilpa_covariance_root(0 * covariance)",src,fixed=TRUE)),envir=ns); environment(zero) <- environment(orig)
unlockBinding(".multilpa_step_one_correction",ns); assign(".multilpa_step_one_correction",zero,ns)
invisible(lapply(c("bch","proportional","ml"),function(m) cat(m,"max |SE(est, V1=0) - SE(fixed)|:",max(abs(three_step(f,d,"y",method=m,step_one="estimated")$standard_error-three_step(f,d,"y",method=m)$standard_error)),
   " categorical:",max(abs(three_step(f,d,"cat",method=m,step_one="estimated")$standard_error-three_step(f,d,"cat",method=m)$standard_error)),"\n")))
assign(".multilpa_step_one_correction",orig,ns); lockBinding(".multilpa_step_one_correction",ns)
