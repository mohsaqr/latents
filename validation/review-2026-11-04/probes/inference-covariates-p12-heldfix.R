suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
n <- 500; z <- rnorm(n); k <- rbinom(n,1,plogis(.2+.8*z))+1L; q <- c(.75,.3)
d <- data.frame(g=seq_len(n),z=z); for (v in paste0("u",1:5)) d[[v]] <- factor(ifelse(runif(n)<q[k],"yes","no")); vars <- paste0("u",1:5)
base <- multilpa(d, vars, "g", n_profiles=2, n_group_classes=1, categorical=vars, n_starts=5, seed=1, tol=1e-10, max_iter=5000)
two <- multilpa(d, vars, "g", n_profiles=2, n_group_classes=1, categorical=vars, profile_covariates="z", start=starting_values(base,what="measurement"), fixed="measurement", n_starts=1, seed=1, tol=1e-10, max_iter=5000)
rp <- two$response_probabilities
ld <- sapply(1:2,function(kk) Reduce(`+`,lapply(vars,function(v) log(rp[[v]][kk,as.character(d[[v]])]))))
f <- function(b){e<-b[1]+b[2]*d$z;P<-plogis(e);sum(log(P*exp(ld[,1])+(1-P)*exp(ld[,2])))}
cf<-as.vector(two$profile_coefficients); hand<-sqrt(diag(solve(-numDeriv::hessian(f,cf))))
cat("before: pkg",parameter_inference(two)$standard_error," hand",hand,"\n"); cat("vcov(unconstrained) dim before:",dim(vcov(two,scale="unconstrained")),"\n")
ns<-asNamespace("latents"); g<-get(".multilpa_cov_held_parameters",ns); src<-deparse(g); src<-sub('"response")','c("response", "response_logit"))',src,fixed=TRUE); h<-eval(parse(text=src),envir=ns); environment(h)<-ns; unlockBinding(".multilpa_cov_held_parameters",ns); assign(".multilpa_cov_held_parameters",h,ns)
.latents_cache<-get(".latents_cache",ns); rm(list=ls(.latents_cache),envir=.latents_cache)
cat("after one-line patch: pkg",parameter_inference(two)$standard_error," hand",hand,"\n"); cat("vcov(unconstrained) dim after:",dim(vcov(two,scale="unconstrained")),"\n")
