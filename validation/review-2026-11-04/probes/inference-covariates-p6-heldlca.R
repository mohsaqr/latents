suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
## two-step LCA with a covariate: weak separation, one-level
n <- 500; z <- rnorm(n); k <- rbinom(n,1,plogis(.2+.8*z))+1L
q <- c(.75,.3)
d <- data.frame(g=seq_len(n),z=z)
for (v in paste0("u",1:5)) d[[v]] <- factor(ifelse(runif(n)<q[k],"yes","no"))
vars <- paste0("u",1:5)
base <- multilpa(d, vars, "g", n_profiles=2, n_group_classes=1, categorical=vars, n_starts=5, seed=1, tol=1e-10, max_iter=5000)
st <- starting_values(base, what="measurement")
two <- multilpa(d, vars, "g", n_profiles=2, n_group_classes=1, categorical=vars, profile_covariates="z", start=st, fixed="measurement", n_starts=1, seed=1, tol=1e-10, max_iter=5000)
cat("n_parameters",two$n_parameters,"fixed",two$fixed,"\n")
p <- parameter_inference(two); print(p[,c("term","estimate","standard_error")])
rp <- two$response_probabilities
ld <- sapply(1:2,function(kk) Reduce(`+`,lapply(vars,function(v) log(rp[[v]][kk,as.character(d[[v]])]))))
f <- function(b){e<-b[1]+b[2]*d$z;P<-plogis(e);sum(log(P*exp(ld[,1])+(1-P)*exp(ld[,2])))}
cf<-as.vector(two$profile_coefficients); cat("hand LL",f(cf),"pkg",two$log_likelihood,"\n")
V<-solve(-numDeriv::hessian(f,cf)); cat("conditional SE (hand):",sqrt(diag(V)),"\n pkg/hand:",p$standard_error/sqrt(diag(V)),"\n")
one<-multilpa(d, vars, "g", n_profiles=2, n_group_classes=1, categorical=vars, profile_covariates="z", n_starts=5, seed=1, tol=1e-10, max_iter=5000)
po<-parameter_inference(one); print(po[po$parameter=="coefficient",c("term","estimate","standard_error")])
## ordinal / count blocks under fixed = "measurement"
d$o<-pmin(3,1+rbinom(n,2,c(.3,.7)[k])); d$cn<-rpois(n,c(1,3)[k])
b2<-multilpa(d,c(vars,"o","cn"),"g",n_profiles=2,n_group_classes=1,categorical=vars,ordinal="o",count="cn",n_starts=3,seed=1)
st2<-starting_values(b2,what="measurement"); cat("start blocks:",names(st2),"\n")
f2<-tryCatch(multilpa(d,c(vars,"o","cn"),"g",n_profiles=2,n_group_classes=1,categorical=vars,ordinal="o",count="cn",profile_covariates="z",start=st2,fixed="measurement",n_starts=1,seed=1),error=function(e)conditionMessage(e))
if(is.character(f2)) cat("ordinal+fixed:",f2,"\n") else {cat("ordinal+fixed: n_par",f2$n_parameters,"\n");print(parameter_inference(f2)[,c("parameter","outcome","term","estimate","standard_error")]);print(coef(b2)[grepl("threshold|count",names(coef(b2)))])}
## catpair bootstrap
d3<-do.call(rbind,lapply(seq_len(50),function(g){nn<-6;kk<-rbinom(nn,1,.5)+1L;data.frame(g=g,y1=rnorm(nn,c(-1.5,1.5)[kk]),a=factor(sapply(kk,function(j)sample(c("p","q"),1,prob=list(c(.8,.2),c(.2,.8))[[j]]))),b=factor(sapply(kk,function(j)sample(c("r","s"),1,prob=list(c(.8,.2),c(.2,.8))[[j]]))))}))
x<-multilpa(d3,c("y1","a","b"),"g",n_profiles=2,n_group_classes=1,categorical=c("a","b"),local_dependence="a:b",n_starts=3,seed=1)
cat("catpair vars:",x$vars," n_par",x$n_parameters,"\n")
bb<-tryCatch(withCallingHandlers(parameter_inference(x,method="bootstrap",iter=8,n_starts=2,seed=2),warning=function(w){message("W: ",conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)conditionMessage(e))
if(is.character(bb)) cat("catpair bootstrap:",bb,"\n") else {cat("rows",nrow(bb)," valid",attr(bb,"n_valid"),"\n"); print(as.data.frame(bb)[,c("outcome","term","parameter","estimate","standard_error")],digits=3)}
## covariate LCA + VEI + categorical chart path
d3$z<-rnorm(nrow(d3)); d3$y2<-d3$y1+rnorm(nrow(d3)); 
y<-withCallingHandlers(multilpa(d3,c("y1","y2","a","b"),"g",n_profiles=2,n_group_classes=1,categorical=c("a","b"),model="VEI",profile_covariates="z",n_starts=3,seed=1),warning=function(w)invokeRestart("muffleWarning"))
pi<-parameter_inference(y); print(pi[,c("outcome","term","parameter","estimate","standard_error")],digits=4)
