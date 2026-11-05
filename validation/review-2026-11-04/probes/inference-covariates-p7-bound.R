suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
## LCA, 2 classes, N=400, item u1 3 categories with class-1 never in cat "c" (non-last floor case impossible to force exactly,
## so category labels are arranged two ways: floored category in the middle, floored category last)
sim <- function(lev){ n<-400; z<-rnorm(n); k<-rbinom(n,1,plogis(.3+.7*z))+1L
  u1<-ifelse(k==1, sample(c("a","b"),n,TRUE,c(.6,.4)), sample(c("a","b","c"),n,TRUE,c(.2,.3,.5)))
  d<-data.frame(g=seq_len(n),z=z,u1=factor(u1,levels=lev))
  for(v in paste0("v",1:4)) d[[v]]<-factor(ifelse(runif(n)<c(.85,.2)[k],"y","n"))
  d}
hand <- function(fit,d,vars,cov=FALSE){
  th0 <- coef(fit,"unconstrained"); nm <- names(th0); bound<-.multilpa_response_bound(fit)
  free <- !(nm %in% bound$unconstrained)
  blocks <- fit$response_probabilities; items<-names(blocks)
  ll <- function(tf){ th<-th0; th[free]<-tf; p<-list(); at<-0L
    rp<-lapply(items,function(i){b<-blocks[[i]]; w<-ncol(b)-1L; mm<-matrix(th[grepl(paste0("response_logit.*\\.",i,":"),nm)],nrow(b),w,byrow=TRUE); e<-cbind(exp(mm),1); e/rowSums(e)})
    names(rp)<-items
    ld<-sapply(1:2,function(k) Reduce(`+`,lapply(items,function(i) log(rp[[i]][k,match(as.character(d[[i]]),colnames(blocks[[i]]))]))))
    if(cov){ cf<-th[grepl("^profile\\.coefficient",nm)]; P<-plogis(cf[1]+cf[2]*d$z)} else { lg<-th[grepl("^profile\\.",nm)]; P<-rep(plogis(lg[1]),nrow(d)) }
    sum(log(P*exp(ld[,1])+(1-P)*exp(ld[,2])))}
  list(ll=ll,free=free,th0=th0,V=tryCatch(solve(-numDeriv::hessian(ll,th0[free])),error=function(e)NULL))}
for(lev in list(c("a","b","c"),c("a","c","b"))){
  d<-sim(lev); set.seed(2); cat("\n=== levels",lev,"\n")
  v<-c("u1",paste0("v",1:4))
  fit<-multilpa(d,v,"g",n_profiles=2,n_group_classes=1,categorical=v,n_starts=5,seed=1,tol=1e-10,max_iter=5000)
  print(round(fit$response_probabilities$u1,8))
  e<-tryCatch(parameter_inference(fit),error=function(e)class(e)[1]); cat("boundary=error:",if(is.character(e)) e else "accepted","\n")
  pf<-parameter_inference(fit,boundary="fix"); h<-hand(fit,d,v); cat("hand LL",h$ll(h$th0[h$free]),"pkg",fit$log_likelihood,"\n")
  uc<-attr(pf,"covariance_unconstrained")
  se_h<-sqrt(diag(h$V)); se_p<-sqrt(diag(uc)); cat("held:",attr(pf,"fixed_at_bound"),"\n"); cat("max rel diff unconstrained SE (fit vs base-R numDeriv):",max(abs(se_p/se_h-1)),"\n")
  print(pf[grepl("u1",pf$term)|pf$parameter=="probability",c("outcome","term","parameter","estimate","standard_error")],digits=4)
  fc<-multilpa(d,v,"g",n_profiles=2,n_group_classes=1,categorical=v,profile_covariates="z",n_starts=5,seed=1,tol=1e-10,max_iter=5000)
  pc<-parameter_inference(fc,boundary="fix"); hc<-hand(fc,d,v,TRUE); cat("cov hand LL",hc$ll(hc$th0[hc$free]),"pkg",fc$log_likelihood,"\n")
  uc2<-vcov(fc,scale="unconstrained",boundary="fix"); keep<-diag(uc2)>0
  cat("cov max rel diff SE:",max(abs(sqrt(diag(uc2))[keep]/sqrt(diag(hc$V))-1)),"\n")
  print(pc[grepl("u1",pc$term),c("outcome","term","estimate","standard_error","conf_low","conf_high")],digits=4)
}
## covariate fit: Wald at a loose tol (no score check on this path)
sim2<-function(G=50,n=6){set.seed(5);do.call(rbind,lapply(seq_len(G),function(g){z<-rnorm(n);k<-rbinom(n,1,plogis(.8*z))+1L;data.frame(g=g,z=z,y1=rnorm(n,c(-1,1)[k]),y2=rnorm(n,c(1,-1)[k]))}))}
d2<-sim2()
for(tol in c(1e-8,1e-4,1e-2)){ f<-suppressWarnings(multilpa(d2,c("y1","y2"),"g",n_profiles=2,n_group_classes=1,profile_covariates="z",n_starts=2,seed=1,tol=tol))
 w<-character(); p<-withCallingHandlers(parameter_inference(f),warning=function(x){w<<-c(w,class(x)[1]);invokeRestart("muffleWarning")})
 th<-.multilpa_cov_encode(f); x<-sweep(f$indicator_data,2,f$center,"-"); g<-colSums(.multilpa_cov_group_scores(th,x,f,NULL)); V<-vcov(f,scale="unconstrained")
 cat("cov tol",tol,"converged",f$converged,"disp(%SE)",100*max(abs(g)*sqrt(diag(V))),"warnings:",w,"\n")
 f0<-suppressWarnings(multilpa(d2,c("y1","y2"),"g",n_profiles=2,n_group_classes=1,n_starts=2,seed=1,tol=tol)); w<-character()
 p<-withCallingHandlers(parameter_inference(f0),warning=function(x){w<<-c(w,class(x)[1]);invokeRestart("muffleWarning")}); cat("  plain tol",tol,"disp(%SE)",100*attr(p,"score_displacement"),"warnings:",w,"\n")}
## covariate + person centring runs; single level covariate
f<-multilpa(d2,c("y1","y2"),"g",n_profiles=2,n_group_classes=1,profile_covariates="z",centering="person",n_starts=2,seed=1); print(attr(parameter_inference(f),"vcov_type"))
## OPG with the influence built from the inverse information: restores invariance?
a<-list(data=d2,vars=c("y1","y2"),id="g",n_profiles=2,n_group_classes=1,n_starts=3,seed=1,tol=1e-10)
f0<-do.call(multilpa,a);fc<-do.call(multilpa,c(a,list(centering="grand")))
po0<-parameter_inference(f0,vcov_type="opg");poc<-parameter_inference(fc,vcov_type="opg");pr<-parameter_inference(fc,vcov_type="observed")
gs<-attr(poc,"group_scores"); Oinv<-solve(crossprod(gs)); Ainv<-solve(attr(pr,"hessian"))
dimnames(Oinv)<-dimnames(Ainv)<-list(colnames(gs),colnames(gs))
V<-.multilpa_transform_covariance(fc,Oinv,Ainv,seq_len(ncol(gs)),NULL,robust=FALSE)
J<-.multilpa_inference_jacobian(fc); se<-sqrt(diag(J%*%V%*%t(J)))
print(data.frame(parameter=poc$parameter,outcome=poc$outcome,term=poc$term,opg_uncentred=po0$standard_error,opg_centred_now=poc$standard_error,opg_centred_fixed=se,observed_centred=pr$standard_error),digits=4,row.names=FALSE)
