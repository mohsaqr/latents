suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(G=50,n=6,seed=1){ set.seed(seed)
  do.call(rbind,lapply(seq_len(G),function(g){ w<-rnorm(1); h<-rbinom(1,1,plogis(.5*w)); z<-rnorm(n); k<-rbinom(n,1,plogis(-.5+1.5*h+.8*z))+1L
    data.frame(g=g,w=w,z=z,wt=rep(runif(1,.5,2),n),
      y1=rnorm(n,c(-1.5,1.5)[k],1),y2=rnorm(n,c(1,-1)[k],1)+.0,y3=rnorm(n,c(0,1)[k],.7),
      cat=factor(sapply(k,function(kk) sample(c("p","q","r"),1,prob=list(c(.7,.2,.1),c(.15,.25,.6))[[kk]]))),
      ord=sapply(k,function(kk) sample(1:3,1,prob=list(c(.6,.3,.1),c(.1,.3,.6))[[kk]])),
      cnt=rpois(n,c(1,4)[k])) }))}
d<-sim(); dm<-d; set.seed(2); dm$y1[sample(nrow(dm),30)]<-NA; dm$y2[sample(nrow(dm),30)]<-NA; dm$cat[sample(nrow(dm),20)]<-NA
cfgs<-list(
 full_complete=list(data=d,vars=c("y1","y2","y3"),n_group_classes=2,profile_covariates="z",group_covariates="w",covariance_model="full"),
 full_fiml=list(data=dm,vars=c("y1","y2","y3"),n_group_classes=2,profile_covariates="z",covariance_model="full",missing="fiml"),
 diag_fiml_cat=list(data=dm,vars=c("y1","y2","y3","cat"),categorical="cat",n_group_classes=2,profile_covariates="z",group_covariates="w",missing="fiml"),
 weighted_fiml_slopes=list(data=dm,vars=c("y1","y2","y3"),n_group_classes=2,profile_covariates="z",profile_slopes="group_class",weights="wt",missing="fiml"),
 ord_count=list(data=d,vars=c("y1","cat","ord","cnt"),categorical="cat",ordinal="ord",count="cnt",n_group_classes=1,profile_covariates="z"),
 nb=list(data=d,vars=c("y1","cnt"),count="cnt",count_model="negative_binomial",n_group_classes=1,profile_covariates="z"),
 VEI=list(data=d,vars=c("y1","y2","y3"),n_group_classes=2,profile_covariates="z",model="VEI"),
 EEV=list(data=d,vars=c("y1","y2","y3"),n_group_classes=1,profile_covariates="z",model="EEV"),
 VVE=list(data=d,vars=c("y1","y2","y3"),n_group_classes=1,profile_covariates="z",model="VVE"),
 equal_full=list(data=d,vars=c("y1","y2","y3"),n_group_classes=1,profile_covariates="z",variance_model="equal",covariance_model="full"),
 pair=list(data=transform(d,c2=factor(paste0("l",pmin(ord,3))),c3=factor(sample(c("u","v"),nrow(d),TRUE))),vars=c("y1","y2","c2","c3"),categorical=c("c2","c3"),local_dependence="c2:c3",n_group_classes=1,profile_covariates="z")
)
res<-lapply(names(cfgs),function(nm){ cfg<-cfgs[[nm]]
  out<-tryCatch({
    fit<-withCallingHandlers(do.call(multilpa,c(list(id="g",n_profiles=2,n_starts=2,seed=1,tol=1e-9,max_iter=3000),cfg)),warning=function(w)invokeRestart("muffleWarning"))
    th<-.multilpa_cov_encode(fit); x<-fit$indicator_data; if(ncol(x)>0) x<-sweep(x,2,fit$center,"-"); codes<-fit$categorical_data
    ll<-function(t){p<-.multilpa_cov_decode(t,fit);.multilpa_cov_expectation(x,fit$group_index,p$parameters,fit$profile_design,fit$group_design,p$beta,p$gamma,codes,unname(fit$sampling_weights),fit$extra_data)$log_likelihood}
    set.seed(3); t1<-th+rnorm(length(th),0,.03)
    an<-colSums(.multilpa_cov_group_scores(t1,x,fit,codes)); nu<-numDeriv::grad(ll,t1)
    gerr<-max(abs(an-nu)/(1+abs(nu)))
    labs<-.multilpa_cov_labels(fit)
    J<-.multilpa_cov_natural_jacobian(fit,th,labs); Jn<-numDeriv::jacobian(function(t).multilpa_cov_natural_estimate(fit,t,labs),th)
    jerr<-max(abs(J-Jn))
    est<-.multilpa_cov_natural_estimate(fit,th,labs); pi<-parameter_inference(fit)
    H<- numDeriv::jacobian(function(t) colSums(.multilpa_cov_group_scores(t,x,fit,codes)), th); V<-solve(-(H+t(H))/2)
    se<-sqrt(pmax(diag(J%*%V%*%t(J)),0)); seerr<-max(abs(se/pi$standard_error-1),na.rm=TRUE)
    data.frame(cfg=nm,npar=fit$n_parameters,ntheta=length(th),conv=fit$converged,grad_err=gerr,jac_err=jerr,se_err=seerr,stringsAsFactors=FALSE)
  },error=function(e) data.frame(cfg=nm,npar=NA,ntheta=NA,conv=NA,grad_err=NA,jac_err=NA,se_err=NA,err=conditionMessage(e)))
  if(is.null(out$err)) out$err<-NA
  out})
print(do.call(rbind,res),digits=3)
# ---- invariance of membership inference to grand centring and scaling
base_args<-list(data=d,vars=c("y1","y2","y3"),id="g",n_profiles=2,n_group_classes=2,profile_covariates="z",group_covariates="w",n_starts=3,seed=1,tol=1e-10,max_iter=5000)
f0<-do.call(multilpa,base_args); fc<-do.call(multilpa,c(base_args,list(centering="grand"))); fs<-do.call(multilpa,c(base_args,list(scaling="grand"))); fcs<-do.call(multilpa,c(base_args,list(centering="grand",scaling="grand")))
tab<-function(f,type) {p<-parameter_inference(f,vcov_type=type); p[p$parameter=="coefficient",c("outcome","term","estimate","standard_error")]}
for(type in c("observed","robust","opg")){ t0<-tab(f0,type);tc<-tab(fc,type);ts<-tab(fs,type);tcs<-tab(fcs,type)
 cat(type,": est diff",max(abs(t0$estimate-tcs$estimate))," SE ratio centred",range(tc$standard_error/t0$standard_error)," scaled",range(ts$standard_error/t0$standard_error)," both",range(tcs$standard_error/t0$standard_error),"\n")}
# mean SE under grand centring vs jackknife-free check: mean_k - xbar variance by robust delta
pc<-parameter_inference(fc,vcov_type="robust"); p0<-parameter_inference(f0,vcov_type="robust")
print(cbind(pc[pc$parameter=="mean",c("outcome","term","estimate","standard_error")],uncentred_se=p0$standard_error[p0$parameter=="mean"]))
# ---- multiplicity family: location shift of an indicator moves p_adjusted of a covariate slope
d2<-d; d2$y1<-d2$y1+1000; d2$y2<-d2$y2+1000; d2$y3<-d2$y3+1000
d3<-d; # centred so that means straddle zero
a1<-do.call(multilpa,modifyList(base_args,list(data=d)))
a2<-do.call(multilpa,modifyList(base_args,list(data=d2)))
for(m in c("BH","holm","bonferroni")){
 p1<-parameter_inference(a1,adjust=m); p2<-parameter_inference(a2,adjust=m)
 s1<-p1[p1$term%in%c("z","w"),c("outcome","term","p_value","p_adjusted")]; s2<-p2[p2$term%in%c("z","w"),c("p_value","p_adjusted")]
 cat(m,"\n"); print(cbind(s1,p_shifted=s2$p_value,adj_shifted=s2$p_adjusted),digits=4)
 cat(" family size:",sum(!is.na(p1$p_value)),"of",nrow(p1)," rows; mean tests in family:",sum(!is.na(p1$p_value)&p1$parameter=="mean"),"\n")
}
# ---- prior with covariates; noise with covariates
pr<-tryCatch(do.call(multilpa,c(base_args,list(prior=prior_control()))),error=function(e)conditionMessage(e)); if(is.character(pr)) cat("prior:",pr,"\n") else {cat("prior accepted\n"); print(utils::head(parameter_inference(pr)))}
