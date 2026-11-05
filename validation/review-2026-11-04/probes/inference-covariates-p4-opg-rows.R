suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(G=50,n=6,seed=1){ set.seed(seed)
  do.call(rbind,lapply(seq_len(G),function(g){ h<-rbinom(1,1,.5); z<-rnorm(n); k<-rbinom(n,1,plogis(-.5+1.5*h+.8*z))+1L
    data.frame(g=g,z=z,wt=rep(runif(1,.5,2),n),
      y1=rnorm(n,c(-1.5,1.5)[k],1),y2=rnorm(n,c(1,-1)[k],1),y3=rnorm(n,c(0,1)[k],.7),
      cat=factor(sapply(k,function(kk) sample(c("p","q","r"),1,prob=list(c(.7,.2,.1),c(.15,.25,.6))[[kk]]))),
      ord=sapply(k,function(kk) sample(1:3,1,prob=list(c(.6,.3,.1),c(.1,.3,.6))[[kk]])),
      cnt=rpois(n,c(1,4)[k])) }))}
d<-sim(); dm<-d; set.seed(2); dm$y1[sample(nrow(dm),30)]<-NA; dm$y2[sample(nrow(dm),30)]<-NA
## 1. OPG under grand transforms, covariate-free path
a<-list(data=d,vars=c("y1","y2","y3"),id="g",n_profiles=2,n_group_classes=2,n_starts=3,seed=1,tol=1e-10,max_iter=5000)
f0<-do.call(multilpa,a); fc<-do.call(multilpa,c(a,list(centering="grand"))); fs<-do.call(multilpa,c(a,list(scaling="grand",centering="grand")))
cmp<-function(type){p0<-parameter_inference(f0,vcov_type=type);pc<-parameter_inference(fc,vcov_type=type);ps<-parameter_inference(fs,vcov_type=type)
 keep<-p0$parameter=="probability"
 data.frame(outcome=p0$outcome[keep],term=p0$term[keep],se0=p0$standard_error[keep],se_centred=pc$standard_error[keep],se_both=ps$standard_error[keep])}
cat("mixing probabilities must be invariant to a grand transform of the indicators\n")
print(cbind(type="observed",cmp("observed")),digits=4);print(cbind(type="robust",cmp("robust")),digits=4);print(cbind(type="opg",cmp("opg")),digits=4)
## what the correct OPG would be: influence built from the inverse information; compare mean SEs
pc_obs<-parameter_inference(fc,vcov_type="observed"); pc_opg<-parameter_inference(fc,vcov_type="opg"); p0_opg<-parameter_inference(f0,vcov_type="opg")
print(data.frame(parameter=pc_opg$parameter,outcome=pc_opg$outcome,term=pc_opg$term,obs=pc_obs$standard_error,opg=pc_opg$standard_error,opg_uncentred=p0_opg$standard_error)[1:6,],digits=4)
## 2. per-group scores: row by row against numDeriv of the group log likelihood
grp_check<-function(name,args){ out<-tryCatch({
  fit<-withCallingHandlers(do.call(multilpa,c(list(id="g",n_profiles=2,n_starts=2,seed=1,tol=1e-9,max_iter=3000),args)),warning=function(w)invokeRestart("muffleWarning"))
  th<-unname(.multilpa_coefficients(fit,"unconstrained"))
  prepared<-.multilpa_inference_matrix(fit,.multilpa_resolve_data(fit,NULL)); x<-prepared$x; codes<-prepared$codes
  obj<-fit; obj$means<-sweep(fit$means,2,prepared$centers,"-"); th<-unname(.multilpa_coefficients(obj,"unconstrained"))
  set.seed(3); t1<-th+rnorm(length(th),0,.02)
  gl<-function(t){e<-.multilpa_expectation(x,fit$group_index,.multilpa_decode(t,fit),codes,unname(fit$sampling_weights),fit$extra_data); w<-fit$sampling_weights%||%rep(1,fit$n_groups); e$group_log_likelihood*w}
  J<-numDeriv::jacobian(gl,t1); S<-.multilpa_group_scores(t1,x,obj,codes)
  data.frame(cfg=name,npar=fit$n_parameters,rows=nrow(S),cols=ncol(S),max_row_err=max(abs(J-S)),max_abs=max(abs(J)))
 },error=function(e)data.frame(cfg=name,npar=NA,rows=NA,cols=NA,max_row_err=NA,max_abs=NA,err=conditionMessage(e)))
 if(is.null(out$err)) out$err<-NA; out}
cfg<-list(diag=list(data=d,vars=c("y1","y2","y3")),full=list(data=d,vars=c("y1","y2","y3"),covariance_model="full"),
 fiml_full=list(data=dm,vars=c("y1","y2","y3"),covariance_model="full",missing="fiml"),
 fiml_diag_wt=list(data=dm,vars=c("y1","y2","y3"),missing="fiml",weights="wt"),
 mixed=list(data=d,vars=c("y1","cat","ord","cnt"),categorical="cat",ordinal="ord",count="cnt"),
 VEI=list(data=d,vars=c("y1","y2","y3"),model="VEI"),VVE=list(data=d,vars=c("y1","y2","y3"),model="VVE",n_group_classes=1),
 ld=list(data=d,vars=c("y1","y2","y3"),local_dependence="y1:y2"),
 direct=list(data=d,vars=c("y1","y2","y3"),direct_effects="y1~z"),
 direct_full=list(data=d,vars=c("y1","y2","y3"),direct_effects="y1~z",covariance_model="full"),
 direct_chart=list(data=d,vars=c("y1","y2","y3"),direct_effects="y1~z",model="VEI"))
print(do.call(rbind,Map(grp_check,names(cfg),cfg)),digits=3,row.names=FALSE)
