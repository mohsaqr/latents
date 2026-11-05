suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(G=50,n=6,seed=1){ set.seed(seed)
  do.call(rbind,lapply(seq_len(G),function(g){ h<-rbinom(1,1,.5); z<-rnorm(n); k<-rbinom(n,1,plogis(-.5+1.5*h+.8*z))+1L
    data.frame(g=g,z=z,wt=rep(runif(1,.5,2),n),
      y1=rnorm(n,c(-1.5,1.5)[k],1),y2=rnorm(n,c(1,-1)[k],1),y3=rnorm(n,c(0,1)[k],.7),
      cat=factor(sapply(k,function(kk) sample(c("p","q","r"),1,prob=list(c(.7,.2,.1),c(.15,.25,.6))[[kk]]))),
      c2=factor(sapply(k,function(kk) sample(c("u","v"),1,prob=list(c(.7,.3),c(.2,.8))[[kk]]))),
      ord=sapply(k,function(kk) sample(1:3,1,prob=list(c(.6,.3,.1),c(.1,.3,.6))[[kk]])),
      cnt=rpois(n,c(1,4)[k])) }))}
d<-sim(); dm<-d; set.seed(2); dm$y1[sample(nrow(dm),30)]<-NA
## 1. refit fidelity: same data, same arguments -> same likelihood and count
cfg<-list(diag=list(data=d,vars=c("y1","y2","y3")),VEI=list(data=d,vars=c("y1","y2","y3"),model="VEI"),
 EVE=list(data=d,vars=c("y1","y2","y3"),model="EVE"),
 fiml=list(data=dm,vars=c("y1","y2","y3"),missing="fiml"),wt=list(data=d,vars=c("y1","y2","y3"),weights="wt"),
 mixed=list(data=d,vars=c("y1","cat","ord","cnt"),categorical="cat",ordinal="ord",count="cnt"),
 nb=list(data=d,vars=c("y1","cnt"),count="cnt",count_model="negative_binomial",count_dispersion="equal"),
 ld=list(data=d,vars=c("y1","y2","y3"),local_dependence="y1:y2"),
 catpair=list(data=d,vars=c("y1","cat","c2"),categorical=c("cat","c2"),local_dependence="cat:c2"),
 direct=list(data=d,vars=c("y1","y2","y3"),direct_effects="y1~z"),
 grand=list(data=d,vars=c("y1","y2","y3"),centering="grand",scaling="grand"),
 person=list(data=d,vars=c("y1","y2","y3"),centering="person"))
one<-function(nm,args){ out<-tryCatch({
  x<-withCallingHandlers(do.call(multilpa,c(list(id="g",n_profiles=2,n_starts=3,seed=1,tol=1e-9,max_iter=3000),args)),warning=function(w)invokeRestart("muffleWarning"))
  a<-.multilpa_refit_arguments(x)
  r<-withCallingHandlers(do.call(multilpa,c(list(data=args$data),a,list(n_starts=3,seed=1,tol=1e-9,max_iter=3000))),warning=function(w)invokeRestart("muffleWarning"))
  r2<-.multilpa_align_labels(r,x)
  b<-withCallingHandlers(parameter_inference(x,method="bootstrap",iter=8,n_starts=2,seed=2),warning=function(w)invokeRestart("muffleWarning"))
  data.frame(cfg=nm,ll=x$log_likelihood,ll_diff=r$log_likelihood-x$log_likelihood,np=x$n_parameters,np_refit=r$n_parameters,coef_diff=max(abs(coef(r2)-coef(x))),boot_rows=nrow(b),valid=attr(b,"n_valid"),struct=paste(x$covariance_structure,r$covariance_structure))
 },error=function(e)data.frame(cfg=nm,ll=NA,ll_diff=NA,np=NA,np_refit=NA,coef_diff=NA,boot_rows=NA,valid=NA,struct=conditionMessage(e)))
 out}
print(do.call(rbind,Map(one,names(cfg),cfg)),digits=5,row.names=FALSE)
## 2. bootstrap of a staged fit with a measurement from ANOTHER sample
other<-sim(seed=99); other$y1<-other$y1+3
meas<-multilpa(other,c("y1","y2","y3"),"g",n_profiles=2,n_group_classes=1,n_starts=3,seed=1)
st<-fit_staged(d,c("y1","y2","y3"),"g",n_profiles=2,n_group_classes=2,measurement=meas,n_starts=2,seed=1)
b<-withCallingHandlers(parameter_inference(st,method="bootstrap",iter=20,n_starts=2,seed=1),warning=function(w)invokeRestart("muffleWarning"))
print(as.data.frame(b)[1:6,c("outcome","term","parameter","estimate","standard_error","conf_low","conf_high")],digits=4)
cat("rows where the estimate is outside its own percentile interval:",sum(b$estimate<b$conf_low|b$estimate>b$conf_high),"of",nrow(b),"\n")
## 3. pooling by hand, with a forced label flip through the internal align
imp<-lapply(1:3,function(i){set.seed(100+i); x<-d; x$y1<-x$y1+rnorm(nrow(x),0,.15); x$z<-x$z+rnorm(nrow(x),0,.1); x})
pl<-pool_imputations(imp,c("y1","y2","y3"),"g",n_profiles=2,n_group_classes=2,profile_covariates="z",n_starts=3,seed=1,tol=1e-10,max_iter=4000)
print(get_results(pl,"fits"))
fits<-lapply(imp,function(x) multilpa(x,c("y1","y2","y3"),"g",n_profiles=2,n_group_classes=2,profile_covariates="z",n_starts=3,seed=1,tol=1e-10,max_iter=4000))
al<-lapply(fits,function(f).multilpa_pool_align(f,fits[[1]])$fit); tb<-lapply(al,parameter_inference)
Q<-sapply(tb,`[[`,"estimate");U<-sapply(tb,function(t)t$standard_error^2)
qb<-rowMeans(Q);ub<-rowMeans(U);b_<-apply(Q,1,var);tt<-ub+(1+1/3)*b_;r<-(1+1/3)*b_/ub;df<-2*(1+1/r)^2
e<-as.data.frame(pl)
cat("pool vs hand: est",max(abs(e$estimate-qb))," se",max(abs(e$standard_error-sqrt(tt)))," df rel",max(abs(e$df/df-1))," ci",max(abs(e$conf_low-(qb-qt(.975,df)*sqrt(tt)))),"\n")
ptest<-!is.na(e$p_value); cat("p",max(abs(e$p_value[ptest]-2*pt(-abs(qb[ptest]/sqrt(tt[ptest])),df[ptest]))),"\n")
ref<-fits[[1]]; flip<-.multilpa_cov_permute(ref,c(2L,1L),c(2L,1L)); back<-.multilpa_pool_align(flip,ref)
cat("flip orders:",back$profile_order,"/",back$group_order," coef diff:",max(abs(coef(back$fit)-coef(ref))),"\n")
v<-vcov(pl); cat("pooled vcov diag vs se^2:",max(abs(sqrt(diag(v))-e$standard_error)),"\n")
## 4. weights + grand centring: membership SE invariant?
a<-list(data=d,vars=c("y1","y2","y3"),id="g",n_profiles=2,n_group_classes=2,profile_covariates="z",n_starts=3,seed=1,tol=1e-10,max_iter=4000,weights="wt")
w0<-do.call(multilpa,a); wc<-do.call(multilpa,c(a,list(centering="grand")));ws<-do.call(multilpa,c(a,list(centering="grand",scaling="grand")))
p0<-parameter_inference(w0);pc<-parameter_inference(wc);ps<-parameter_inference(ws); k<-p0$parameter=="coefficient"
cat("weighted, coefficient SE ratio centred:",range(pc$standard_error[k]/p0$standard_error[k])," scaled:",range(ps$standard_error[k]/p0$standard_error[k]),"\n")
cat("grand offset:",wc$centering_offsets[1,],"\n weighted mean:",sapply(c("y1","y2","y3"),function(v)weighted.mean(d[[v]],d$wt)),"\n unweighted:",colMeans(d[c("y1","y2","y3")]),"\n")
