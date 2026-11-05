suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(G=50,n=6,seed=1){ set.seed(seed)
  do.call(rbind,lapply(seq_len(G),function(g){ w<-rnorm(1); h<-rbinom(1,1,plogis(.8*w)); z<-rnorm(n); k<-rbinom(n,1,plogis(-.8+2.2*h+.8*z))+1L
    data.frame(g=g,w=w,z=z,wt=rep(runif(1,.5,2),n),y1=rnorm(n,c(-1.8,1.8)[k],1),y2=rnorm(n,c(1.5,-1.5)[k],1)) }))}
d<-sim()
## 1. units and origin in the covariate path
a<-list(vars=c("y1","y2"),id="g",n_profiles=2,n_group_classes=2,profile_covariates="z",group_covariates="w",n_starts=4,seed=1,tol=1e-10,max_iter=5000)
f0<-do.call(multilpa,c(list(data=d),a))
d2<-transform(d,y1=y1*1e4,y2=y2/1e3,z=z*1000+5e4,w=w/100-7)
f1<-do.call(multilpa,c(list(data=d2),a)); f1a<-align_labels(f1,f0)
for(tp in c("observed","robust","opg")){ p0<-parameter_inference(f0,vcov_type=tp);p1<-parameter_inference(f1a,vcov_type=tp)
 k<-p0$term%in%c("z","w"); cat(tp,"LL shift",f1$log_likelihood-f0$log_likelihood," slope z-stat diff",max(abs(p0$statistic[k]-p1$statistic[k])),
  " y1 mean SE ratio/1e4 range",range(p1$standard_error[p0$parameter=="mean"&p0$term=="y1"]/p0$standard_error[p0$parameter=="mean"&p0$term=="y1"]/1e4),"\n")}
## 2. weighted covariate sandwich, by hand from per-group numeric scores
fw<-do.call(multilpa,c(list(data=d,weights="wt"),a))
th<-.multilpa_cov_encode(fw); x<-sweep(fw$indicator_data,2,fw$center,"-")
gl<-function(t){p<-.multilpa_cov_decode(t,fw);e<-.multilpa_cov_expectation(x,fw$group_index,p$parameters,fw$profile_design,fw$group_design,p$beta,p$gamma,NULL,unname(fw$sampling_weights),NULL);e$group_log_likelihood*unname(fw$sampling_weights)}
J<-numDeriv::jacobian(gl,th); H<-numDeriv::jacobian(function(t)colSums(.multilpa_cov_group_scores(t,x,fw,NULL)),th); H<-(H+t(H))/2
Ai<-solve(-H); G<-50; V<-Ai%*%crossprod(J)%*%Ai*G/(G-1); pw<-parameter_inference(fw); labs<-.multilpa_cov_labels(fw); Jn<-.multilpa_cov_natural_jacobian(fw,th,labs)
cat("weighted robust max rel diff:",max(abs(sqrt(diag(Jn%*%V%*%t(Jn)))/pw$standard_error-1)),"  vcov_type:",attr(pw,"vcov_type"),"\n")
## 3. structural constants in the pool
imp<-lapply(1:3,function(i){set.seed(100+i); x<-d; x$y1<-x$y1+rnorm(nrow(x),0,.1); x})
pl<-pool_imputations(imp,c("y1","y2"),"g",n_profiles=2,n_group_classes=1,n_starts=3,seed=1); print(as.data.frame(pl)[,c("outcome","term","parameter","estimate","standard_error","df","within","between","riv","fmi")],digits=4)
print(summary(pl)[,c("relabelled","fmi_max","fmi_median")])
## 4. control: staged bootstrap with the measurement estimated on the same data
st<-fit_staged(d,c("y1","y2"),"g",n_profiles=2,n_group_classes=2,n_starts=2,seed=1)
b<-withCallingHandlers(parameter_inference(st,method="bootstrap",iter=30,n_starts=2,seed=1),warning=function(w)invokeRestart("muffleWarning"))
cat("same-data staged: estimate outside its interval in",sum(b$estimate<b$conf_low|b$estimate>b$conf_high),"of",nrow(b),"rows; valid",attr(b,"n_valid"),"\n")
## 5. measurement from a sample of the same population (no shift)
other<-sim(seed=77); m<-multilpa(other,c("y1","y2"),"g",n_profiles=2,n_group_classes=1,n_starts=3,seed=1)
st2<-fit_staged(d,c("y1","y2"),"g",n_profiles=2,n_group_classes=2,measurement=m,n_starts=2,seed=1)
b2<-withCallingHandlers(parameter_inference(st2,method="bootstrap",iter=30,n_starts=2,seed=1),warning=function(w)invokeRestart("muffleWarning"))
print(as.data.frame(b2)[1:4,c("outcome","term","parameter","estimate","standard_error","conf_low","conf_high")],digits=4)
cat("other-sample staged: outside in",sum(b2$estimate<b2$conf_low|b2$estimate>b2$conf_high),"of",nrow(b2),"\n")
cat("first-stage identity check available? stage_one n obs",st2$stage_one$n_observations,"; identical indicator data:",identical(unname(st2$stage_one$indicator_data),unname(st2$indicator_data)),"\n")
