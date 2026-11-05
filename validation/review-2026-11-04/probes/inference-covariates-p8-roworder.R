suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(G=40,n=6,seed=1){ set.seed(seed)
  do.call(rbind,lapply(seq_len(G),function(g){ w<-rnorm(1); h<-rbinom(1,1,plogis(.8*w)); z<-rnorm(n); k<-rbinom(n,1,plogis(-.8+2.2*h+.8*z))+1L
    data.frame(g=sprintf("u%03d",G+1-g),w=w,z=z,wt=rep(runif(1,.5,2),n),y1=rnorm(n,c(-1.8,1.8)[k],1),y2=rnorm(n,c(1.5,-1.5)[k],1)) }))}
d<-sim(); set.seed(9); shuffled<-d[sample(nrow(d)),]; rownames(shuffled)<-NULL
dm<-d; set.seed(2); dm$y1[sample(nrow(dm),25)]<-NA; sm<-dm[as.integer(rownames(shuffled)),]
set.seed(9); idx<-sample(nrow(d)); shuffled<-d[idx,]; sm<-dm[idx,]; rownames(shuffled)<-rownames(sm)<-NULL
cat("first ids original:",unique(d$g)[1:3]," shuffled:",unique(shuffled$g)[1:3],"\n")
cmp<-function(name,data1,data2,...){ a<-list(vars=c("y1","y2"),id="g",n_profiles=2,n_group_classes=2,n_starts=4,seed=1,tol=1e-10,max_iter=5000,...)
  f1<-do.call(multilpa,c(list(data=data1),a)); f2<-do.call(multilpa,c(list(data=data2),a)); f2a<-align_labels(f2,f1)
  out<-lapply(c("observed","robust","opg"),function(tp){ r<-tryCatch({p1<-parameter_inference(f1,vcov_type=tp);p2<-parameter_inference(f2a,vcov_type=tp)
     data.frame(cfg=name,type=tp,ll_diff=f1$log_likelihood-f2$log_likelihood,est=max(abs(p1$estimate-p2$estimate)),se_ratio_lo=min(p2$standard_error/p1$standard_error,na.rm=TRUE),se_ratio_hi=max(p2$standard_error/p1$standard_error,na.rm=TRUE))},
     error=function(e)data.frame(cfg=name,type=tp,ll_diff=NA,est=NA,se_ratio_lo=NA,se_ratio_hi=NA)); r})
  do.call(rbind,out)}
res<-rbind(cmp("plain",d,shuffled),cmp("full",d,shuffled,covariance_model="full"),cmp("cov",d,shuffled,profile_covariates="z",group_covariates="w"),
 cmp("cov_full",d,shuffled,profile_covariates="z",covariance_model="full"),cmp("cov_fiml_full",dm,sm,profile_covariates="z",covariance_model="full",missing="fiml"),
 cmp("plain_fiml",dm,sm,missing="fiml"),cmp("grand",d,shuffled,centering="grand"),cmp("grand_cov",d,shuffled,centering="grand",profile_covariates="z"),
 cmp("VEI",d,shuffled,model="VEI"),cmp("wt",d,shuffled,weights="wt"),cmp("person",d,shuffled,centering="person"))
print(res,digits=4,row.names=FALSE)
