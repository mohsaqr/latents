pkgload::load_all('/tmp/astra-core-release', quiet=TRUE)
set.seed(1027)
d<-data.frame(y=rep(2,60),x=seq(-2,3,length.out=60),o=rnorm(60))
invisible(lapply(list(y~0+x,y~1+offset(o)),function(f){
 ref<-lm(f,d)
 fit<-mixture_regression(f,d,1,n_starts=0,vcov_type='none')
 print(c(beta=coef(ref),sigma2=mean(residuals(ref)^2),refll=as.numeric(logLik(ref)),fitll=as.numeric(logLik(fit)),converged=fit$converged)); print(fit$params)
}))
