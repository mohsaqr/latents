Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
library(latents,lib.loc='/private/tmp/latents-release-audit-20261028/lib')
set.seed(531)
d<-data.frame(x=rnorm(500,2,1.5), y=rnorm(500,-3,.7), w=exp(rnorm(500)))
d$x[sample.int(500,60)]<-NA_real_;d$y[sample.int(500,75)]<-NA_real_;d$w[1:10]<-0
print(str(d));print(summary(d))
f<-suppressWarnings(lpa(d,c('x','y'),1,weights='w',missing='fiml',n_starts=1,tol=1e-12))
print(f$means);print(f$variances)
w<-f$sampling_weights
ref<-lapply(c('x','y'),function(v) {
 rows<-!is.na(d[[v]]) & w>0;x<-d[[v]][rows];ww<-w[rows]
 mu<-weighted.mean(x,ww);variance<-sum(ww*(x-mu)^2)/sum(ww)
 c(mean=mu,variance=variance,loglik=sum(ww*dnorm(x,mu,sqrt(variance),log=TRUE)))
})
print(ref)
stopifnot(max(abs(unname(f$means)-vapply(ref,`[[`,numeric(1),'mean')))<1e-10,
 max(abs(unname(f$variances)-vapply(ref,`[[`,numeric(1),'variance')))<1e-10,
 abs(f$log_likelihood-sum(vapply(ref,`[[`,numeric(1),'loglik')))<1e-9)
cat('PASS: weighted one-class FIML means, ML variances and log likelihood match analytical values.\n')
