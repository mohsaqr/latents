Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
library(latents,lib.loc='/private/tmp/latents-release-audit-20261028/lib')
set.seed(312)
g<-60L;m<-8L;h<-rep(0:1,length.out=g);z<-rbinom(g*m,1,rep(ifelse(h==0,.15,.85),each=m))
d<-data.frame(id=rep(seq_len(g),each=m),x=rnorm(g*m,5*z,.8),y=rnorm(g*m,-4*z,1))
fit <- function(scale, family) {
 dd<-d;dd[c('x','y')]<-dd[c('x','y')]*scale
 args<-list(data=dd,vars=c('x','y'),id='id',family=family,n_group_classes=2,n_starts=3,seed=2,tol=1e-10)
 if(family %in% c('profile','full_cross_level','restricted_cross_level')) args$n_profiles<-2
 tryCatch(suppressWarnings(do.call(multilpa,args)),error=identity)
}
results<-lapply(c('profile','additive','dispersion','additive_dispersion','full_cross_level','restricted_cross_level'),function(family) {
 a<-fit(1,family); b<-fit(1e-6,family)
 inf<-lapply(list(a,b),function(x) tryCatch(suppressWarnings(parameter_inference(x)),error=identity))
 list(family=family,base=a,small=b,inference=inf)
})
saveRDS(results,'validation/release-audit-2026-10-28/core/multilevel-units.rds')
invisible(lapply(results,function(x) {
 cat('\n',x$family,'\n')
 invisible(lapply(c('base','small'),function(n) {
  a<-x[[n]];cat(n,if(inherits(a,'error')) conditionMessage(a) else paste('logLik',a$log_likelihood,'converged',a$converged),'\n')
 }))
 invisible(lapply(x$inference,function(i) cat(if(inherits(i,'error')) conditionMessage(i) else paste('inference rows',nrow(i)),'\n')))
}))
