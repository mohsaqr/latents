Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
pkgload::load_all('.',quiet=TRUE)
set.seed(8)
n <- 80L
path <- unlist(lapply(seq_len(n), function(i) Reduce(function(previous,step) if(runif(1)<.85) previous else 3L-previous,seq_len(4),sample(1:2,1),accumulate=TRUE)))
d <- data.frame(id=rep(seq_len(n),each=5),time=rep(1:5,n),z=rnorm(n*5),y1=rnorm(n*5,c(0,2.5)[path]),y2=rnorm(n*5,c(0,1.5)[path],.8))
str(d)
rescale <- function(f,s){
 f$x <- f$x*s;f$centers<-f$centers*s;f$min_variance<-f$min_variance*s^2
 f$measurement<-lapply(f$measurement,function(b){b$means<-b$means*s;b$variances<-b$variances*s^2;if(!is.null(b$covariances))b$covariances<-b$covariances*s^2;b})
 f
}
results <- lapply(c('EEE','VVV','VEE','EVE'),function(model){
 cat('MODEL',model,'\n')
 f <- lta(d,c('y1','y2'),'id',2,time='time',transition_covariates='z',model=model,n_starts=1,seed=2,tol=1e-10)
 scales <- latents:::.lta_inference_scales(f)
 stopifnot(identical(names(scales),names(coef(f))),length(scales)==length(coef(f)),all(is.finite(scales)),all(scales>0))
 cat('coefficient count',length(coef(f)),'scale count',length(scales),'\n')
 print(scales)
 ref <- lapply(c('observed','robust','opg'),function(type)vcov(f,type=type));names(ref)<-c('observed','robust','opg')
 rows <- lapply(c(1e-5,1e6),function(s){
  g <- rescale(f,s)
  gs <- latents:::.lta_inference_scales(g)
  stopifnot(identical(names(gs),names(coef(g))),length(gs)==length(coef(g)),all(is.finite(gs)),all(gs>0))
  map <- ifelse(startsWith(names(coef(f)),'mean.') | startsWith(names(coef(f)),'cholesky.'),s,1)
  do.call(rbind,lapply(names(ref),function(type){
   v <- vcov(g,type=type)
   err <- max(abs(v/tcrossprod(map)-ref[[type]])/tcrossprod(sqrt(diag(ref[[type]]))))
   cat('model',model,'scale',s,'type',type,'max standardized covariance error',err,'\n')
   data.frame(model=model,scale=s,type=type,error=err)
  }))
 })
 do.call(rbind,rows)
})
results <- do.call(rbind,results)
print(results)
write.csv(results,'validation/release-fixes-2026-10-28/lta/covariance-charts.csv',row.names=FALSE)
stopifnot(all(results$error<1e-3))
