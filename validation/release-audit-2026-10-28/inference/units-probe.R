Sys.setenv(OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1)
pkgload::load_all('.', quiet=TRUE)
set.seed(8)
n <- 80L
path <- unlist(lapply(seq_len(n), function(i) Reduce(function(previous, step) if(runif(1) < .85) previous else 3L-previous, seq_len(4), sample(1:2,1), accumulate=TRUE)))
d <- data.frame(id=rep(seq_len(n),each=5),time=rep(1:5,n),z=rnorm(n*5), y1=rnorm(n*5,c(0,2.5)[path]),y2=rnorm(n*5,c(0,1.5)[path],.8))
str(d)
fits <- lapply(c(1,1e-3,1e-5,1e3,1e5), function(s) {
 ds <- d; ds[c('y1','y2')] <- ds[c('y1','y2')]*s
 f <- lta(ds,c('y1','y2'),'id',2,time='time',transition_covariates='z',n_starts=1,seed=2,tol=1e-10)
 v <- tryCatch(vcov(f),error=function(e)e)
 cat('SCALE',s,'LL', f$log_likelihood,'CONVERGED',f$converged,'VCOV',class(v),'\n')
 if (inherits(v,'error')) print(v) else print(sqrt(diag(v)))
 list(scale=s,fit=f,vcov=v)
})
saveRDS(fits,'validation/release-audit-2026-10-28/inference/units-probe.rds')
