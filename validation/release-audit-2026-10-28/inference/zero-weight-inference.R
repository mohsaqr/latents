Sys.setenv(OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
pkgload::load_all('.',quiet=TRUE)
set.seed(72)
d <- data.frame(x=rnorm(100),w=1)
d$y <- 2+1.5*d$x+rnorm(100)
str(d)
f <- mixture_regression(y~x,d,1,weights='w',n_starts=1,seed=2)
print(as.data.frame(get_results(f,'coefficients')))
results <- lapply(c(1e2,1e4,1e5,1e6,1e7,1e8),function(v){
 dp <- rbind(d,data.frame(x=v,w=0,y=v))
 g <- mixture_regression(y~x,dp,1,weights='w',n_starts=1,seed=2)
 cat('zero row x=',v,'LLgap=',g$log_likelihood-f$log_likelihood,'\n')
 print(as.data.frame(get_results(g,'coefficients')))
 list(x=v,fit=g)
})
saveRDS(list(base=f,results=results),'validation/release-audit-2026-10-28/inference/zero-weight-inference.rds')
g <- results[[length(results)]]$fit
g$spec <- latents:::.mixture_carrying_spec(g$spec)
fixed <- latents:::.mixture_inference(g,'robust')
cat('Recomputing only inference on positive-weight units restores covariance, max absolute error:', max(abs(fixed$vcov-f$inference$vcov)),'\n')
stopifnot(max(abs(fixed$vcov-f$inference$vcov))<1e-10)

comparison <- do.call(rbind,lapply(results,function(r) {
 ref <- as.data.frame(get_results(f,'coefficients'))
 got <- as.data.frame(get_results(r$fit,'coefficients'))
 data.frame(zero_predictor=r$x,ll_gap=r$fit$log_likelihood-f$log_likelihood,
 coefficient_gap=max(abs(got$estimate-ref$estimate)),
 intercept_se=got$std_error[1],slope_se=got$std_error[2],
 max_relative_se_error=max(abs(got$std_error/ref$std_error-1)))
}))
print(comparison,digits=12,row.names=FALSE)
write.csv(comparison,'validation/release-audit-2026-10-28/inference/zero-weight-inference.csv',row.names=FALSE)
