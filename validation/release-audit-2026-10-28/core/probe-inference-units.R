Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
library(latents, lib.loc='/private/tmp/latents-release-audit-20261028/lib')
set.seed(842)
n <- 450L
z <- rep(0:1,length.out=n)
d <- data.frame(x=rnorm(n,5*z,0.7+0.5*z), y=rnorm(n,-4*z,1.2-0.4*z))
print(str(d))
models <- c('EII','VII','EEI','VEI','EVI','VVI','EEE','VEE','EVE','VVE','EEV','VEV','EVV','VVV')
fits <- lapply(models,function(model) {
  out <- lapply(c(1,1e-6,1e6),function(scale) {
    fit <- suppressWarnings(lpa(d*scale,names(d),2,model=model,n_starts=2,seed=4,tol=1e-10))
    inf <- tryCatch(suppressWarnings(parameter_inference(fit)),error=identity)
    list(fit=fit, inf=inf,scale=scale)
  }); names(out)<-c('base','small','large'); out
}); names(fits)<-models
print(names(fits[[1]]$base$inf)); print(head(fits[[1]]$base$inf$table))
saveRDS(fits,'validation/release-audit-2026-10-28/core/inference-units.rds')
print(do.call(rbind,lapply(names(fits),function(model) do.call(rbind,lapply(fits[[model]],function(result) data.frame(model,scale=result$scale,converged=result$fit$converged,inference=if(inherits(result$inf,'error')) conditionMessage(result$inf) else 'ok'))))))
