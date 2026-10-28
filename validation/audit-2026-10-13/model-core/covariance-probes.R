pkgload::load_all('.', quiet=TRUE)
ns <- asNamespace('latents')
f <- function(name) get(name, ns)
diagonals <- rbind(c(.001,10),c(1,1))
weights <- c(10,10)
solve_variance <- f('.multilpa_structure_variances')
results <- lapply(c('EVI','VEI'), function(code) {
  variances <- solve_variance(diagonals,weights,code,.1)
  list(code=code,variances=variances,determinants=apply(variances,1,prod), ratios=variances[1,]/variances[2,])
})
print(results)
set.seed(7)
d <- data.frame(x=c(rnorm(100, -10,.005),rnorm(100,10,.3)),
                y=c(rnorm(100,0,1),rnorm(100,0,.3)))
fits <- lapply(c('EVI','VEI','EVV','VEV'), function(code) {
  warnings <- character()
  tryCatch({
    fit <- withCallingHandlers(lpa(d,c('x','y'),2,model=code,min_variance=.05,n_starts=1,seed=1,tol=1e-10),
                               warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
    matrices <- if (!is.null(fit$covariances)) lapply(1:2,function(k) matrix(fit$covariances[,,k],2,2)) else
      lapply(1:2,function(k) diag(fit$variances[k,]))
    list(code=code,converged=fit$converged,variance=fit$variances, determinants=vapply(matrices,det,numeric(1)),
         eigenvalues=lapply(matrices,eigen,symmetric=TRUE,only.values=TRUE),warnings=warnings)
  },error=function(e) list(code=code,error=conditionMessage(e),class=class(e),warnings=warnings))
})
print(fits)
saveRDS(list(internal=results, public=fits),'validation/audit-2026-10-13/model-core/covariance-probes.rds')
