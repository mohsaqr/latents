Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet = TRUE)
noise_values <- c(.03,.08,.12)
out<-setNames(lapply(noise_values, function(noise) {
 set.seed(12);n<-180L;cls<-sample(c(-5,5),n,replace=TRUE);z<-rnorm(n)
 d<-data.frame(id=rep(1:60,each=3),t=rep(1:3,60),a=cls+z,b=cls+z+rnorm(n,sd=noise))
 warnings <- list()
 f<-withCallingHandlers(lta(d,c('a','b'),'id',2,time='t',model='VVV',transitions='occasion',n_starts=1,seed=1,max_iter=300,min_variance=.01), warning=function(condition) {
   warnings[[length(warnings)+1L]] <<- list(classes=class(condition), message=conditionMessage(condition))
   invokeRestart('muffleWarning')
 })
 eigenvalues <- apply(f$measurement[[1]]$covariances,3,function(s)eigen(s,symmetric=TRUE,only.values=TRUE)$values)
 cat('\nNOISE',noise,'\n');print(f[c('converged','boundary','probability_boundary')]);print(eigenvalues);print(warnings)
 inf<-tryCatch(parameter_inference(f),error=function(e)list(error_class=class(e),message=conditionMessage(e)))
 print(inf)
 stopifnot(isTRUE(f$converged), identical(f$boundary,FALSE), all(abs(eigenvalues[2L,]-.01)<1e-10), length(warnings)==0L, is.data.frame(inf), all(is.finite(inf$standard_error)), all(is.finite(inf$p_value)))
 reference <- f
 reference$boundary <- TRUE
 boundary_error <- tryCatch(parameter_inference(reference),error=identity)
 stopifnot(inherits(boundary_error,'latents_boundary_fit'))
 cat('ASSERTIONS PASSED; corrected boundary flag raises',class(boundary_error)[1L],'\n')
 list(fit=f,inference=inf,data=d,warnings=warnings,corrected_flag_error=boundary_error)
}), as.character(noise_values))
saveRDS(out,'validation/review-fit-search-2026-10-22/lta-boundary-wald-probe.rds')
