Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet=TRUE)
capture_fit <- function(...) {
  warnings <- list()
  fit <- tryCatch(withCallingHandlers(multilpa(...), warning=function(notice) {
    warnings[[length(warnings)+1L]] <<- list(class=class(notice)[1L], message=conditionMessage(notice))
    invokeRestart('muffleWarning')
  }), error=function(e) list(error=conditionMessage(e)))
  list(fit=fit,warnings=warnings)
}
set.seed(7)
weighted <- data.frame(g=rep(1:40,each=2), a=rnorm(80), b=rnorm(80))
weighted[1:6,c('a','b')] <- weighted[1:6,c('a','b')] + 20
weighted$w <- ifelse(weighted$g<=3,.01,1)
fits <- lapply(c('full_cross_level','restricted_cross_level'), function(family) capture_fit(weighted,c('a','b'),'g',2,1,family=family,weights='w',n_starts=1,seed=1,max_iter=300))
set.seed(19)
between <- data.frame(g=rep(1:40,each=3), a=rnorm(120), b=rnorm(120))
between[between$g<=3,c('a','b')] <- between[between$g<=3,c('a','b')] + 20
fits <- c(fits,lapply(c('full_cross_level','restricted_cross_level'), function(family) capture_fit(between,c('a','b'),'g',1,2,family=family,n_starts=1,seed=1,max_iter=300)))
lapply(fits, function(item) { if (!is.null(item$fit$error)) return(print(item$fit$error)); print(item$warnings); print(item$fit[c('family','converged','boundary')]); print(colSums(item$fit$subject_posteriors)); print(colSums(item$fit$group_posteriors)); print(dim(item$fit$profile_means)) })
saveRDS(fits,'/private/tmp/latents-cross-level-fixes-20261022/after.rds')
