Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet=TRUE)
capture <- function(expr) {
  warnings <- list()
  fit <- withCallingHandlers(expr, warning=function(w) {
    warnings[[length(warnings)+1L]] <<- w
    invokeRestart('muffleWarning')
  })
  list(fit=fit, warnings=warnings)
}
set.seed(1)
data <- data.frame(g=rep(1:30, each=3), t=rep(c(2,4,9),30),
 a=c(rnorm(87),rnorm(3,20)), b=c(rnorm(87),rnorm(3,20)))
empty <- data.frame(g=rep(31:130,each=3),t=rep(c(2,4,9),100),a=NA_real_,b=NA_real_)
fits <- lapply(c('homogeneous','occasion'), function(transitions) {
 lapply(c(TRUE,FALSE), function(drop_empty) capture(lta(rbind(data,empty),c('a','b'),'g',2,time='t', missing='fiml', transitions=transitions,drop_empty=drop_empty,n_starts=1,seed=1,max_iter=500)))
})
str(fits[[1]][[1]]$fit[c('informative','effective_profile_counts','subject_posteriors')])
print(lapply(fits,function(blocks) lapply(blocks,function(x) vapply(x$warnings,conditionMessage,character(1)))))
saveRDS(fits,'/private/tmp/lta-thin-before.rds')
