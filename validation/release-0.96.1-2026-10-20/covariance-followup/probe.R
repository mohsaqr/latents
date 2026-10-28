pkgload::load_all('.', quiet=TRUE)
set.seed(7)
d <- data.frame(x=c(rnorm(100,-10,.005),rnorm(100,10,.3)),
                y=c(rnorm(100,0,1),rnorm(100,0,.3)))
saveRDS(d,'validation/release-0.96.1-2026-10-20/covariance-followup/data.rds')
invisible(lapply(c('EVI','VEI','EVV','VEV','VEE','EVE','VVE'),function(code) {
  started <- proc.time()[['elapsed']]
  fit <- tryCatch(withCallingHandlers(lpa(d,c('x','y'),2,model=code,min_variance=.05,
    n_starts=1,seed=1,tol=1e-10),warning=function(w)invokeRestart('muffleWarning')),
    error=identity)
  cat(code,class(fit)[1],if(inherits(fit,'condition')) conditionMessage(fit) else
    paste('LL',fit$log_likelihood,'iterations',fit$iterations,
      'minimum increment',min(diff(fit$log_likelihood_history))),
    'seconds',proc.time()[['elapsed']]-started,'\n')
}))
