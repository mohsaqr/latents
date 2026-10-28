pkgload::load_all('.',quiet=TRUE)
d<-readRDS('validation/release-0.96.1-2026-10-20/covariance-followup/data.rds')
set.seed(43)
orders<-lapply(1:30,function(i)sample(nrow(d)))
invisible(lapply(seq_along(orders),function(i) {
  invisible(lapply(c('EVI','VEI','EVV','VEV','VEE','EVE','VVE'),function(code) {
    data<-d[orders[[i]],]
    fit<-tryCatch(withCallingHandlers(lpa(data,c('x','y'),2,model=code,
      min_variance=.05,n_starts=1,seed=1,tol=1e-10),
      warning=function(w)invokeRestart('muffleWarning')),error=identity)
    if(inherits(fit,'condition')) {
      cat(i,code,conditionMessage(fit),'\n');flush.console()
      saveRDS(list(data=data,code=code,permutation=i),paste0(
        'validation/release-0.96.1-2026-10-20/covariance-followup/permutation-',i,'-',code,'.rds'))
    }
  }))
}))
cat('DONE\n')
