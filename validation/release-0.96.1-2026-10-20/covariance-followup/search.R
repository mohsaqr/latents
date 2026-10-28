pkgload::load_all('.',quiet=TRUE)
d <- readRDS('validation/release-0.96.1-2026-10-20/covariance-followup/data.rds')
cases <- expand.grid(code=c('EVI','VEI','EVV','VEV','VEE','EVE','VVE'),
                     digits=c(15,14,13,12,11,10,9,8),stringsAsFactors=FALSE)
results <- lapply(seq_len(nrow(cases)),function(i) {
  code <- cases$code[i];data <- signif(d,cases$digits[i])
  fit <- tryCatch(withCallingHandlers(lpa(data,c('x','y'),2,model=code,min_variance=.05,
    n_starts=1,seed=1,tol=1e-10),warning=function(w)invokeRestart('muffleWarning')),error=identity)
  if(inherits(fit,'condition')) {
    cat(code,cases$digits[i],conditionMessage(fit),'\n');flush.console()
    saveRDS(list(data=data,code=code,digits=cases$digits[i]),
      paste0('validation/release-0.96.1-2026-10-20/covariance-followup/failure-',i,'.rds'))
  }
  data.frame(cases[i,],failed=inherits(fit,'condition'))
})
print(do.call(rbind,results))
