pkgload::load_all('.',quiet=TRUE)
# The old course fixture under alternate matrix multiplication implementations.
invisible(lapply(c('default','internal','blas'),function(mode) {
  options(matprod=mode)
  warnings<-character()
  fit<-withCallingHandlers(lta(course_engagement,c('browse','lectures','forum_read'),
    'student',2,time='sequence',transitions='occasion',n_starts=1,seed=1),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  cat('OLD',mode,'LL',fit$log_likelihood,'boundary',fit$boundary,'prob_boundary',fit$probability_boundary,
      'converged',fit$converged,'warnings',length(warnings),'\n')
}))
# Every four-occasion binary state path occurs equally often. Each path has
# the same symmetric within-state residual quantiles, avoiding sample accidents.
paths<-as.matrix(expand.grid(rep(list(1:2),4)))
people<-rep(seq_len(nrow(paths)),each=10)
residual<-rep(stats::qnorm((seq_len(10)-.5)/10),nrow(paths))
d<-data.frame(id=rep(seq_along(people),each=4),time=rep(seq_len(4),length(people)),
              y=as.vector(t(6*paths[people,]-9))+.5*rep(residual,each=4))
saveRDS(d,'validation/release-0.96.1-2026-10-20/engine-fixture-followup/data.rds')
invisible(lapply(c('default','internal','blas'),function(mode) {
  options(matprod=mode)
  fit<-lta(d,'y','id',2,time='time',transitions='occasion',n_starts=1,seed=1)
  a<-parameter_inference(fit,step=1e-2);b<-parameter_inference(fit,step=1e-4)
  cat('NEW',mode,'LL',fit$log_likelihood,'boundary',fit$boundary,'prob_boundary',fit$probability_boundary,
      'converged',fit$converged,'SE max difference',max(abs(a$standard_error-b$standard_error)),'\n')
  print(names(fit));print(fit$params$transition)
  print(tryCatch(get_tna(fit),error=function(e)class(e)[1]))
}))
