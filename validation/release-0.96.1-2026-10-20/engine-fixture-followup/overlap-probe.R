pkgload::load_all('.',quiet=TRUE)
paths<-as.matrix(expand.grid(rep(list(1:2),4)))
people<-rep(seq_len(nrow(paths)),each=10)
residual<-rep(qnorm((seq_len(10)-.5)/10),nrow(paths))
invisible(lapply(c(.5,.7,1),function(spread) {
  d<-data.frame(id=rep(seq_along(people),each=4),time=rep(seq_len(4),length(people)),
     y=.1*(as.vector(t(2*paths[people,]-3))+spread*rep(residual,each=4)))
  invisible(lapply(c('default','internal','blas'),function(mode) {
    options(matprod=mode)
    f<-lta(d,'y','id',2,time='time',transitions='occasion',n_starts=1,seed=1)
    a<-parameter_inference(f,step=1e-2);b<-parameter_inference(f,step=1e-4)
    cat('SPREAD',spread,'MODE',mode,'LL',f$log_likelihood,'boundary',f$boundary,
      'converged',f$converged,'transitions',range(plogis(f$transition_coefficients)),
      'initial',range(plogis(f$initial_coefficients)),
      'SEdiff',max(abs(a$standard_error-b$standard_error)),'\n')
  }))
}))
