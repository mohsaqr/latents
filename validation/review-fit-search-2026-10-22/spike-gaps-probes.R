Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet = TRUE)
cap <- function(expr) {
  ws<-list()
  f<-withCallingHandlers(expr,warning=function(w){ws[[length(ws)+1L]]<<-list(class=class(w)[1L],message=conditionMessage(w));invokeRestart('muffleWarning')},message=function(m)invokeRestart('muffleMessage'))
  list(fit=f,warnings=ws)
}
cat('\nKNOWN MODERATE FLAT SPIKE\n')
a<-cap(lpa(srl,names(srl),4,model='VVV',n_starts=10,seed=1))
print(a$warnings);print(a$fit[c('converged','boundary','log_likelihood','effective_profile_counts')])
source('equivalence/search/spike-measure.R')
print(spike_measure(a$fit,apply(srl,2,var)*(nrow(srl)-1)/nrow(srl)))
cat('\nCROSS LEVEL BETWEEN DENSITY THINNESS\n')
set.seed(19);nn<-120;d<-data.frame(g=rep(1:40,each=3),a=rnorm(nn),b=rnorm(nn));d[d$g<=3,c('a','b')]<-d[d$g<=3,c('a','b')]+20
ans<-lapply(c('full_cross_level','restricted_cross_level'),function(fam){
  r<-cap(multilpa(d,c('a','b'),'g',n_profiles=1,n_group_classes=2,family=fam,between_variance='varying',n_starts=1,seed=1,max_iter=300))
  print(fam);print(r$warnings);print(r$fit[c('converged','boundary','group_mean_variances')]);print(colSums(r$fit$group_posteriors));r
})
saveRDS(list(known_flat=a,between=ans),'validation/review-fit-search-2026-10-22/spike-gaps-probes.rds')
