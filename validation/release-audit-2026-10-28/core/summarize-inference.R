f <- readRDS('validation/release-audit-2026-10-28/core/inference-units.rds')
print(f[[1]]$base$inf)
summaries <- do.call(rbind,lapply(names(f),function(model) {
 base<-f[[model]]$base$inf
 do.call(rbind,lapply(f[[model]][c('small','large')],function(q) {
  d<-q$inf; power<-ifelse(d$parameter=='mean',1,ifelse(d$parameter %in% c('variance','covariance'),2,0))
  stopifnot(identical(d$parameter,base$parameter),identical(d$outcome,base$outcome))
  d$estimate<-d$estimate/q$scale^power
  d$standard_error<-d$standard_error/q$scale^power
  data.frame(model,scale=q$scale,max_estimate_error_in_SE=max(abs(d$estimate-base$estimate)/base$standard_error,na.rm=TRUE),max_relative_SE_error=max(abs(d$standard_error/base$standard_error-1),na.rm=TRUE))
 }))
}))
print(summaries)
write.csv(summaries,'validation/release-audit-2026-10-28/core/inference-units-summary.csv',row.names=FALSE)
stopifnot(all(summaries$max_estimate_error_in_SE<0.001),all(summaries$max_relative_SE_error<0.001))
