pkgload::load_all('.',quiet=TRUE)
active_solver <- .multilpa_bounded_variances
expressions <- parse('validation/release-0.96.1-2026-10-20/covariance-followup/baseline-source.R')
declaration <- Filter(function(z) is.call(z) && identical(z[[1]], as.name('<-')) && identical(z[[2]],as.name('.multilpa_bounded_variances')),as.list(expressions))[[1]]
old_solver <- eval(declaration[[3]],asNamespace('latents'))
set.seed(92)
cases<-expand.grid(k=c(2,3),d=c(2,3,4),structure=c('EVI','VEI'),trial=1:5,
                   stringsAsFactors=FALSE)
results<-lapply(seq_len(nrow(cases)),function(i) {
  z<-cases[i,];w<-runif(z$k,10,100);s<-matrix(exp(rnorm(z$k*z$d,0,4)),z$k,z$d)*w
  floor<-.1;v<-active_solver(s,w,z$structure,floor)
  old<-tryCatch(old_solver(s,w,z$structure,floor),error=identity)
  cost<-function(v)sum(w*rowSums(log(v))+rowSums(s/v))/sum(w)
  stopifnot(min(v)>=floor*(1-1e-10))
  lv<-rowMeans(log(v));shape<-log(v)-lv
  if(z$structure=='EVI')stopifnot(max(lv)-min(lv)<1e-8)
  if(z$structure=='VEI')stopifnot(max(abs(sweep(shape,2,shape[1,])))<1e-8)
  gap<-if(inherits(old,'condition'))NA_real_ else cost(v)-cost(old)
  if(is.finite(gap))stopifnot(gap<1e-7)
  data.frame(z,cost_gap=gap)
})
r<-do.call(rbind,results);print(summary(r$cost_gap));cat('Cases:',nrow(r),'\n')
# Analytic active-bound EVI reference, k=d=2, weights=100 each.
s<-matrix(c(.0023,115,8,9.6),2,byrow=TRUE);w<-c(100,100);floor<-.05
a<-sqrt(prod(s[2,]));volume<-(a+sqrt(a*a+4*sum(w)*s[1,2]*floor))/(2*sum(w))
reference<-rbind(c(floor,volume^2/floor),volume*s[2,]/a)
v<-active_solver(s,w,'EVI',floor)
print(list(actual=v,reference=reference,error=max(abs(v-reference))))
stopifnot(max(abs(v-reference))<1e-9)
cat('REFERENCE CHECKS COMPLETE\n')
