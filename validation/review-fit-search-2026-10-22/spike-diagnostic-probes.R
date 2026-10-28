Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet = TRUE)
cap <- function(expr) {
  ws <- list()
  f <- withCallingHandlers(expr, warning=function(w) {ws[[length(ws)+1L]] <<- list(class=class(w)[1L],message=conditionMessage(w));invokeRestart('muffleWarning')}, message=function(m) invokeRestart('muffleMessage'))
  list(fit=f, warnings=ws)
}
show <- function(z, gaussian=NULL) {
  print(z$warnings)
  print(z$fit[c('converged','boundary','log_likelihood','effective_profile_counts')])
  if(!is.null(gaussian)) print(colSums(z$fit$subject_posteriors[gaussian,,drop=FALSE]))
}
set.seed(73)
d <- as.data.frame(matrix(rnorm(80*5),80,5));names(d)<-letters[1:5]
d[1:3,]<-d[1:3,]+10
start <- list(means=rbind(rep(10,5),rep(0,5)),variances=matrix(1,2,5),profile_probabilities=matrix(c(3,77)/80,1),group_probabilities=1)
f <- function(frame,...) cap(lpa(frame,letters[1:5],2,start=start,n_starts=1,seed=1,missing='fiml',...))
cat('\nORDINARY EMPTY ROWS\n')
a<-f(d,drop_empty=FALSE);show(a)
empty<-d[rep(1L,300L),];empty[,]<-NA_real_
b<-f(rbind(d,empty),drop_empty=FALSE);show(b,c(rep(TRUE,80),rep(FALSE,300)))
cat('\nORDINARY CATEGORICAL-ONLY ROWS\n')
d$categ<-c(rep('rare',3),rep('common',77));empty$categ<-'rare'
mix <- rbind(d,empty)
stcat <- start; stcat$response_probabilities<-list(categ=rbind(c(1e-10,1-1e-10),c(1-1e-10,1e-10)))
c<-cap(lpa(mix,c(letters[1:5],'categ'),2,categorical='categ',start=stcat,n_starts=1,seed=1,missing='fiml'));show(c,c(rep(TRUE,80),rep(FALSE,300)))
cat('\nLTA EIGENVALUE BOUNDARY\n')
set.seed(12);nn<-180L; cls<-sample(c(-5,5),nn,replace=TRUE);z<-rnorm(nn)
lt<-data.frame(id=rep(1:60,each=3),t=rep(1:3,60),a=cls+z,b=cls+z+rnorm(nn,sd=1e-5))
results<-lapply(c('homogeneous','occasion'),function(tr) {
  ans<-cap(lta(lt,c('a','b'),'id',2,time='t',model='VVV',transitions=tr,n_starts=1,seed=1,max_iter=300,min_variance=1e-6))
  print(tr);show(ans)
  blocks<-if(is.null(ans$fit[['measurement']])) list(ans$fit) else ans$fit[['measurement']]
  print(lapply(blocks,function(block) list(diagonals=block$variances,eigenvalues=apply(block$covariances,3,function(s) eigen(s,symmetric=TRUE,only.values=TRUE)$values))))
  ans
})
saveRDS(list(empty_base=a,empty_appended=b,categorical=c,lta=results),'validation/review-fit-search-2026-10-22/spike-diagnostic-probes.rds')
