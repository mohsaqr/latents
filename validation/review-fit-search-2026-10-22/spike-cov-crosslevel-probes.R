Sys.setenv(VECLIB_MAXIMUM_THREADS=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
pkgload::load_all(quiet = TRUE)
cap <- function(expr) {
  ws <- list()
  f <- withCallingHandlers(expr,warning=function(w){ws[[length(ws)+1L]] <<- list(class=class(w)[1L],message=conditionMessage(w));invokeRestart('muffleWarning')},message=function(m)invokeRestart('muffleMessage'))
  list(fit=f,warnings=ws)
}
set.seed(73)
d <- as.data.frame(matrix(rnorm(80*5),80,5));names(d)<-letters[1:5];d[1:3,]<-d[1:3,]+10;d$z<-rnorm(80)
cat('\nCOVARIATES: GAUSSIAN-ONLY\n')
a<-cap(lpa(d,letters[1:5],2,profile_covariates='z',n_starts=1,seed=1,missing='fiml',max_iter=1000));print(a$warnings);print(a$fit[c('effective_profile_counts','converged','boundary')])
d$categ<-c(rep('rare',3),rep('common',77));empty<-d[rep(1L,300L),];empty[,letters[1:5]]<-NA_real_;empty$z<-rnorm(300)
cat('\nCOVARIATES: CATEGORICAL-ONLY ROWS\n')
b<-cap(lpa(rbind(d,empty),c(letters[1:5],'categ'),2,profile_covariates='z',categorical='categ',n_starts=1,seed=1,missing='fiml',max_iter=1000));print(b$warnings);print(b$fit[c('effective_profile_counts','converged','boundary')]);print(colSums(b$fit$subject_posteriors[1:80,]))
cat('\nFULL CROSS LEVEL: SCALED WEIGHT COUNTS\n')
set.seed(7);cross<-data.frame(g=rep(1:40,each=2),a=rnorm(80),b=rnorm(80));cross[1:6,c('a','b')]<-cross[1:6,c('a','b')]+20;cross$weight<-ifelse(cross$g<=3,.01,1)
c<-cap(multilpa(cross,c('a','b'),'g',2,1,family='full_cross_level',weights='weight',n_starts=1,seed=1,max_iter=300));print(c$warnings);print(c$fit[c('converged','boundary','sampling_weights')]);print(colSums(c$fit$subject_posteriors));print(colSums(c$fit$subject_posteriors*c$fit$sampling_weights[c$fit$sufficient_statistics$index]));
saveRDS(list(cov_base=a,cov_cat=b,cross_level=c),'validation/review-fit-search-2026-10-22/spike-cov-crosslevel-probes.rds')
