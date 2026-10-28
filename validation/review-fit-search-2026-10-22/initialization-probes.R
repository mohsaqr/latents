suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
quiet <- function(expr) withCallingHandlers(expr, warning=function(w) {cat('WARNING:',conditionMessage(w),'\n');invokeRestart('muffleWarning')},message=function(m)invokeRestart('muffleMessage'))
set.seed(803)
n <- 100L
d <- data.frame(y=c(rnorm(n/2,-3,.5),rnorm(n/2,3,.5)),w=1)
a <- quiet(lpa(d,'y',2,weights='w',n_starts=10,seed=1))
cat('BASE LL',a$log_likelihood,'MEANS',a$means,'VARS',a$variances,'\n')
for (extreme in c(1e4,1e10,1e16,1e20)) {
  dd <- rbind(d, data.frame(y=extreme,w=0))
  b <- tryCatch(quiet(lpa(dd,'y',2,weights='w',n_starts=10,seed=1)),error=function(e)e)
  if(inherits(b,'error')) cat('ZERO',extreme,'ERROR',conditionMessage(b),'\n') else {
    cat('ZERO',extreme,'LL',b$log_likelihood,'MEANS',b$means,'VARS',b$variances,'\n')
  }
}
# Search for violated nonempty-cluster guarantee in Ward's reassignment repair.
set.seed(20261022)
for (trial in 1:30000) {
  nr <- sample(6:40,1); nc <- sample(1:5,1); k <- sample(3:min(nr-1,12),1)
  x <- matrix(rnorm(nr*nc),nr,nc)
  if (trial %% 3L == 0) x <- matrix(rt(nr*nc,df=1),nr,nc)
  assigned <- latents:::.multilpa_ward_assignments(x,k)
  if (length(unique(assigned))<k) {
    saveRDS(list(trial=trial,x=x,k=k,assigned=assigned),'validation/review-fit-search-2026-10-22/ward-empty.rds')
    cat('WARD EMPTY',trial,nr,nc,k,'sizes',tabulate(assigned,k),'\n');break
  }
  if(trial %% 5000L == 0)cat('Ward checked',trial,'\n')
}
