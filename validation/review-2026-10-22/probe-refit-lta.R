# Read-only public API probes of HEAD 12e0145; run from the repository root.
# Third probe reproduces the bootstrap's EXACT ARGUMENT RECONSTRUCTION;
# it refits the original data, rather than performing a stochastic bootstrap.
Sys.setenv(VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
pkgload::load_all('.', quiet=TRUE)
with_notices <- function(expr) {
  notices <- character()
  fit <- withCallingHandlers(expr, warning=function(w) {
    notices <<- c(notices, paste(class(w)[1], conditionMessage(w), sep=': '))
    invokeRestart('muffleWarning')
  })
  list(fit=fit, notices=notices)
}
cat('Probe 1: all-missing subjects should not change thin-profile warning\n')
set.seed(1)
d <- data.frame(g=rep(1:30,each=3),t=rep(1:3,30),
                a=c(rnorm(87),rnorm(3,20)),b=c(rnorm(87),rnorm(3,20)))
args <- list(vars=c('a','b'),id='g',n_profiles=2,time='t',missing='fiml',
             n_starts=1,seed=1,max_iter=500)
x <- with_notices(do.call(lta,c(list(data=d),args)))
e <- data.frame(g=rep(31:130,each=3),t=rep(1:3,100),a=NA_real_,b=NA_real_)
y <- with_notices(do.call(lta,c(list(data=rbind(d,e)),args)))
invisible(lapply(list(original=x, with_empty=y), function(pair) {
  cat('LL: ',format(pair$fit$log_likelihood,digits=16),'\n',sep='')
  cat('Effective counts: ',paste(pair$fit$effective_profile_counts,collapse=', '),'\n',sep='')
  cat('Raw sums: ',paste(colSums(pair$fit$subject_posteriors),collapse=', '),'\n',sep='')
  cat('Warnings: ',paste(pair$notices,collapse=' | '),'\n',sep='')
}))
cat('Expected: both should warn for the same 3-member profile (4 own parameters).\n')
stopifnot(isTRUE(all.equal(x$fit$log_likelihood,y$fit$log_likelihood,tolerance=1e-8)),
          isTRUE(all.equal(x$fit$effective_profile_counts,y$fit$effective_profile_counts,tolerance=1e-8)))

cat('\nProbe 2: occasion-specific measurement should assess each separate block\n')
set.seed(1)
d2 <- data.frame(g=rep(1:30,each=3),t=rep(1:3,30),
                 a=c(rnorm(84),rnorm(6,20)),b=c(rnorm(84),rnorm(6,20)))
z <- with_notices(lta(d2,c('a','b'),'g',2,time='t',measurement='occasion',
                       n_starts=1,seed=1,max_iter=200))
cat('Pooled counts: ',paste(colSums(z$fit$subject_posteriors),collapse=', '),'\n',sep='')
cat('Memberships by occasion:\n')
print(rowsum(z$fit$subject_posteriors,d2$t))
cat('Measurement blocks: ',length(z$fit$measurement),'\n',sep='')
cat('Warnings: ',paste(z$notices,collapse=' | '),'\n',sep='')
cat('Expected: rare profile warns; each block has 2 effective members for 4 own parameters.\n')

cat('\nProbe 3: exact staged-bootstrap argument reconstruction, without resampling\n')
set.seed(6)
d3 <- data.frame(g=rep(1:12,each=4),a=rnorm(48),b=rnorm(48))
staged <- with_notices(fit_staged(d3,c('a','b'),'g',2,2,n_starts=1,n_screen=8,
                                 search='emEM',seed=1,max_iter=20))$fit
args <- latents:::.multilpa_staged_refit_arguments(staged)
control <- latents:::.latents_refit_controls(staged,n_starts=1)
# Mirrors R/bootstrap-inference.R:134-137, with resampled := original d3.
refit <- with_notices(do.call(fit_staged,c(list(data=d3),args,
  list(n_starts=1,max_iter=20,tol=1e-8,n_screen=control$n_screen))))$fit
cat('Original first-stage pool: ',staged$stage_one$n_screen,
    '; screened rows: ',nrow(staged$stage_one$starts),'\n',sep='')
cat('Returned second-stage pool: ',if(is.null(staged$n_screen)) 'NULL' else staged$n_screen,'\n',sep='')
cat('Inherited bootstrap pool: ',if(is.null(control$n_screen)) 'NULL' else control$n_screen,'\n',sep='')
cat('Reconstructed first-stage pool: ',refit$stage_one$n_screen,
    '; screened rows: ',nrow(refit$stage_one$starts),'\n',sep='')
cat('Expected: reconstructed first-stage pool remains 8.\n')
