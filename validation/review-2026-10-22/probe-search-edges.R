Sys.setenv(VECLIB_MAXIMUM_THREADS = 1, OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
pkgload::load_all('.', quiet = TRUE)
ns <- asNamespace('latents')
run_starts <- get('.latents_run_starts',ns)
set.seed(32)
d <- data.frame(id=rep(1:20, each=4), a=rnorm(80),b=rnorm(80),z=rnorm(80), item=sample(letters[1:3],80,TRUE))
d$a[1:30] <- d$a[1:30] + 3
d$b[1:30] <- d$b[1:30] + 3
cases <- list(diagonal=list(), full=list(covariance_model='full'), pattern=list(local_dependence='a:b'), direct=list(direct_effects='a ~ z'), prior=list(prior=prior_control()), noise=list(noise=TRUE), mixed=list(vars=c('a','b','item'),categorical='item'), missing=list(missing='fiml'))
results <- lapply(names(cases),function(nm) {
  args <- utils::modifyList(list(data=d,vars=c('a','b'),id='id',n_profiles=2,n_group_classes=1,n_starts=3,max_iter=55,tol=1e-13,seed=9),cases[[nm]])
  if(nm=='missing') args$data$a[c(1,2,15)] <- NA
  if(nm=='noise') {args['id'] <- list(NULL)}
  a <- suppressWarnings(do.call(multilpa,c(args,list(search='kmeans'))))
  b <- suppressWarnings(do.call(multilpa,c(args,list(search='emEM',n_screen=6))))
  both <- subset(b$starts,stage=='final' & start<=3)
  data.frame(case=nm, n_shared=nrow(both), exact_loglik=identical(both$log_likelihood,a$starts$log_likelihood[both$start]),exact_iterations=identical(both$iterations,a$starts$iterations[both$start]))
})
print(do.call(rbind,results))
# Screening includes a failure during prepare, a failure during resume, a protected boundary start,
# and an interior survivor. These should never make a screened-only attempt eligible.
prepare <- function(i) {if(i==2L) stop('prepare failed');i}
run <- function(i,values,pause=Inf,resume=NULL) {
  if(!is.null(resume) && i==3L) stop('resume failed')
  list(state=i,parameters=list(i=i),expectation=list(log_likelihood=-i),iterations=if(is.null(resume)) 2L else 10L,converged=!is.null(resume),checkpoint=list(state=i))
}
f <- suppressWarnings(run_starts(3,run,'converged',prepare=prepare,interior=function(a)a$parameters$i!=1L,screen=list(n_pool=5,iterations=2L,protect=1L)))
print(f[c('valid','stage','best_start')])
# Evaluate-only preserves supplied start and ignores the screening pool deliberately.
g <- suppressWarnings(multilpa(d,c('a','b'),NULL,2,n_starts=3,seed=9))
h <- suppressWarnings(multilpa(d,c('a','b'),NULL,2,n_starts=3,max_iter=0,start=starting_values(g),search='emEM',n_screen=12,seed=4))
print(h[c('n_starts','n_screen','iterations','search')])
print(all.equal(g$means,h$means,tolerance=1e-12))
