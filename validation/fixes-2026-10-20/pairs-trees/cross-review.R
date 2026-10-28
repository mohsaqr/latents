pkgload::load_all('.', quiet=TRUE)
set.seed(20261020)
n <- 600
cl <- rep(1:2, each=n/2)
d <- data.frame(g=rep(1:60,each=10), z=rnorm(n), w=rnorm(n))
d$y <- c(-3,3)[cl] + .4*d$z + rnorm(n)
d$a <- rbinom(n,1,c(.25,.75)[cl]); d$b <- rbinom(n,1,c(.3,.7)[cl])
d$o <- ordered(sample(1:3,n,TRUE)); d$c <- rpois(n,c(1.8,4)[cl])
fit <- function(data) multilpa(data,c('y','a','b','o','c'),'g',2,1,
  categorical=c('a','b'),ordinal='o',count='c',local_dependence='a:b',
  direct_effects='y ~ z',n_starts=2,seed=17,tol=1e-12)
a <- fit(d); b <- fit(transform(d,z=z+2026))
cat('Mixed pair/direct/ordinal/count LL difference:',a$log_likelihood-b$log_likelihood,'\n')
invisible(lapply(c('observed','opg','robust'),function(type) {
  ia <- parameter_inference(a,vcov_type=type)
  ib <- parameter_inference(b,vcov_type=type)
  select <- ia$parameter=='direct_effect'
  cat(type,' slope SEs:',ia$standard_error[select],ib$standard_error[select],'\n')
  stopifnot(max(abs(ia$standard_error[select]-ib$standard_error[select]))<1e-7)
}))
# Scores and information under a mixed model must occupy one consistent chart.
parts <- .latents_vuong_parts(a)
stopifnot(ncol(parts$scores)==ncol(parts$information),
          identical(dim(parts$information),rep(ncol(parts$scores),2)),
          parts$n_units==60)
cat('Mixed conditioned chart dimensions:',dim(parts$scores),dim(parts$information),'\n')
# Provenance in rbind subsets and id renaming; sources are regenerated.
pairs <- .latents_parse_pairs('a:b')
joined <- .latents_join_pairs(d,pairs)
combined <- rbind(joined[1:100,],joined[101:600,])
stopifnot(identical(.latents_join_pairs(combined,pairs)[['a:b']],joined[['a:b']]))
names(combined)[names(combined)=='g'] <- 'school'
stopifnot(identical(.latents_join_pairs(combined,pairs)[['a:b']],joined[['a:b']]))
cat('rbind/subset/id-renaming provenance: PASS\n')
# Arbitrary nonsingular linear changes of numerical helper coordinates.
set.seed(18)
scores <- matrix(rnorm(120*7),120,7)
j <- matrix(rnorm(7*4),7,4)
i1 <- crossprod(scores); s0 <- scores%*%j; i0 <- crossprod(s0)
a <- list(scores=scores,information=i1); b <- list(scores=s0,information=i0)
q1 <- diag(7);q1[1,5] <- 40;q1[2,7] <- -100
q0 <- diag(4);q0[1,3] <- 25
transform_parts <- function(x,q) list(scores=x$scores%*%q,information=t(q)%*%x$information%*%q)
ref <- .latents_vuong_weights(a,b)
changed <- .latents_vuong_weights(transform_parts(a,q1),transform_parts(b,q0))
cat('Vuong shear-transform maximum eigenvalue error:',max(abs(sort(ref)-sort(changed))),'\n')
stopifnot(max(abs(sort(ref)-sort(changed)))<1e-7)
cat('CROSS REVIEW COMPLETE\n')
# Categorical and continuous local dependence coexist with ordinal/count
# blocks, and prediction/refitting preserve the prepared joint table.
d$y2 <- .6*d$y + rnorm(n)
pattern <- multilpa(d,c('y','y2','a','b','o','c'),'g',2,1,
  categorical=c('a','b'),ordinal='o',count='c',
  local_dependence=c('a:b','y:y2'),n_starts=2,seed=17,tol=1e-12)
frame <- .multilpa_resolve_data(pattern,NULL)
stopifnot(max(abs(predict(pattern,d,type='posterior')-
                   predict(pattern,frame,type='posterior')))<1e-12)
refit <- do.call(multilpa,c(list(data=frame,n_starts=2,seed=17,tol=1e-12),
                          .multilpa_refit_arguments(pattern)))
cat('Mixed continuous/categorical pattern+ordinal/count refit LL error:',
    refit$log_likelihood-pattern$log_likelihood,'\n')
stopifnot(abs(refit$log_likelihood-pattern$log_likelihood)<1e-6,
          identical(refit$covariance_pattern,pattern$covariance_pattern))
cat('INTEGRATION CROSS REVIEW COMPLETE\n')
