Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
pkgload::load_all('.', quiet = TRUE)
set.seed(271)
d <- data.frame(id=rep(seq_len(30), each=4), x=rep(rnorm(30),each=4)+rnorm(120), w=1)
print(str(d))
probe <- function(family, magnitude) {
  e <- data.frame(id=rep(31:32,each=4), x=magnitude+rep(c(-1,1),4), w=0)
  dd <- rbind(d,e)
  run <- function(data) tryCatch(withCallingHandlers(do.call(multilpa, c(list(data=data, vars='x', id='id', n_group_classes=2, family=family, weights='w', n_starts=2, seed=3), if (family != 'additive') list(n_profiles=2))), warning=function(w) invokeRestart('muffleWarning')), error=identity)
  a <- run(d); b <- run(dd)
  data.frame(family, magnitude, base=if(inherits(a,'error')) conditionMessage(a) else a$log_likelihood, extra=if(inherits(b,'error')) conditionMessage(b) else b$log_likelihood)
}
print(do.call(rbind,lapply(c('profile','additive','full_cross_level'), function(f) do.call(rbind,lapply(c(1e3,1e6,1e9), function(m) probe(f,m))))))
