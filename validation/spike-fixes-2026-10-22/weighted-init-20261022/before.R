suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
set.seed(803)
d <- data.frame(y=c(rnorm(50,-3,.5),rnorm(50,3,.5)),w=1)
q <- function(x) suppressWarnings(suppressMessages(x))
a <- q(lpa(rbind(d,data.frame(y=0,w=0)),'y',2,weights='w',n_starts=10,seed=1))
b <- q(lpa(rbind(d,data.frame(y=1e4,w=0)),'y',2,weights='w',n_starts=10,seed=1))
print(c(plain=a$log_likelihood,outlier=b$log_likelihood));str(b$means)
plain <- q(lpa(d,'y',2,n_starts=3,seed=4)); rng <- .Random.seed
set.seed(91); initial <- .multilpa_initialize(matrix(d$y,ncol=1),rep(1:20,each=5),2,2,'varying',1e-6,3,hierarchical=FALSE); init_rng <- .Random.seed
saveRDS(list(data=d,fit=plain,rng=rng,initial=initial,init_rng=init_rng),'/private/tmp/latents-weighted-init-20261022/before.rds')
