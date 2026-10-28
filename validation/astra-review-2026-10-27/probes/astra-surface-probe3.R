pkgload::load_all('.', quiet=TRUE)
set.seed(1027)
d <- data.frame(id=rep(1:20,each=5),x=sample(0:1,100,TRUE),y=sample(0:1,100,TRUE))
calls <- 0L
make_data <- function() {calls <<- calls+1L;d}
f <- multilca(make_data(),vars=c('x','y'),id=NULL,n_profiles=1,n_starts=1,seed=1)
cat('MULTILCA DATA EVALUATIONS',calls,'\n')
calls <- 0L
f <- lca(make_data(),vars=c('x','y'),n_classes=1,n_starts=1,seed=1)
cat('LCA DATA EVALUATIONS',calls,'\n')
# A one-class regression should explicitly say no classification exists.
r <- mixture_regression(y~1,d,n_classes=1,n_starts=1,seed=1)
print(get_results(diagnostics(r,plots=FALSE),'checks'))
# Emulate a saved pre-check diagnostics object, as the package test does.
x <- diagnostics(f,plots=FALSE);x$checks<-NULL;x$replication<-NULL
cat('LEGACY DIAGNOSTICS\n'); print(x);print(identical(as.data.frame(x),get_results(x,'entropy')))
