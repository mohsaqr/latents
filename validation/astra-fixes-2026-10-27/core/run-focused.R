Sys.setenv(NOT_CRAN='true', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
pkgload::load_all('.',quiet=TRUE)
library(testthat)
files <- c('astra-core', 'categorical', 'missing-covariance', 'local-dependence',
 'covariate-inference', 'covariate-missing', 'predict', 'get-data',
 'mixture-regression', 'growth-mixture', 'growth-correlated', 'growth-predict',
 'review-1026-codex', 'review-1026-growth-r2', 'review-1026-lta-r2',
 'review-1026-lta-open', 'lta-bootstrap', 'lta-inference')
files <- paste0('tests/testthat/test-',files,'.R')
files <- files[file.exists(files)]
results <- lapply(files,function(path){cat('\nFILE ',path,'\n');r<-test_file(path, reporter='summary',package='latents', stop_on_failure=FALSE);saveRDS(r,file.path('validation/astra-fixes-2026-10-27/core',paste0(basename(path),'.rds')));d<-as.data.frame(r);data.frame(file=basename(path),tests=nrow(d),failed=sum(d$failed),errors=sum(d$error),warnings=sum(d$warning),passed=sum(d$passed))})
out<-do.call(rbind,results);write.csv(out,'validation/astra-fixes-2026-10-27/core/focused-summary.csv',row.names=FALSE);print(out);stopifnot(sum(out$failed)==0,sum(out$errors)==0)
