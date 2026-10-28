pkgload::load_all('.',quiet=TRUE)
library(testthat)
test_that('EVI retains equal volumes when its floor is active', {
 v<-latents:::.multilpa_structure_variances(rbind(c(.001,10),c(1,1)),c(10,10),'EVI',.1)
 expect_equal(prod(v[1,]),prod(v[2,]),tolerance=1e-10)
})
test_that('VEI retains common shape when its floor is active', {
 v<-latents:::.multilpa_structure_variances(rbind(c(.01,.1),c(1,10)),c(10,10),'VEI',.1)
 expect_equal(v[1,1]/v[2,1],v[1,2]/v[2,2],tolerance=1e-10)
})
test_that('growth random-intercept quadratic retains within-person scatter', {
 y<-c(1e8-1,1e8,1e8+1)
 spec<-list(x=matrix(1,3,1),z=matrix(numeric(),3,0),random_design=matrix(1,3,1),y=y,
            offset=rep(0,3),n=3,group_index=rep(1L,3))
 s<-latents:::.growth_statistics(spec)
 actual<-latents:::.growth_class_terms(s,0,matrix(1e8,1,1),1)$log_density
 q<-sum((y-mean(y))^2)+3*mean(y)^2/(1+3e16)
 expected<--.5*(3*log(2*pi)+log(1+3e16)+q)
 expect_equal(actual,expected,tolerance=1e-8)
})
test_that('requested person centering is applied to a covariate fit', {
 p<-readRDS(file.path(pkgload::pkg_path(),'validation/audit-2026-10-13/model-core/dispatch-probes.rds'))
 expect_equal(unname(p$person$means),unname(p$manual$means),tolerance=1e-8)
})
