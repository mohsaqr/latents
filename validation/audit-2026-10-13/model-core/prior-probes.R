pkgload::load_all('.',quiet=TRUE)
set.seed(20)
x<-matrix(rnorm(200),100,2)
post<-matrix(1,100,1)
run<-function(code,kappa=100,scale=1,dof=4) {
 p<-prior_control(shrinkage=kappa,mean=c(20,20),scale=if(code=='VEV')diag(scale,2)else scale,dof=dof)
 h<-latents:::.multilpa_prior_parameters(p,x,1,code)
 latents:::.multilpa_prior_maximize(x,post,code,h,1e-12)
}
results<-list(eei_means=run('EEI')$means,vvi_means=run('VVI')$means,
 eii_scale_1=run('EII',scale=1)$variances,eii_scale_1e6=run('EII',scale=1e6)$variances,
 vev_dof_4=run('VEV',dof=4)$variances,vev_dof_400=run('VEV',dof=400)$variances)
print(results)
stopifnot(identical(results$eii_scale_1,results$eii_scale_1e6),identical(results$vev_dof_4,results$vev_dof_400))
saveRDS(results,'validation/audit-2026-10-13/model-core/prior-probes.rds')
