pkgload::load_all('.',quiet=TRUE)
s<-readRDS('/tmp/astra-lta-1e+08.rds')
print(c(base=s$a$log_likelihood,padded=s$b$log_likelihood))
print(c(base_boundary=s$a$boundary,padded_boundary=s$b$boundary))
original_chart<-latents:::.lta_design_chart
assignInNamespace('.lta_design_chart', function(designs,layout) {
 keep<-seq_len(35)
 designs$initial<-designs$initial[keep,,drop=FALSE]
 designs$transition<-lapply(designs$transition,function(d)d[keep,,drop=FALSE])
 layout$within<-layout$within[keep,,drop=FALSE]
 original_chart(designs,layout)
},ns='latents')
f<-suppressWarnings(lta(s$padded,c('y1','y2'),'id',2,'time',weights='sampling',transition_covariates='w',initial_covariates='z',n_starts=1,seed=1,max_iter=500))
print(c(correct_chart_ll=f$log_likelihood,gap=f$log_likelihood-s$a$log_likelihood,converged=f$converged))
saveRDS(f,'/tmp/astra-lta-chart-corrected.rds')
