Sys.setenv(OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1)
pkgload::load_all('.', quiet=TRUE)
a <- readRDS('validation/release-audit-2026-10-28/inference/units-probe.rds')
f0 <- a[[1]]$fit
v0 <- a[[1]]$vcov
mean_i <- startsWith(rownames(v0),'mean.')
# Map the fitted object exactly, removing optimisation from the experiment.
rescale_fit <- function(f,s) {
 f$x <- f$x*s; f$centers <- f$centers*s
 f$measurement <- lapply(f$measurement,function(b){b$means<-b$means*s;b$variances<-b$variances*s^2;b})
 f$min_variance <- f$min_variance*s^2
 f
}
results <- lapply(c(1,1e-2,1e-3,3e-4,1e-4,3e-5,1e-5,1e5,1e6),function(s){
 f <- rescale_fit(f0,s)
 v <- tryCatch(vcov(f),error=function(e)e)
 inf <- tryCatch(parameter_inference(f),error=function(e)e)
 out <- data.frame(scale=s,vcov_status=if(inherits(v,'error')) conditionMessage(v) else 'OK',
    max_se_relative_error=if(inherits(v,'error')) NA_real_ else max(abs(sqrt(diag(v)/diag(v0))/ifelse(mean_i,s,1)-1)),
    inference_status=if(inherits(inf,'error')) conditionMessage(inf) else 'OK')
 print(out,digits=12,row.names=FALSE)
 out
})
write.csv(do.call(rbind,results),'validation/release-audit-2026-10-28/inference/units-inference.csv',row.names=FALSE)
# Inspect original vs standardized curvature at tiny units with a safe step.
f <- rescale_fit(f0,1e-5)
theta <- latents:::.lta_theta(f)
score <- function(x) colSums(latents:::.lta_group_scores(x,f))
info <- latents:::.inference_score_information(score,theta,1e-8,scheme='five_point',h_rule='max1')
D <- diag(1/sqrt(diag(info)))
scaled_info <- D %*% info %*% D
fine <- D %*% solve(scaled_info) %*% D
cat('raw information eigenvalues:',range(eigen(info,symmetric=TRUE)$values),'\n')
cat('scaled information eigenvalues:',range(eigen(scaled_info,symmetric=TRUE)$values),'\n')
cat('stable inverse relative SE error vs original units:',max(abs(sqrt(diag(fine)/diag(v0))/ifelse(mean_i,1e-5,1)-1)),'\n')
# coef() reports natural measurement locations, covariance is translation-invariant.
cat('maximum coef/profile mean mismatch:',max(abs(unname(coef(f0)[mean_i])-as.vector(t(f0$measurement[[1]]$means)))),'\n')
cat('coef and vcov names identical:',identical(names(coef(f0)),rownames(v0)),'\n')
f <- rescale_fit(f0,1e-4)
r0 <- as.data.frame(parameter_inference(f0))
r1 <- as.data.frame(parameter_inference(f))
cat('transition/initial standard errors at scale 1e-4:\n')
print(data.frame(block=r0$block,from=r0$from,to=r0$to,term=r0$term,
 original_se=r0$standard_error,scaled_se=r1$standard_error,
 relative_difference=r1$standard_error/r0$standard_error-1),row.names=FALSE)
