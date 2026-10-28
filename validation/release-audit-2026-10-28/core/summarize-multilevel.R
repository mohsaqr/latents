x <- readRDS('validation/release-audit-2026-10-28/core/multilevel-units.rds')
out <- do.call(rbind,lapply(x,function(q) {
 n <- if(grepl('cross_level',q$family)) 1080 else 960
 data.frame(family=q$family,likelihood_gap=q$small$log_likelihood-q$base$log_likelihood+n*log(1e-6))
}))
print(out)
stopifnot(max(abs(out$likelihood_gap))<1e-6)
