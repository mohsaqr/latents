source("/Users/mohammedsaqr/Documents/Github/latents/validation/review-2026-11-04/probes/three-step-mutants-common.R")
setwd("/Users/mohammedsaqr/Documents/Github/latents")
assign("root", "/Users/mohammedsaqr/Documents/Github/latents/equivalence", envir=globalenv())
f <- "test-latentgold-step3-markov.R"
print(run_files(f))
mutate("lta: D transposed in emissions", ".lta_three_step_inputs",
  "t(errors)[modal, , drop = FALSE]", "errors[modal, , drop = FALSE]", f)
mutate("lta: chain posteriors, not occasion", ".lta_occasion_posteriors",
  "joint/rowSums(joint)", "posteriors", f)
mutate("lta: robust meat doubled", ".lta_three_step_fit",
  ".inference_sandwich(inverse, crossprod(scores), nrow(scores))", ".inference_sandwich(inverse, 2 * crossprod(scores), nrow(scores))", f)
mutate("lta: adjust includes intercepts", ".lta_three_step",
  "tested <- fitted$labels$term != \"(Intercept)\"", "tested <- rep(TRUE, length(p_value))", f)
# rescaling: is the 0.1% gap of step-one SE tol-related?
set.seed(1); n <- 500; tr <- sample.int(3,n,TRUE,c(.4,.35,.25)); mus <- rbind(c(-1.5,-1),c(1.5,-1),c(0,1.5))
d <- data.frame(a=rnorm(n,mus[tr,1]),b=rnorm(n,mus[tr,2]),y=rnorm(n,c(0,1,3)[tr],c(1,1.3,.8)[tr]))
d2 <- transform(d,a=a*1e-4,b=b*1e4)
invisible(lapply(c(1e-10,1e-13),function(tol){
  f1 <- lpa(d,c("a","b"),n_profiles=3,n_starts=5,seed=1,tol=tol,max_iter=50000); f2 <- lpa(d2,c("a","b"),n_profiles=3,n_starts=5,seed=1,tol=tol,max_iter=50000)
  cat("tol",tol," LL gap",f1$log_likelihood-f2$log_likelihood-(-n*log(1)),"\n")
  a <- three_step(f1,d,"y"); b <- three_step(f2,d2,"y"); ae <- three_step(f1,d,"y",step_one="estimated"); be <- three_step(f2,d2,"y",step_one="estimated")
  cat(" fixed",format(sort(a$standard_error),digits=7),"|",format(sort(b$standard_error),digits=7),"\n estimated",format(sort(ae$standard_error),digits=7),"|",format(sort(be$standard_error),digits=7),"\n")
}))
