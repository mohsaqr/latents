suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
G <- 60; m <- 8; n <- G*m
g <- rep(seq_len(G), each=m)
tr <- 1L+(runif(n)<c(.3,.7)[rep(rep(1:2,length.out=G),each=m)])
d <- data.frame(g=g, cov=rnorm(n), a=rnorm(n,c(-1.5,1.5)[tr]), b=rnorm(n,c(-1.3,1.3)[tr]), c=rnorm(n,c(1,-1)[tr]),
                y=rnorm(n,c(0,1)[tr]), yb=factor(rbinom(n,1,c(.3,.7)[tr]),labels=c("n","y")),
                u1=rbinom(n,1,c(.2,.8)[tr]), u2=rbinom(n,1,c(.8,.3)[tr]), u3=rbinom(n,1,c(.25,.75)[tr]), u4=rbinom(n,1,c(.3,.9)[tr]),
                cnt=rpois(n,c(1,5)[tr]), ord=sample(1:4,n,TRUE), w=rep(runif(G,.5,2),each=m), time=rep(1:m,G))
d$ordf <- pmin(pmax(round(c(1,3)[tr]+rnorm(n,0,.8)),1),4)
d$groupcov <- rep(rnorm(G),each=m)
try_it <- function(label, expr) {
  out <- withCallingHandlers(tryCatch({ v <- expr(); paste("OK", paste(format(head(v$standard_error %||% v$statistic,3),digits=4),collapse=" ")) },
    error=function(e) paste0("REFUSED [", paste(class(e)[1],collapse="/"),"]: ", substr(conditionMessage(e),1,90))),
    warning=function(w){ cat("   warn",class(w)[1],":",substr(conditionMessage(w),1,80),"\n"); invokeRestart("muffleWarning")})
  cat(sprintf("%-46s %s\n",label,out))
}
fits <- list(
 lpa = quote(lpa(d,c("a","b","c"),n_profiles=2,n_starts=2,seed=1)),
 lpa_VVV = quote(lpa(d,c("a","b","c"),n_profiles=2,model="VVV",n_starts=2,seed=1)),
 lca = quote(lca(d,c("u1","u2","u3","u4"),n_classes=2,n_starts=2,seed=1)),
 mlpa2 = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=2,n_starts=2,seed=1)),
 mlpa_direct = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,direct_effects="cov",n_starts=2,seed=1)),
 mlpa_count = quote(multilpa(d,c("a","b","cnt"),"g",n_profiles=2,n_group_classes=1,count="cnt",n_starts=2,seed=1)),
 mlpa_ord = quote(multilpa(d,c("a","b","ordf"),"g",n_profiles=2,n_group_classes=1,ordinal="ordf",n_starts=2,seed=1)),
 mlpa_prior = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,prior=prior_control(),n_starts=2,seed=1)),
 mlpa_w = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,weights="w",n_starts=2,seed=1)),
 mlpa_noise = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,noise=TRUE,n_starts=2,seed=1)),
 mlpa_cov = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,profile_covariates="cov",n_starts=2,seed=1)),
 add = quote(multilpa(d,c("a","b","c"),"g",n_group_classes=2,family="additive",n_starts=2,seed=1)),
 cross = quote(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=2,family="restricted_cross_level",n_starts=2,seed=1)),
 staged = quote(fit_staged(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=2,n_starts=2,seed=1)),
 lta = quote(lta(d,c("a","b","c"),id="g",time="time",n_profiles=2,n_starts=2,seed=1)),
 lta_occ = quote(lta(d,c("a","b","c"),id="g",time="time",n_profiles=2,transitions="occasion",n_starts=2,seed=1)),
 lta_stayer = quote(lta(d,c("a","b","c"),id="g",time="time",n_profiles=2,mover_stayer=TRUE,n_starts=2,seed=1))
)
built <- lapply(names(fits), function(nm) {
  tryCatch(withCallingHandlers(eval(fits[[nm]]), warning=function(w) invokeRestart("muffleWarning")),
    error=function(e){ cat("could not build",nm,":",substr(conditionMessage(e),1,120),"\n"); NULL })
})
names(built) <- names(fits)
print(lapply(built, class))
invisible(lapply(names(built), function(nm) {
  f <- built[[nm]]; if (is.null(f)) return(NULL)
  cat("\n==",nm,"\n")
  try_it("three_step bch", function() three_step(f,d,"y"))
  try_it("three_step ml", function() three_step(f,d,"y",method="ml"))
  try_it("three_step bch estimated", function() three_step(f,d,"y",step_one="estimated"))
  try_it("three_step modal estimated", function() three_step(f,d,"y",method="modal",step_one="estimated"))
  try_it("three_step categorical ml estimated", function() three_step(f,d,"yb",method="ml",step_one="estimated"))
  try_it("r3step covariates", function() r3step(f,d,"cov"))
  try_it("r3step by_group_class", function() r3step(f,d,"cov",by_group_class=TRUE))
  try_it("r3step initial+transition", function() r3step(f,d,initial_covariates="groupcov",transition_covariates="cov"))
}))
