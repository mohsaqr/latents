suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(5)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {invokeRestart("muffleWarning")})
J <- 50; nj <- sample(3:7, J, TRUE); gid <- rep(seq_len(J), nj); n <- length(gid)
gc <- rbinom(J,1,.5)[gid]; z <- rbinom(n,1,ifelse(gc==1,.75,.25))+1L
dat <- data.frame(g=gid, y1=rnorm(n,c(0,1.2)[z],1), y2=rnorm(n,c(1,0)[z],c(.8,1.2)[z]), y3=rnorm(n,c(0,.8)[z],1),
  cat=vapply(z,\(k) sample(1:3,1,prob=if(k==1) c(.5,.3,.2) else c(.2,.3,.5)),1L),
  ord=vapply(z,\(k) sample(1:4,1,prob=if(k==1) c(.4,.3,.2,.1) else c(.1,.2,.3,.4)),1L),
  cnt=rnbinom(n,size=3,mu=c(2,4)[z]), cv=rnorm(n), w=c(1,2,3)[(gid%%3)+1])
miss <- dat; for (v in c("y1","y2","y3","cat","ord","cnt")) miss[[v]][sample(n, round(.08*n))] <- NA
base <- list(data=miss, vars=c("y1","y2","y3"), id="g", n_profiles=3, n_group_classes=2, missing="fiml", tol=1e-10, max_iter=3000)
cases <- list(
  list(model="VVI"), list(model="EEI"), list(model="EVI"), list(model="VEI"), list(model="EEE"),list(model="VEE"),list(model="EVE"), list(model="VVE"), list(model="EEV"),list(model="VEV"), list(model="EVV"), list(model="VVV"),
  list(vars=c("y1","cat","ord","cnt"),categorical="cat",ordinal="ord",count="cnt",count_model="negative_binomial"),
  list(weights="w"), list(direct_effects="y1 ~ cv"), list(local_dependence="y1:y2"),
  list(vars=c("cat","ord"),categorical=c("cat","ord")))
tab <- do.call(rbind, lapply(cases, function(cc) {
  args <- modifyList(base, cc); args$n_starts <- 1L; args$seed <- 11L
  runs <- lapply(c("squarem","none"), function(a) tryCatch(quiet(do.call(multilpa, c(args, list(acceleration=a)))), error=function(e) e))
  s <- runs[[1]]; p <- runs[[2]]
  if (inherits(s,"error")||inherits(p,"error")) return(data.frame(case=paste(cc,collapse="/"),err=paste(conditionMessage(if(inherits(s,"error")) s else p)),LLs=NA,LLp=NA,minstep=NA,lastdiff=NA,convS=NA,convP=NA,itS=NA,itP=NA,lenH=NA))
  h <- s$log_likelihood_history
  data.frame(case=paste(unlist(cc),collapse="/"),err="",LLs=s$log_likelihood,LLp=p$log_likelihood,minstep=min(diff(h)),
    lastdiff=h[length(h)]-s$log_likelihood,convS=s$converged,convP=p$converged,itS=s$iterations,itP=p$iterations,lenH=length(h))
}))
print(tab,digits=9)
# Small max_iter: does the count respect the budget, and is `converged` honest?
for (mi in c(1L,2L,3L,4L,5L,7L)) { f<-quiet(multilpa(miss,c("y1","y2","y3"),"g",3,2,missing="fiml",n_starts=1,seed=11,max_iter=mi))
  cat("max_iter",mi,"iterations",f$iterations,"converged",f$converged,"hist",length(f$log_likelihood_history),"\n") }
# Evaluate the final point: LL recomputed by max_iter=0 equals the reported one, and one more EM step gains < threshold
f<-quiet(multilpa(miss,c("y1","y2","y3"),"g",3,2,missing="fiml",n_starts=1,seed=11))
g<-quiet(multilpa(miss,c("y1","y2","y3"),"g",3,2,missing="fiml",n_starts=1,seed=11,start=starting_values(f),max_iter=1L,acceleration="none"))
cat("default tol: one more EM step gains", g$log_likelihood-f$log_likelihood, " threshold", 1e-8*sum(!is.na(miss[c("y1","y2","y3")])),"\n")
# gap to the maximum at default tol
h<-quiet(multilpa(miss,c("y1","y2","y3"),"g",3,2,missing="fiml",n_starts=1,seed=11,start=starting_values(f),tol=1e-13,max_iter=50000))
cat("default tol stopped short of the maximum by", h$log_likelihood-f$log_likelihood, "iters",h$iterations,"\n")
# Reported optimizer flags for a held fit under SQUAREM
st<-starting_values(quiet(multilpa(dat,c("y1","y2","cat"),"g",3,1,categorical="cat",n_starts=2,seed=2)))
fx<-quiet(multilpa(dat,c("y1","y2","cat"),"g",3,2,categorical="cat",start=st,fixed="measurement",n_starts=1,seed=2))
cat("held exact after SQUAREM: means", identical(unname(fx$means),unname(st$means)), " var", identical(unname(fx$variances),unname(st$variances)), " resp", all.equal(unname(fx$response_probabilities[[1]]),unname(st$response_probabilities[[1]])), "\n")
print(max(abs(unname(fx$means)-unname(st$means))))
print(max(abs(unname(fx$variances)-unname(st$variances))))
