suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(6)
caught <- function(e) { w <- character(); v <- withCallingHandlers(e, warning=function(c){w<<-c(w,class(c)[1]);invokeRestart("muffleWarning")}); attr(v,"w") <- unique(w); v }
n <- 300; z <- rep(1:2, c(150,150))
d3 <- data.frame(y=rnorm(n,c(0,6)[z]), o=ifelse(z==2,4L,sample(1:4,n,TRUE)))
f3 <- caught(lpa(d3,c("y","o"),2,ordinal="o",n_starts=3,seed=1,tol=1e-10))
cat("separated ordinal: locations",f3$ordinal_locations," intercepts",unlist(f3$ordinal_intercepts)," conv",f3$converged," boundary",f3$boundary," warnings:",attr(f3,"w"),"\n")
r <- tryCatch(caught(parameter_inference(f3)),error=function(e) e)
cat(" Wald:",if(inherits(r,"error")) paste("REFUSED",class(r)[1],conditionMessage(r)) else "given","\n")
if(!inherits(r,"error")) print(subset(as.data.frame(r),grepl("ordinal",parameter))[c("term","parameter","estimate","standard_error","p_value","conf_low","conf_high")])
print(caught(get_results(f3,"ordinal"))[c("profile","category","probability","location","location_standard_error")])
# Same data, tol and seed with the profile order swapped via the data: is the reported location a stable number?
sapply(c(1e-6,1e-8,1e-10,1e-12),function(t){f<-caught(lpa(d3,c("y","o"),2,ordinal="o",n_starts=1,seed=1,tol=t));c(tol=t,loc=f$ordinal_locations[1,1],iter=f$iterations,ll=f$log_likelihood)})|>print()
# Sanity check on the M-step: Newton with zero counts for categories
cc <- rbind(c(40,50,30,40),c(0,0,0,150))
print(.latents_ordinal_maximize(cc,c(0,0,0),c(0,0)))
# Prior: does best-by-likelihood choose a lower posterior?
logpost <- function(f) { h <- f$prior_parameters; k<-f$n_profiles
  f$log_likelihood+sum(vapply(seq_len(k),function(p) sum(dnorm(f$means[p,],h$mean,sqrt(f$variances[p,]/h$shrinkage),log=TRUE) + (-(h$dof/2+1)*log(f$variances[p,])-h$scale/(2*f$variances[p,]))),0)) }
set.seed(9); xx <- matrix(rt(120*3,4),ncol=3); colnames(xx)<-c("a","b","c"); dx<-as.data.frame(xx)
sols <- do.call(rbind,lapply(1:30,function(s){f<-caught(lpa(dx,c("a","b","c"),4,model="VVI",prior=prior_control(),n_starts=1,seed=s,search="kmeans"));data.frame(seed=s,ll=round(f$log_likelihood,4),lp=round(logpost(f),4),conv=f$converged)}))
u<-unique(sols[c("ll","lp")]); print(u[order(-u$ll),][seq_len(min(8,nrow(u))),])
cat("LL-best has highest log posterior?",which.max(sols$ll)==which.max(sols$lp)," top LL",max(sols$ll)," its lp",sols$lp[which.max(sols$ll)],"; top lp",max(sols$lp)," its ll",sols$ll[which.max(sols$lp)],"\n")
