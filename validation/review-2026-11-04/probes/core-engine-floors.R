suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(6)
caught <- function(e) { w <- character(); v <- withCallingHandlers(e, warning=function(c){w<<-c(w,class(c)[1]);invokeRestart("muffleWarning")}); attr(v,"w") <- unique(w); v }
n <- 300; z <- rep(1:2, c(150,150))
# (1) Poisson: one class has only zero counts. Mean at its floor
d1 <- data.frame(y=rnorm(n,c(0,4)[z]), ct=ifelse(z==1,0L,rpois(n,3)))
f1 <- caught(lpa(d1,c("y","ct"),2,count="ct",n_starts=3,seed=1))
cat("Poisson class of zeros: means",round(f1$count_means[,1],12)," boundary",f1$boundary," converged",f1$converged," warnings:",attr(f1,"w"),"\n")
r <- tryCatch(parameter_inference(f1),error=function(e) e); cat(" Wald:",if(inherits(r,"error")) paste("REFUSED",class(r)[1]) else "given","\n")
if(!inherits(r,"error")) print(subset(as.data.frame(r),grepl("count",parameter)))
print(get_results(f1,"count_means"))
# (2) Ordinal: one class always in the top category. Location diverges
d2 <- data.frame(y=rnorm(n,c(0,4)[z]), o=ifelse(z==1,4L,sample(1:4,n,TRUE)))
f2 <- caught(lpa(d2,c("y","o"),2,ordinal="o",n_starts=3,seed=1))
cat("Ordinal class in top category: locations",f2$ordinal_locations," intercepts",unlist(f2$ordinal_intercepts)," boundary",f2$boundary," converged",f2$converged," warnings:",attr(f2,"w"),"\n")
r <- tryCatch(parameter_inference(f2),error=function(e) e); cat(" Wald:",if(inherits(r,"error")) paste("REFUSED",class(r)[1],conditionMessage(r)) else "given","\n")
if(!inherits(r,"error")) print(subset(as.data.frame(r),grepl("ordinal",parameter)))
print(subset(get_results(f2,"ordinal"),select=c(profile,category,probability,location,location_standard_error)))
print(f2$iterations); print(tail(diff(f2$log_likelihood_history),3))
# (2b) Same, ordinal 2 cats, class-separated exactly, profile order swapped: reference class in the extreme
d3 <- data.frame(y=rnorm(n,c(0,4)[z]), o=ifelse(z==2,4L,sample(1:4,n,TRUE)))
f3 <- caught(lpa(d3,c("y","o"),2,ordinal="o",n_starts=3,seed=1))
cat("swap: locations",f3$ordinal_locations," intercepts",unlist(f3$ordinal_intercepts)," warnings:",attr(f3,"w"),"\n")
# LL hand-check of both
pr <- exp(.latents_ordinal_log_probabilities(f3$ordinal_intercepts$o, f3$ordinal_locations[,"o"])); print(round(pr,6))
# (3) Numeric ordinal with a gap: the scores skip nothing
d4 <- data.frame(y=rnorm(n,c(0,3)[z]), o=sample(c(1L,2L,5L),n,TRUE))
f4 <- caught(lpa(d4,c("y","o"),2,ordinal="o",n_starts=2,seed=1))
d5 <- d4; d5$o <- factor(d5$o,levels=1:5,ordered=TRUE)
f5 <- caught(lpa(d5,c("y","o"),2,ordinal="o",n_starts=2,seed=1))
cat("numeric 1,2,5: n_par",f4$n_parameters," levels",f4$extra_data$ordinal_levels$o," LL",f4$log_likelihood,"\n")
cat("factor levels 1:5: n_par",f5$n_parameters," levels",f5$extra_data$ordinal_levels$o," LL",f5$log_likelihood," warnings",attr(f5,"w"),"\n")
# (4) Row permutation with interleaved groups: evaluate-only LL, hand group posterior
J<-25; gid<-rep(seq_len(J),each=6); nn<-length(gid); zz<-rbinom(nn,1,ifelse(rep(rbinom(J,1,.5),each=6)==1,.8,.2))
dg<-data.frame(g=gid,y=rnorm(nn,zz*2),u=rnorm(nn,-zz))
fg<-caught(multilpa(dg,c("y","u"),"g",2,2,n_starts=3,seed=1,tol=1e-10))
perm<-sample(nn); dp<-dg[perm,]
fp<-caught(multilpa(dp,c("y","u"),"g",2,2,start=starting_values(fg),max_iter=0L))
cat("permuted LL diff",fp$log_likelihood-fg$log_likelihood,"\n")
cat("row posteriors aligned",max(abs(unname(fp$subject_posteriors)-unname(fg$subject_posteriors)[perm,])),"\n")
ord<-match(rownames(fg$group_posteriors),rownames(fp$group_posteriors))
cat("group posteriors aligned by id",max(abs(fp$group_posteriors[ord,]-fg$group_posteriors)),"\n")
hg<-t(vapply(unique(dp$g),function(gg){rows<-dp[dp$g==gg,];l<-vapply(1:2,function(h){log(fg$group_probabilities[h])+sum(apply(rows[c("y","u")],1,function(r){log(sum(fg$profile_probabilities[h,]*vapply(1:2,\(k) prod(dnorm(r,fg$means[k,],sqrt(fg$variances[k,]))),0)))}))},0);exp(l-max(l))/sum(exp(l-max(l)))},numeric(2)))
cat("hand group posterior vs fit",max(abs(hg-unname(fp$group_posteriors))),"\n")
# (5) prior + several starts: selected by LL, whereas the iteration climbs the posterior
logpost <- function(f,x) { # VVI mclust-form prior: Normal-InverseGamma per indicator with kappa, nu, s^2
  k <- f$n_profiles; d <- ncol(x); h <- f$prior_parameters
  lp <- sum(vapply(seq_len(k),function(p) sum(dnorm(f$means[p,],h$mean,sqrt(f$variances[p,]/h$shrinkage),log=TRUE) +
      (-(h$dof/2+1)*log(f$variances[p,])-h$scale/(2*f$variances[p,]))),0))
  f$log_likelihood+lp }
xx <- rbind(matrix(rnorm(160*3),ncol=3),matrix(rnorm(20*3,3,.3),ncol=3),matrix(rnorm(6*3,-4,.05),ncol=3)); colnames(xx)<-c("a","b","c"); dx<-as.data.frame(xx)
sols <- do.call(rbind,lapply(1:25,function(s){f<-caught(lpa(dx,c("a","b","c"),4,model="VVI",prior=prior_control(),n_starts=1,seed=s));data.frame(seed=s,ll=f$log_likelihood,lp=logpost(f,xx),conv=f$converged)}))
sols$ll<-round(sols$ll,4); sols$lp<-round(sols$lp,4); u<-unique(sols[c("ll","lp")]); print(u[order(-u$ll),][1:min(8,nrow(u)),])
cat("best LL row vs best posterior row same?",which.max(sols$ll)==which.max(sols$lp)," LLtop",max(sols$ll)," its lp",sols$lp[which.max(sols$ll)]," lp top",max(sols$lp),"\n")
