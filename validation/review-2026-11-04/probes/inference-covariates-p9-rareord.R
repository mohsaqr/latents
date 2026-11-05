suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
G<-50;n<-6
d<-do.call(rbind,lapply(seq_len(G),function(g){k<-rbinom(n,1,.5)+1L; data.frame(g=g,y=rnorm(n,c(-1.5,1.5)[k]),o=sapply(k,function(j) sample(1:3,1,prob=list(c(.6,.3,.1),c(.1,.3,.6))[[j]])))}))
d$o[d$g %in% c(3,17) & d$y>0][1:2]<-4L   # rare top level in two groups only
print(table(d$o)); print(tapply(d$o==4,d$g,any)[c(3,17)])
for(kind in c("numeric","ordered")){
  dd<-d; if(kind=="ordered") dd$o<-factor(dd$o,levels=1:4,ordered=TRUE)
  fit<-suppressWarnings(multilpa(dd,c("y","o"),"g",n_profiles=2,n_group_classes=1,ordinal="o",n_starts=3,seed=1))
  cat("\n---",kind,": n_par",fit$n_parameters," n coef",length(coef(fit)),"\n")
  warn<-character()
  b<-tryCatch(withCallingHandlers(parameter_inference(fit,method="bootstrap",iter=40,n_starts=2,seed=3),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)conditionMessage(e))
  if(is.character(b)) {cat("error:",b,"\n")} else {
   print(as.data.frame(b)[,c("outcome","term","parameter","estimate","standard_error","conf_low","conf_high")],digits=4)
   cat("n_valid",attr(b,"n_valid"),"of 40\n"); print(table(attr(b,"messages"),useNA="ifany")) }
  print(unique(warn))
  w<-parameter_inference(fit); print(w[,c("outcome","term","parameter","estimate","standard_error")],digits=4)
}
## Direct effect with covariates
d$z<-rnorm(nrow(d)); r<-tryCatch(multilpa(d,c("y"),"g",n_profiles=2,n_group_classes=1,profile_covariates="z",direct_effects="y~z",n_starts=2,seed=1),error=function(e)conditionMessage(e)); cat("direct+cov:",if(is.character(r)) r else "accepted","\n")
## additive bootstrap claim
