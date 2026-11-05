suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
library(testthat)
ns<-asNamespace("latents"); root<-"/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"
patch<-function(name,from,to){ f<-get(name,ns); src<-paste(deparse(f,width.cutoff=500L),collapse="\n"); stopifnot(grepl(from,src,fixed=TRUE)); g<-eval(parse(text=sub(from,to,src,fixed=TRUE)),envir=ns); environment(g)<-ns
  unlockBinding(name,ns); assign(name,g,ns); f}
restore<-function(name,f){assign(name,f,ns)}
run<-function(files){ r<-lapply(files,function(fl){ tt<-as.data.frame(test_file(file.path(root,fl),package="latents",load_package="none",reporter="silent",stop_on_failure=FALSE)); data.frame(file=fl,tests=sum(tt$nb),failed=sum(tt$failed),errors=sum(tt$error))}); do.call(rbind,r)}
sab<-function(label,name,from,to,files){ f<-patch(name,from,to); t0<-Sys.time(); r<-tryCatch(run(files),error=function(e)data.frame(file=NA,tests=NA,failed=NA,errors=conditionMessage(e))); restore(name,f); cat("\n==",label,sprintf("(%.0fs)",as.numeric(difftime(Sys.time(),t0,units="secs"))),"\n"); print(r,row.names=FALSE)}
pool<-"test-pool-imputations.R"
print(run(pool))
sab("Rubin total without (1+1/m)",".multilpa_rubin","total <- within + (1 + 1/m) * between","total <- within + between",pool)
sab("Rubin df (1+1/r) not squared",".multilpa_rubin","(1 + 1/riv)^2","(1 + 1/riv)",pool)
sab("pooled interval uses qnorm",".multilpa_rubin","stats::qt((1 + level)/2, df)","stats::qnorm((1 + level)/2)",pool)
sab("fmi wrong",".multilpa_rubin","(riv + 2/(df + 3))/(riv + 1)","riv/(riv + 1)",pool)
sab("CR1 factor gone",".inference_cluster_factor","n_units/(n_units - 1)","1",c("test-cluster-robust-cr1.R"))
sab("covariate variance delta wrong",".multilpa_cov_front_jacobian","exp(theta[is_variance])","1",c("test-covariate-inference.R"))
sab("response simplex Jacobian: drops first instead of last",".inference_simplex_jacobian","block[, -1L, drop = FALSE] else","block[, -1L, drop = FALSE] else block[, -1L, drop = FALSE]; if (FALSE) NULL else",c("test-covariate-measurement.R"))
