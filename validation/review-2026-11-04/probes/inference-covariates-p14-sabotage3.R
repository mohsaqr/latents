suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
library(testthat)
ns<-asNamespace("latents"); root<-"/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"; eq<-"/Users/mohammedsaqr/Documents/Github/latents/tests/equivalence"
patch<-function(name,from,to){ f<-get(name,ns); src<-paste(deparse(f,width.cutoff=500L),collapse="\n"); stopifnot(grepl(from,src,fixed=TRUE)); g<-eval(parse(text=sub(from,to,src,fixed=TRUE)),envir=ns); environment(g)<-ns; unlockBinding(name,ns); assign(name,g,ns); f}
run<-function(files){ do.call(rbind,lapply(files,function(fl){ tt<-as.data.frame(test_file(fl,package="latents",load_package="none",reporter="silent",stop_on_failure=FALSE)); data.frame(file=basename(fl),tests=sum(tt$nb),failed=sum(tt$failed),errors=sum(tt$error),skipped=sum(tt$skipped))}))}
sab<-function(label,name,from,to,files){ f<-patch(name,from,to); r<-tryCatch(run(files),error=function(e)data.frame(file=NA,tests=NA,failed=NA,errors=conditionMessage(e),skipped=NA)); assign(name,f,ns); cat("\n==",label,"\n"); print(r,row.names=FALSE)}
cat("baseline\n"); print(run(file.path(root,"test-audit-1101-accessors.R")),row.names=FALSE)
sab("Rubin total without (1+1/m), equivalence file",".multilpa_rubin","total <- within + (1 + 1/m) * between","total <- within + between",file.path(eq,"test-pool-imputations.R"))
sab("simplex reference wrong, other covariate files",".multilpa_cov_natural_jacobian",'.inference_simplex_jacobian(block[profile, ], "last")','.inference_simplex_jacobian(block[profile, ], "first")',file.path(root,c("test-categorical-inference.R","test-covariate-structure-inference.R","test-audit-1101-covariates-diagnostics.R")))
