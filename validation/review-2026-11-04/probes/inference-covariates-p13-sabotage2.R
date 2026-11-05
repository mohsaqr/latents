suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
library(testthat)
ns<-asNamespace("latents"); root<-"/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"
patch<-function(name,from,to){ f<-get(name,ns); src<-paste(deparse(f,width.cutoff=500L),collapse="\n"); stopifnot(grepl(from,src,fixed=TRUE)); g<-eval(parse(text=sub(from,to,src,fixed=TRUE)),envir=ns); environment(g)<-ns; unlockBinding(name,ns); assign(name,g,ns); f}
run<-function(files){ do.call(rbind,lapply(files,function(fl){ tt<-as.data.frame(test_file(file.path(root,fl),package="latents",load_package="none",reporter="silent",stop_on_failure=FALSE)); data.frame(file=fl,tests=sum(tt$nb),failed=sum(tt$failed),errors=sum(tt$error),skipped=sum(tt$skipped))}))}
sab<-function(label,name,from,to,files){ f<-patch(name,from,to); r<-tryCatch(run(files),error=function(e)data.frame(file=NA,tests=NA,failed=NA,errors=conditionMessage(e),skipped=NA)); assign(name,f,ns); cat("\n==",label,"\n"); print(r,row.names=FALSE)}

files<-grep("pool|rubin",list.files(root),value=TRUE); print(files)
hits<-list.files(root,pattern="^test.*\\.R$",full.names=TRUE); hits<-hits[vapply(hits,function(p) any(grepl("multilpa_rubin|pool_imputations",readLines(p,warn=FALSE))),TRUE)]; print(basename(hits))
sab("Rubin total without (1+1/m), every testthat file naming pool",".multilpa_rubin","total <- within + (1 + 1/m) * between","total <- within + between",basename(hits))
sab("Rubin df not squared",".multilpa_rubin","(1 + 1/riv)^2","(1 + 1/riv)",basename(hits))
sab("response simplex Jacobian reference wrong",".multilpa_cov_natural_jacobian",'.inference_simplex_jacobian(block[profile, ], "last")','.inference_simplex_jacobian(block[profile, ], "first")',c("test-covariate-measurement.R","test-covariate-inference.R"))
