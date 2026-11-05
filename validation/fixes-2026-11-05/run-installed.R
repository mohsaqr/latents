args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
options(mc.cores=2L)
cat(R.version.string,"\nVersion",as.character(utils::packageVersion("latents")),"\nLibrary",find.package("latents"),"\n")
run <- function(filter, name) {
  start<-proc.time()
  results<-testthat::test_dir(file.path(args[[1]],"tests/testthat"),filter=filter,
    package="latents",load_package="installed",reporter="summary",stop_on_failure=FALSE)
  saveRDS(results,file.path(args[[2]],paste0(name,"-results.rds")))
  table<-as.data.frame(results)
  write.csv(table[,setdiff(names(table),"result"),drop=FALSE],file.path(args[[2]],paste0(name,"-results.csv")),row.names=FALSE)
  cat("\n",name,"TOTALS\n"); print(colSums(table[,c("passed","failed","error","warning","skipped")]))
  print(proc.time()-start)
  invisible(table)
}
preflight<-run("^fixes-1105-","installed-preflight")
stopifnot(!any(preflight$failed),!any(preflight$error))
run(NULL,"full-suite")
