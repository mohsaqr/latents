suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
root <- "/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"
ns <- asNamespace("latents")
run_files <- function(files) {
  vapply(files, function(f) {
    res <- as.data.frame(testthat::test_file(file.path(root,f), package="latents", load_package="none", reporter="silent", stop_on_failure=FALSE))
    sum(res$failed) + sum(res$error)
  }, 0)
}
mutate <- function(name, fun, pattern, replacement, files) {
  orig <- get(fun, envir=ns)
  src <- paste(deparse(orig, width.cutoff=500L, control=c("keepNA","keepInteger","niceNames","showAttributes")), collapse="\n")
  stopifnot(grepl(pattern, src, fixed=TRUE))
  new <- eval(parse(text=sub(pattern, replacement, src, fixed=TRUE)), envir=ns)
  environment(new) <- environment(orig)
  unlockBinding(fun, ns); assign(fun, new, envir=ns)
  on.exit({ assign(fun, orig, envir=ns); lockBinding(fun, ns) }, add=TRUE)
  t0 <- proc.time()[[3]]
  out <- run_files(files)
  cat(sprintf("%-34s %-5s failures: %s  (%.0fs)\n", name, if (sum(out)>0) "CAUGHT" else "MISSED",
              paste(sprintf("%s=%d", sub("test-","",sub(".R$","",names(out))), out), collapse=" "), proc.time()[[3]]-t0))
}
