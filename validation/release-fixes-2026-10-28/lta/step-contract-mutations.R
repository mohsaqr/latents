Sys.setenv(NOT_CRAN='true')
pkgload::load_all('.',quiet=TRUE)
mode <- commandArgs(trailingOnly=TRUE)[1]
ns <- asNamespace('latents')
if (mode=='ignored_step') {
 name <- 'parameter_inference.multilpa_lta'
 mutant <- get(name,ns)
 body(mutant) <- parse(text=sub('.lta_inference(x, vcov_type, step)', '.lta_inference(x, vcov_type)', paste(deparse(body(mutant)),collapse='\n'),fixed=TRUE))[[1]]
 unlockBinding(name,ns);assign(name,mutant,ns);lockBinding(name,ns)
 if (mode=='ignored_step') registerS3method('parameter_inference','multilpa_lta',mutant,envir=ns)
}
if (mode=='recomputed_table') {
 name <- '.lta_coefficient_table'
 mutant <- get(name,ns)
 body(mutant) <- bquote({inference <- NULL; .(body(mutant))})
 unlockBinding(name,ns);assign(name,mutant,ns);lockBinding(name,ns)
 if (mode=='ignored_step') registerS3method('parameter_inference','multilpa_lta',mutant,envir=ns)
}
expressions <- parse('tests/testthat/test-engine-fixes.R')
selected <- Filter(function(x) {
 text <- paste(deparse(x),collapse='\n')
 startsWith(text,'fixes_transition_fit <-') || startsWith(text,'test_that("transition parameter_inference() uses')
},as.list(expressions))
stopifnot(length(selected)==2L)
path <- file.path('validation/release-fixes-2026-10-28/lta',paste0('test-step-',mode,'.R'))
writeLines(unlist(lapply(selected,deparse)),path)
results <- testthat::test_file(path,reporter='summary')
print(as.data.frame(results)[c('test','failed','error','passed')])
failed <- sum(as.data.frame(results)$failed)+sum(as.data.frame(results)$error)
stopifnot(if(mode=='baseline') failed==0 else failed>0)
