files <- c(Sys.glob('R/kernel*.R'),Sys.glob('R/gaussian*.R'),Sys.glob('R/covariance*.R'), paste0('R/', c('categorical','extra-indicators','multilpa','lpa','multilca','covariates','covariate-inference','lta-engine','lta-inference','lta-bootstrap','transition-inference','growth-mixture','growth-mixture-engine','growth-mixture-inference','mixture-regression','mixture-regression-engine','mixture-regression-inference','mixture-regression-ordinal','mixture-regression-starts','additive','additive-inference','cross-level','weights','noise','prior','fixed-measurement'),'.R'))
test_files <- list.files('tests',pattern='\\.R$',recursive=TRUE,full.names=TRUE)
test_text <- vapply(test_files,function(path)paste(readLines(path,warn=FALSE),collapse='\n'),character(1))
pkgload::load_all('.', quiet=TRUE)
ns <- asNamespace('latents')
all <- lapply(files, function(file) {
  parsed <- parse(file,keep.source=TRUE)
  functions <- vapply(parsed,function(e) is.call(e)&& identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')),logical(1))
  expressions <- parsed[functions]
  source_refs <- attr(parsed,'srcref')[functions]
  records <- Map(function(expression,source_ref) {
    name <- as.character(expression[[2]])
    fn <- eval(expression[[3]], envir=ns)
    text <- capture.output(codetools::checkUsage(fn,name=name))
    refs <- test_files[vapply(test_text,grepl,logical(1),pattern=name,fixed=TRUE)]
    data.frame(file=file,function_name=name,line_start=source_ref[1],line_end=source_ref[3],
               parameters=paste(names(formals(fn)),collapse=';'), test_files=paste(refs,collapse=';'),
               direct_test_reference=length(refs)>0,usage_diagnostics=paste(text,collapse=' | '),
               review='static review; targeted runtime probes and full-suite evidence recorded separately')
  },expressions,source_refs)
  do.call(rbind,records)
})
inventory <- do.call(rbind,all)
write.csv(inventory,'validation/audit-2026-10-13/model-core/function-inventory.csv',row.names=FALSE)
cat(nrow(inventory),'functions in',length(files),'files\n')
# Compact semantic review source for larger files, without prose comments.
large <- c('R/multilpa.R','R/covariates.R','R/covariate-inference.R','R/lta-engine.R','R/transition-inference.R','R/extra-indicators.R','R/fixed-measurement.R','R/prior.R','R/mixture-regression-starts.R')
invisible(lapply(large,function(file){
  p<-parse(file,keep.source=FALSE)
  writeLines(unlist(lapply(p,function(e)deparse(e,width.cutoff=100))),file.path('validation/audit-2026-10-13/model-core',paste0(basename(file),'.static.txt')))
}))
# Explicit review/probe categories added after the synthetic audit.
findings <- list(
 'MC-1'=c('multilpa','.multilpa_covariate_model','.multilpa_fit_covariates','.multilpa_cov_designs'),
 'MC-2'=c('.multilpa_structure_variances','.multilpa_evi_variances','.multilpa_vei_variances','.multilpa_structure_covariances'),
 'MC-3'=c('.growth_class_terms','.growth_residuals'),
 'MC-4'=c('prior_control','.multilpa_prior_maximize','.multilpa_prior_parameters','.multilpa_prior_vev'))
probed <- c(unlist(findings), '.latents_sampling_weights','.growth_statistics','.growth_patterns',
            'mixture_regression','.mixture_spec','.mixture_response')
inventory$static_status <- 'source reviewed; namespace usage checked'
inventory$runtime_status <- ifelse(inventory$function_name %in% probed,'targeted synthetic probe','root integration suite evidence; individual execution not claimed')
inventory$finding_ids <- vapply(inventory$function_name,function(name)paste(names(findings)[vapply(findings,function(group)name %in% group,logical(1))],collapse=';'),character(1))
inventory$review <- NULL
write.csv(inventory,'validation/audit-2026-10-13/model-core/function-inventory.csv',row.names=FALSE)
