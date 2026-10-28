inputs <- readRDS('validation/review-fit-search-2026-10-22/benchmark-inputs.rds')
out <- tempfile('latents-archived-pilot-')
dir.create(out)
invisible(lapply(names(inputs),function(batch) {
  invisible(lapply(seq_along(inputs[[batch]]),function(i) {
    saveRDS(inputs[[batch]][[i]],file.path(out,sprintf('%s_%03d.rds',batch,i)))
  }))
}))
log <- system2(file.path(R.home('bin'),'Rscript'),
  c('equivalence/search/analyse-pilot.R',shQuote(out)),stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(log,'status')))
writeLines(log,'validation/spike-fixes-2026-10-22/pilot-metric-coverage-archived.log')
summary <- readRDS(file.path(out,'summary.rds'))
print(subset(summary$metric_counts,Freq>0))
stopifnot(sum(summary$metric_counts$Freq[summary$metric_counts$source=='finalist' & summary$metric_counts$metric_status=='unclassified']) == 1400L,
  !anyNA(summary$rows$status), !any(summary$rows$status=='unclassified'))
print(summary$overall[,c('fit','fits','best','miss','spike','unclassified')])
cat('All archived Gaussian selected fits measured; 1,400 unmeasured Gaussian finalists excluded from reference.\n')
