library(latents)
z <- readRDS('/tmp/astra-lta-bootstrap-probe.rds')
f <- z$f
keep <- rowSums(!is.na(f$x)) > 0
p <- f$subject_posteriors[keep, , drop = FALSE]
checks <- get_results(diagnostics(f, plots = FALSE), 'checks')
cat('n_observations:',f$n_observations,'n_informative:',f$n_informative,'drop_empty:',f$drop_empty,'\n')
cat('informative expected class counts:',colSums(p),'\n')
cat('diagnostics expected class counts:',latents:::.latents_level_classification(f,'individuals')$estimated_n,'\n')
expected_entropy <- 1 + sum(ifelse(p>0,p*log(p),0))/(nrow(p)*log(ncol(p)))
cat('direct informative-row entropy:',expected_entropy,'reported:',checks$value[checks$check=='relative entropy'],'\n')
dg <- diagnostics(f,plots=FALSE)
cat('as.data.frame columns:',names(as.data.frame(dg)),'\n')
print(tryCatch(get_results(dg,'entropy'),error=conditionMessage))
