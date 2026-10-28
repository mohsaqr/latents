pkgload::load_all(quiet = TRUE)
output <- file.path("validation", "spike-fixes-2026-10-22")
before <- readRDS(file.path("validation", "review-fit-search-2026-10-22",
                           "spike-gaps-probes.rds"))$known_flat$fit
set.seed(741)
stream <- .Random.seed
after <- lpa(srl, names(srl), 4, model = "VVV", n_starts = 10, seed = 1)
stopifnot(identical(.Random.seed, stream))
fields <- names(before)
equal <- vapply(fields, function(field) identical(before[[field]], after[[field]]), logical(1))
print(data.frame(field = fields, identical = equal))
stopifnot(all(equal))
print(after)
print(get_results(after, "gaussian_diagnostics"))
utils::write.csv(get_results(after, "gaussian_diagnostics"),
                 file.path(output, "diagnostic-example.csv"), row.names = FALSE)
saveRDS(after, file.path(output, "diagnostic-example.rds"))
cat("Every preceding field and the caller's RNG stream are identical.\n")
