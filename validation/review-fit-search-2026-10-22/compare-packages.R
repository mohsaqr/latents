# Recompute the saved matched comparison; no model fitting is performed.
output <- file.path("validation", "review-fit-search-2026-10-22")
mclust_rows <- readRDS(file.path(output, "mclust-benchmark-inputs.rds"))
latents_rows <- utils::read.csv(file.path(output, "recomputed.csv"))
latents_rows <- subset(latents_rows, condition %in% unique(mclust_rows$condition))
stopifnot(nrow(mclust_rows) == 160L, nrow(latents_rows) == 160L,
          all(is.finite(mclust_rows$ratio)),
          all(is.finite(mclust_rows$mclust_ll)))
mclust_summary <- aggregate(
  cbind(spikes = ratio < .01, fits = rep(1L, nrow(mclust_rows))) ~ prior,
  mclust_rows, sum)
latents_summary <- aggregate(
  cbind(spikes = status == "spike", fits = rep(1L, nrow(latents_rows))) ~ fit,
  latents_rows, sum)
stopifnot(nrow(mclust_summary) == 2L, nrow(latents_summary) == 2L,
          all(mclust_summary$spikes == c(1L, 0L)),
          all(latents_summary$spikes == c(25L, 12L)),
          all(mclust_summary$fits == 80L), all(latents_summary$fits == 80L))
print(mclust_summary)
print(latents_summary)
reference <- unique(latents_rows[c("condition", "replicate", "k", "best")])
matched <- merge(mclust_rows, reference, by = c("condition", "replicate", "k"))
matched$gap <- matched$best - matched$mclust_ll
stopifnot(nrow(matched) == 160L, all(is.finite(matched$gap)))
print(aggregate(gap ~ condition + k + prior, matched, median))
utils::write.csv(matched, file.path(output, "mclust-comparison.csv"),
                 row.names = FALSE)
