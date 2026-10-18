# Compare the starts study's latents reference (best of 500 starts) and its
# emulated 10- and 20-start fits with Mplus (STARTS = 400 80) on the
# replicate-1 datasets. Mplus's best final-stage log likelihood is read from
# the ranked list of random-start results; "replicated" means that value was
# reached by more than one final-stage start.
#
# Run from the project root after analyse.R and mplus/run.sh:
#   Rscript equivalence/starts-study/compare-mplus.R

starts_mplus_read <- function(path) {
  lines <- readLines(path, warn = FALSE)
  first <- grep("Final stage loglikelihood values at local maxima", lines)
  ranked <- lines[seq(first + 2L, length(lines))]
  ranked <- ranked[seq_len(which(!nzchar(trimws(ranked)))[1L] - 1L)]
  values <- as.numeric(vapply(strsplit(trimws(ranked), "[[:space:]]+"), `[`, character(1), 1L))
  stem <- sub("[.]out$", "", basename(path))
  data.frame(dataset = sub("-k[0-9]+$", "", stem), k = as.integer(sub(".*-k", "", stem)),
             mplus_best = values[1L],
             mplus_replicated = sum(abs(values - values[1L]) < 1e-3) > 1L,
             mplus_not_replicated_warning = any(grepl("NOT REPLICATED", lines)),
             terminated = any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
}

if (sys.nframe() == 0L) {
  folder <- file.path("equivalence", "starts-study")
  outputs <- list.files(file.path(folder, "mplus"), pattern = "[.]out$", full.names = TRUE)
  # Runs stopped before they finished are reported and left out.
  finished <- vapply(outputs, function(path)
    any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", readLines(path, warn = FALSE))),
    logical(1))
  message(sprintf("unfinished Mplus runs left out: %s",
                  paste(basename(outputs[!finished]), collapse = ", ")))
  mplus <- do.call(rbind, lapply(outputs[finished], starts_mplus_read))
  blocks <- utils::read.csv(file.path(folder, "blocks.csv"))
  emulated <- aggregate(cbind(reached_10 = reached) ~ dataset + k,
                        data = subset(blocks, n_starts == 10L), FUN = mean)
  emulated_20 <- aggregate(cbind(reached_20 = reached) ~ dataset + k,
                           data = subset(blocks, n_starts == 20L), FUN = mean)
  reference <- unique(subset(blocks, select = c(dataset, family, true_k, separation, n, k,
                                                reference)))
  compared <- Reduce(function(a, b) merge(a, b, by = c("dataset", "k")),
                     list(mplus, reference, emulated, emulated_20))
  compared$latents_minus_mplus <- compared$reference - compared$mplus_best
  compared$relative_k <- ifelse(compared$k == compared$true_k, "K", "K+1")
  utils::write.csv(compared, file.path(folder, "mplus-comparison.csv"), row.names = FALSE)
  print(compared, digits = 8)
  print(aggregate(latents_minus_mplus ~ family + relative_k, data = compared,
                  FUN = \(x) c(min = min(x), median = stats::median(x), max = max(x))))
}
