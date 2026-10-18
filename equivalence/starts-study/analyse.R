# Analyse the starts study (see simulate.R for the design).
#
# For every dataset and k the reference is the best log likelihood of the 500
# starts. An s-start fit is emulated by the best of the Ward start and one
# disjoint block of s - 1 k-means starts: 55 blocks for 10 starts, 26 for 20.
# A fit "reached" the reference when it is within 1e-3 of it.
#
# Selection: block b at every k plays the role of one enumeration run, and BIC
# (-2 LL + p log n) picks k among the fitted values; this is compared with the
# choice made from the reference fits and with the true K.
#
# Run from the project root: Rscript equivalence/starts-study/analyse.R

starts_read <- function(dir = file.path("equivalence", "starts-study", "fits")) {
  do.call(rbind, lapply(list.files(dir, pattern = "[.]rds$", full.names = TRUE), readRDS))
}

starts_blocks <- function(fits, n_starts, tolerance = 1e-3) {
  cells <- split(fits, list(fits$dataset, fits$k), drop = TRUE)
  do.call(rbind, lapply(cells, function(cell) {
    cell <- cell[order(cell$start), ]
    ward <- cell$log_likelihood[cell$start == 1L]
    random <- cell$log_likelihood[cell$start > 1L]
    size <- n_starts - 1L
    n_blocks <- length(random) %/% size
    block <- rep(seq_len(n_blocks), each = size)
    best <- vapply(split(random[seq_along(block)], block), max, numeric(1))
    emulated <- pmax(best, ward)
    reference <- max(cell$log_likelihood)
    data.frame(cell[1L, c("dataset", "family", "true_k", "separation", "n", "rep", "k",
                          "n_parameters")],
               n_starts = n_starts, block = seq_len(n_blocks), log_likelihood = emulated,
               reference = reference, gap = emulated - reference,
               reached = emulated - reference > -tolerance,
               ward_reached = ward - reference > -tolerance,
               single_start_rate = mean(cell$log_likelihood - reference > -tolerance),
               row.names = NULL)
  }))
}

starts_selection <- function(blocks) {
  blocks$bic <- -2 * blocks$log_likelihood + blocks$n_parameters * log(blocks$n)
  blocks$reference_bic <- -2 * blocks$reference + blocks$n_parameters * log(blocks$n)
  runs <- split(blocks, list(blocks$dataset, blocks$n_starts, blocks$block), drop = TRUE)
  do.call(rbind, lapply(runs, function(run) {
    reference_k <- run$k[which.min(run$reference_bic)]
    chosen_k <- run$k[which.min(run$bic)]
    sorted <- sort(run$reference_bic)
    data.frame(run[1L, c("dataset", "family", "true_k", "separation", "n", "n_starts",
                         "block")],
               chosen_k = chosen_k, reference_k = reference_k,
               agrees = chosen_k == reference_k,
               direction = sign(chosen_k - reference_k),
               reference_margin = sorted[2L] - sorted[1L], row.names = NULL)
  }))
}

starts_summaries <- function(blocks, selection) {
  blocks$relative_k <- sprintf("K%+d", blocks$k - blocks$true_k)
  blocks$relative_k[blocks$k == blocks$true_k] <- "K"
  reached <- aggregate(cbind(reached, ward_reached) ~ family + relative_k + n_starts,
                       data = blocks, FUN = mean)
  worst <- aggregate(gap ~ family + relative_k + n_starts, data = blocks, FUN = min)
  by_separation <- aggregate(reached ~ family + relative_k + separation + n_starts,
                             data = blocks, FUN = mean)
  by_n <- aggregate(reached ~ family + relative_k + n + n_starts, data = blocks, FUN = mean)
  selection$fewer <- selection$direction < 0
  selection$more <- selection$direction > 0
  selection$true_by_reference <- selection$reference_k == selection$true_k
  selection$true_by_run <- selection$chosen_k == selection$true_k
  choice <- aggregate(cbind(agrees, fewer, more, true_by_reference, true_by_run) ~
                        family + n_starts, data = selection, FUN = mean)
  choice_by_k <- aggregate(cbind(agrees, fewer, more) ~ family + true_k + n_starts,
                           data = selection, FUN = mean)
  list(reached = merge(reached, worst), by_separation = by_separation, by_n = by_n,
       choice = choice, choice_by_k = choice_by_k)
}

if (sys.nframe() == 0L) {
  fits <- starts_read()
  blocks <- rbind(starts_blocks(fits, 10L), starts_blocks(fits, 20L))
  selection <- starts_selection(blocks)
  out <- file.path("equivalence", "starts-study")
  utils::write.csv(blocks, file.path(out, "blocks.csv"), row.names = FALSE)
  utils::write.csv(selection, file.path(out, "selection.csv"), row.names = FALSE)
  summaries <- starts_summaries(blocks, selection)
  message(sprintf("datasets: %d, cells: %d", length(unique(blocks$dataset)),
                  nrow(unique(subset(blocks, select = c(dataset, k))))))
  print(summaries)
  saveRDS(summaries, file.path(out, "summaries.rds"))
}
