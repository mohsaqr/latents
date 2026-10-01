# Shared helpers for the latent transition extension comparisons. Run from
# the project root; JStats is expected beside this repository.
devtools::load_all(".", quiet = TRUE)
jstats <- file.path("..", "JStats")
stopifnot("JStats must sit beside latents" = dir.exists(jstats))

permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  smaller <- permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], nrow(smaller)))
  }))
}

# The permutation `perm` with ours[perm[m], ] closest to reference[m, ].
align <- function(ours, reference) {
  options <- permutations(nrow(reference))
  loss <- vapply(seq_len(nrow(options)), function(k) {
    sum((ours[options[k, ], , drop = FALSE] - reference)^2)
  }, numeric(1))
  options[which.min(loss), ]
}

# Wide panel (subject x occasion x item array, or columns) to long rows.
to_long <- function(array3, extra = NULL) {
  dims <- dim(array3)
  do.call(rbind, lapply(seq_len(dims[2L]), function(t) {
    items <- as.data.frame(matrix(array3[, t, ], dims[1L]))
    names(items) <- paste0("y", seq_len(dims[3L]))
    base <- data.frame(subject = seq_len(dims[1L]), occasion = t)
    if (!is.null(extra)) base <- cbind(base, extra)
    cbind(base, items)
  }))
}

row_result <- function(case, quantity, ours, reference, tolerance) {
  data.frame(case = case, quantity = quantity, latents = ours,
             reference = reference, difference = abs(ours - reference),
             tolerance = tolerance, agree = abs(ours - reference) <= tolerance,
             stringsAsFactors = FALSE)
}
