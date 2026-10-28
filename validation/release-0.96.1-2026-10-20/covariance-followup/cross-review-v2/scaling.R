pkgload::load_all(quiet = TRUE)
scatter <- rbind(c(1e-5, 1), c(1e-4, 100))
reference <- matrix(c(1, 1, 50.5, 50.5), 2)
settings <- expand.grid(units = c(1e-100, 1e-50, 1, 1e50, 1e100),
                        counts = c(1e-100, 1, 1e100))
result <- lapply(seq_len(nrow(settings)), function(index) {
  unit <- settings$units[index]
  count <- settings$counts[index]
  actual <- .multilpa_bounded_variances(scatter * unit * count,
    c(1, 1) * count, "EVI", unit)
  relative_error <- max(abs(actual / unit - reference)) / max(reference)
  stopifnot(relative_error < 1e-10)
  data.frame(units = unit, counts = count, relative_error = relative_error)
})
print(do.call(rbind, result), digits = 16)
