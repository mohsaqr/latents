# Compare latents with Latent GOLD 6.1 on the ordinal/count kit. Run from the
# package root after make-kit.R and the Latent GOLD batch runs (listings
# copied to equivalence/latentgold-ordinal/returned/):
#   Rscript equivalence/latentgold-ordinal/compare.R
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-ordinal")
cases <- readRDS(file.path(base, "kit", "cases.rds"))

listing <- function(name) {
  iconv(readLines(file.path(base, "returned", paste0(name, ".lst")), warn = FALSE),
        "latin1", "UTF-8")
}
# The Profile table: class sizes, then per indicator its category
# probabilities (ordinal) or its mean (count, continuous), one value column
# per class (each followed by its standard error).
lg_profile <- function(lines, n_classes, n_group_classes = 0L) {
  lines <- sub("\r$", "", lines)
  start <- which(lines == "Profile")[1L]
  # The table runs to its first blank line, the end of the listing, or the
  # iteration detail that follows it when Latent GOLD printed warnings.
  stop_at <- min(which(seq_along(lines) > start &
                         (lines == "Iteration Detail" | trimws(lines) == "")),
                 length(lines) + 1L)
  rows <- strsplit(lines[(start + 3L):(stop_at - 1L)], "\t")
  label <- vapply(rows, function(row) trimws(row[1L]), character(1))
  # Class value columns are 2, 4, ... (each followed by its standard error);
  # a two-level table prints the group classes' columns first.
  values <- lapply(rows, function(row) {
    columns <- 2L * n_group_classes + 1L + 2L * seq_len(n_classes) - 1L
    cells <- row[columns]
    if (length(row) < max(columns)) return(NULL)
    as.numeric(cells)
  })
  header <- vapply(values, is.null, logical(1))
  data_rows <- which(!header)[-1L]
  out <- list(sizes = values[[1L]])
  indicator_of <- label[header][cumsum(header)[data_rows]]
  invisible(lapply(unique(indicator_of), function(name) {
    taken <- data_rows[indicator_of == name]
    out[[name]] <<- stats::setNames(values[taken], label[taken])
  }))
  out
}

compare_case <- function(name) {
  case <- cases[[name]]
  data <- utils::read.delim(file.path(base, "kit", paste0(name, ".dat")))
  types <- vapply(case$indicators, `[[`, character(1), "type")
  vars <- names(case$indicators)
  two_level <- !is.null(case$group_size)
  fit <- withCallingHandlers(
    multilpa(data, vars, if (two_level) "g" else NULL, length(case$sizes),
             n_group_classes = if (two_level) length(case$group_classes) else 1L,
             ordinal = vars[types == "ordinal"],
             count = vars[types %in% c("count", "negbin")],
             count_model = if (any(types == "negbin")) "negative_binomial" else "poisson",
             count_dispersion = case$dispersion %||% "varying",
             # tol 1e-10: a weakly separated negative-binomial case (o07)
             # converges at EM's linear rate and does not reach 1e-12.
             n_starts = 30, seed = 1, tol = 1e-10, max_iter = 5000),
    latents_single_level = function(notice) invokeRestart("muffleMessage"))
  lines <- listing(name)
  lg_ll <- as.numeric(strsplit(grep("^Log-likelihood \\(LL\\)", lines, value = TRUE)[1L],
                               "\t")[[1L]][2L])
  profile <- lg_profile(lines, length(case$sizes),
                        if (two_level) length(case$group_classes) else 0L)
  # Ours, in the same layout, for every relabelling of the classes.
  # Two-level class sizes are shares within group classes, which Latent GOLD
  # prints differently; there the measurement parameters are compared.
  ours <- function(order) {
    c(if (!two_level) fit$profile_probabilities[1L, order],
      unlist(lapply(vars, function(v) {
        switch(types[[v]],
          continuous = fit$means[order, v],
          count = fit$count_means[order, v],
          negbin = fit$count_means[order, v],
          ordinal = {
            j <- match(v, colnames(fit$extra_data$ordinal))
            exp(.latents_ordinal_log_probabilities(fit$ordinal_intercepts[[j]],
                                                   fit$ordinal_locations[, j]))[order, ]
          })
      })))
  }
  theirs <- c(if (!two_level) profile$sizes, unlist(lapply(vars, function(v) {
    block <- profile[[v]]
    if (types[[v]] == "ordinal") {
      do.call(cbind, block[setdiff(names(block), "Mean")])
    } else block[["Mean"]]
  })))
  orders <- if (length(case$sizes) == 2L) list(1:2, 2:1) else
    list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  gaps <- vapply(orders, function(order) max(abs(ours(order) - theirs)), numeric(1))
  data.frame(case = name, latents_ll = fit$log_likelihood, latent_gold_ll = lg_ll,
             ll_difference = fit$log_likelihood - lg_ll,
             n_parameters = fit$n_parameters,
             max_parameter_difference = min(gaps))
}
comparison <- do.call(rbind, lapply(names(cases), compare_case))
print(comparison, digits = 10, row.names = FALSE)
utils::write.csv(comparison, file.path(base, "comparison.csv"), row.names = FALSE)
