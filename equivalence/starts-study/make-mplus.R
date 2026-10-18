# Mplus reference for the starts study: replicate 1 of every condition, fitted
# at the true K and at K + 1 with STARTS = 400 80 (as in the enumeration kit).
# This is a reference from a different start scheme, so the study's own
# reference (the best of 500 latents starts) is not judged only against itself.
#
# Run from the project root:
#   Rscript equivalence/starts-study/make-mplus.R
#   sh equivalence/starts-study/mplus/run.sh

source(file.path("equivalence", "starts-study", "simulate.R"))

starts_mplus_input <- function(row, k) {
  vars <- if (identical(row$family, "lpa")) paste0("y", 1:4) else paste0("u", 1:6)
  stem <- sprintf("%s-k%d", row$dataset, k)
  model <- if (identical(row$family, "lpa")) {
    c(" %OVERALL%", unlist(lapply(seq_len(k), function(class)
      c(sprintf(" %%c#%d%%", class), sprintf(" [%s];", paste(vars, collapse = " ")),
        sprintf(" %s;", paste(vars, collapse = " "))))))
  } else " %OVERALL%"
  c(sprintf("TITLE: latents starts study %s, %d classes;", row$dataset, k),
    sprintf("DATA: FILE = %s.dat;", row$dataset),
    "VARIABLE:",
    sprintf(" NAMES = id %s;", paste(vars, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", paste(vars, collapse = " ")),
    " IDVARIABLE = id;",
    if (identical(row$family, "lca")) sprintf(" CATEGORICAL = %s;", paste(vars, collapse = " ")),
    sprintf(" CLASSES = c(%d);", k),
    "ANALYSIS:", " TYPE = MIXTURE;", " ESTIMATOR = ML;", " STARTS = 400 80;",
    " STITERATIONS = 40;", " STSEED = 20261018;", " MITERATIONS = 10000;",
    " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
    " MCONVERGENCE = 0.000001;", "MODEL:", model,
    "SAVEDATA:", sprintf(" RESULTS = %s-results.dat;", stem))
}

if (sys.nframe() == 0L) {
  out <- file.path("equivalence", "starts-study", "mplus")
  dir.create(out, showWarnings = FALSE)
  chosen <- subset(starts_design(), rep == 1L)
  invisible(lapply(seq_len(nrow(chosen)), function(i) {
    row <- chosen[i, ]
    data <- starts_simulate(row)$data
    columns <- c("id", setdiff(names(data), "id"))
    utils::write.table(subset(data, select = columns), file.path(out, paste0(row$dataset, ".dat")),
                       row.names = FALSE, col.names = FALSE)
    lapply(row$true_k + 0:1, function(k)
      writeLines(starts_mplus_input(row, k), file.path(out, sprintf("%s-k%d.inp", row$dataset, k))))
  }))
  message(sprintf("%d inputs written", length(list.files(out, pattern = "[.]inp$"))))
}
