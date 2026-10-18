# Enumeration against Mplus: two Monte Carlo datasets fitted by Mplus at every
# number of classes, so that an enumeration grid, its minima and the model
# comparison can be checked candidate by candidate.
#
# Run from the project root:
#   Rscript equivalence/montecarlo/enumeration/make-kit.R
#   sh equivalence/montecarlo/enumeration/run.sh

directory <- file.path("equivalence", "montecarlo", "enumeration")
cases <- readRDS(file.path("equivalence", "montecarlo", "kit", "cases.rds"))
designs <- list(lpa_vvi_01 = 1:5, lca_binary_01 = 1:4)

#' The Mplus input for one dataset at k classes
enumeration_input <- function(name, case, k) {
  spec <- case$spec
  vars <- paste(spec$vars, collapse = " ")
  class_block <- if (isTRUE(spec$categorical)) character() else
    unlist(lapply(seq_len(k), \(j) c(sprintf(" %%c#%d%%", j), sprintf(" [%s];", vars),
                                     sprintf(" %s;", vars))))
  c(sprintf("TITLE: latents enumeration %s, %d classes;", name, k),
    sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(names(case$data), collapse = " ")),
    sprintf(" USEVARIABLES = %s;", vars), " IDVARIABLE = id;",
    if (isTRUE(spec$categorical)) sprintf(" CATEGORICAL = %s;", vars),
    sprintf(" CLASSES = c(%d);", k),
    "ANALYSIS:", " TYPE = MIXTURE;", " ESTIMATOR = ML;", " STARTS = 400 80;",
    " STITERATIONS = 40;", " STSEED = 20261018;", " MITERATIONS = 10000;",
    " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
    " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", class_block,
    "SAVEDATA:", sprintf(" RESULTS = %s-k%d-results.dat;", name, k))
}

invisible(lapply(names(designs), function(name) {
  case <- cases[[name]]
  file.copy(file.path("equivalence", "montecarlo", "kit", paste0(name, ".dat")),
            directory, overwrite = TRUE)
  lapply(designs[[name]], function(k)
    writeLines(enumeration_input(name, case, k),
               file.path(directory, sprintf("%s-k%d.inp", name, k))))
}))
saveRDS(designs, file.path(directory, "designs.rds"))
cat("written", sum(lengths(designs)), "inputs\n")
