# Where does the three-class ordinal regression stop on this platform?
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("equivalence/argument-sweep/helpers.R")
source("equivalence/argument-sweep/datasets.R")
data <- sweep_regression()
# The sweep runs its cases under this generator; the starts differ from Mersenne-Twister's.
RNGkind("L'Ecuyer-CMRG")
outcome <- withCallingHandlers(
  tryCatch(mixture_regression(grade ~ x, data, n_classes = 3L, family = "ordinal", seed = 81L,
                              n_starts = 10L, tol = 1e-10, max_iter = 5000L),
           error = function(e) e),
  error = function(e) {
    calls <- vapply(sys.calls(), function(call) paste(deparse(call)[[1L]], collapse = ""), "")
    cat(sprintf("%2d: %s\n", seq_along(calls), substr(calls, 1L, 110L)), sep = "")
  },
  warning = function(w) {
    cat("warning:", class(w)[[1L]], substr(conditionMessage(w), 1L, 120L), "\n")
    invokeRestart("muffleWarning")
  })
cat("outcome:", paste(class(outcome), collapse = "/"), "\n")
if (inherits(outcome, "error")) cat(conditionMessage(outcome), "\n")
print(utils::sessionInfo()$BLAS); print(utils::sessionInfo()$LAPACK)
tables <- c("coefficients", "classes", "membership", "fit", "assignments", "fitted", "starts",
            "classification")
invisible(lapply(tables, function(name) {
  withCallingHandlers(
    tryCatch({
      suppressMessages(get_results(outcome, name))
      cat(name, ": ok\n")
    }, error = function(e) cat(name, ": ERROR", conditionMessage(e), "\n")),
    error = function(e) {
      calls <- vapply(sys.calls(), function(call) paste(deparse(call)[[1L]], collapse = ""), "")
      cat(sprintf("   %2d: %s\n", seq_along(calls), substr(calls, 1L, 120L)), sep = "")
    },
    warning = function(w) invokeRestart("muffleWarning"))
}))
