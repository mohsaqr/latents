# Break the ML three-step in memory, one defect per process, and confirm the
# unit tests notice. Source files are never touched.
#
# From the project root:
#   Rscript equivalence/features/ml-three-step/sabotage.R            # every defect
#   Rscript equivalence/features/ml-three-step/sabotage.R transposed  # one
sabotages <- list(
  transposed = list(name = ".multilpa_error_matrix",
                    make = \(original) function(pieces) t(original(pieces))),
  no_cr1 = list(name = ".inference_cluster_factor", make = \(original) function(n) 1),
  inflated_se = list(name = ".multilpa_step_ml_influence",
                     make = \(original) function(...) 1.05 * original(...)),
  equal_ignored = list(name = ".multilpa_step_ml_fit",
                       make = \(original) function(pieces, values, variance, ...)
                         original(pieces, values, "varying", ...)),
  no_polish = list(name = ".multilpa_step_ml_polish",
                   make = \(original) function(fit, values, variance) fit))
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 1L && requested %in% names(sabotages)) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  sabotage <- sabotages[[requested]]
  namespace <- asNamespace("latents")
  # Evaluated before the replacement: a lazy argument would read the broken
  # function itself and recurse.
  original <- get(sabotage$name, envir = namespace)
  force(original)
  broken <- sabotage$make(original)
  utils::assignInNamespace(sabotage$name, broken, ns = "latents")
  # load_all() also exported the original into the package environment.
  exported <- as.environment("package:latents")
  if (exists(sabotage$name, envir = exported, inherits = FALSE)) {
    unlockBinding(sabotage$name, exported)
    assign(sabotage$name, broken, envir = exported)
  }
  stopifnot("the sabotage did not replace the function" =
              !identical(get(sabotage$name, envir = namespace), original))
  Sys.setenv(NOT_CRAN = "true")
  result <- as.data.frame(testthat::test_file(
    file.path("tests", "testthat", "test-features-ml-three-step.R"),
    reporter = if (nzchar(Sys.getenv("SABOTAGE_VERBOSE"))) "summary" else "silent",
    stop_on_failure = FALSE))
  failed <- result$test[result$failed > 0L | result$error]
  cat(sprintf("RESULT %s: %d of %d tests failed%s\n", requested, length(failed), nrow(result),
              if (length(failed)) paste0(": ", paste(failed, collapse = " | ")) else ""))
} else {
  outcome <- vapply(names(sabotages), \(name) {
    lines <- system2(file.path(R.home("bin"), "Rscript"),
                     c(file.path("equivalence", "features", "ml-three-step", "sabotage.R"), name),
                     stdout = TRUE, stderr = TRUE)
    line <- grep("^RESULT", lines, value = TRUE)
    if (length(line)) line else paste("RESULT", name, "did not run:", utils::tail(lines, 1L))
  }, character(1))
  writeLines(outcome)
  caught <- grepl(": [1-9][0-9]* of", outcome)
  if (!all(caught)) stop("a sabotage went unnoticed", call. = FALSE)
}
