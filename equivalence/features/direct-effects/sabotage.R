# Break the direct-effect implementation in memory, one defect at a time, and
# check that tests/testthat/test-features-direct-effects.R notices. No source
# file is changed: each defect replaces one function in the loaded namespace
# and the original is put back afterwards.
#
# From the project root:
#   NOT_CRAN=true Rscript equivalence/features/direct-effects/sabotage.R
# Writes sabotage.csv next to this file. The unit file has three known
# failures (FINDINGS.md); a defect counts as detected when it fails more.
pkgload::load_all(".", quiet = TRUE)
Sys.setenv(NOT_CRAN = "true")
namespace <- asNamespace("latents")
unit_file <- file.path("tests", "testthat", "test-features-direct-effects.R")

#' A copy of a namespace function with one expression substituted in its body
rewritten <- function(name, edit) {
  original <- get(name, envir = namespace)
  changed <- original
  body(changed) <- edit(body(original))
  stopifnot("the edit did not change the function" = !identical(body(changed), body(original)))
  changed
}

substitute_call <- function(expr, from, to) {
  if (identical(expr, from)) return(to)
  if (is.call(expr)) {
    replaced <- lapply(as.list(expr), \(part) if (is.null(part)) part else
      substitute_call(part, from, to))
    return(as.call(replaced))
  }
  expr
}

defects <- list(
  control = list(name = ".latents_direct_shift", edit = identity),
  # The effect is subtracted at half its size in the likelihood.
  half_shift = list(name = ".latents_direct_shift", edit = \(b) substitute_call(
    b, quote(x - tcrossprod(direct$covariates, .latents_direct_matrix(direct,
      parameters$direct_effects, ncol(x)))),
    quote(x - 0.5 * tcrossprod(direct$covariates, .latents_direct_matrix(direct,
      parameters$direct_effects, ncol(x)))))),
  # The M-step ignores the class variances (unweighted pooled regression).
  unweighted_step = list(name = ".latents_direct_maximize", edit = \(b) substitute_call(
    b, quote(sweep(posteriors[observed, , drop = FALSE], 2L, variances[, column], "/")),
    quote(posteriors[observed, , drop = FALSE]))),
  # Means are reported at the covariates' centre instead of at zero.
  no_uncentering = list(name = ".latents_uncenter_direct", edit = \(b) {
    as.call(c(as.name("{"), quote(return(parameters))))
  }),
  # The score of an effect forgets its covariate.
  score_without_covariate = list(name = ".latents_direct_scores", edit = \(b) substitute_call(
    b, quote(direct$covariates[, effects$covariate_column[[row]]]), quote(1))))

failures_of <- function() {
  results <- as.data.frame(testthat::test_file(unit_file, reporter = "silent",
                                               package = "latents", load_package = "none"))
  c(failed = sum(results$failed), errors = sum(results$error),
    passed = sum(results$passed))
}

rows <- lapply(names(defects), function(label) {
  defect <- defects[[label]]
  original <- get(defect$name, envir = namespace)
  broken <- if (identical(defect$edit, identity)) original else
    rewritten(defect$name, defect$edit)
  utils::assignInNamespace(defect$name, broken, ns = "latents")
  on.exit(utils::assignInNamespace(defect$name, original, ns = "latents"), add = TRUE)
  counts <- failures_of()
  data.frame(defect = label, target = defect$name, failed = counts[["failed"]],
             errors = counts[["errors"]], passed = counts[["passed"]])
})
table <- do.call(rbind, rows)
control <- subset(table, defect == "control")
table$detected <- table$failed + table$errors > control$failed + control$errors
print(table)
stopifnot("the control changed the result" = !any(subset(table, defect == "control")$detected))
utils::write.csv(table, file.path("equivalence", "features", "direct-effects", "sabotage.csv"),
                 row.names = FALSE)
