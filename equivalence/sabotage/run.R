# Sabotage battery (Equivalence 2.0, E2.7).
#
# A test that cannot fail proves nothing. Each defect below is a small,
# realistic bug injected into the loaded package, in memory only (the source
# files are not touched): a wrong denominator, a transposed matrix, a row or
# group counted that should not be, a dropped weight. With the defect in
# place, the tests named for it are run; the defect counts as detected when
# at least one expectation fails. A defect that survives marks a missing test.
# Each filter names the tests that exercise the target in this process: the
# scripts run by test-external-scripts.R start their own R and never see an
# injected defect, so no filter relies on them.
#
# Run from the project root:
#   Rscript equivalence/sabotage/run.R            # every defect
#   Rscript equivalence/sabotage/run.R entropy    # defects whose name matches
# Writes equivalence/sabotage/results.csv.

#' Each defect: the internal function it replaces, how (a function of the
#' original), the suite and test-file filter that should catch it.
sabotage_defects <- function() {
  list(
    # The control: the target replaced by itself. It must not be detected; if
    # it is, the tests fail without any defect and no detection below counts.
    list(name = "control_no_defect",
         description = "No defect (the original function)",
         target = ".multilpa_relative_entropy",
         replace = function(original) original,
         suite = "equivalence",
         filter = "missing-patterns|mplus-mixture-weights|mplus-standardize|mplus-person-scaling"),
    list(name = "entropy_denominator",
         description = "Relative entropy divides by 1.01 n log K",
         target = ".multilpa_relative_entropy",
         replace = function(original) function(probabilities) {
           if (ncol(probabilities) == 1L) return(NA_real_)
           1 - .multilpa_entropy_sum(probabilities) /
             (nrow(probabilities) * log(ncol(probabilities)) * 1.01)
         },
         suite = "equivalence", filter = "montecarlo|mplus-classification"),
    list(name = "empty_rows_counted",
         description = "Rows with nothing observed counted in classification statistics",
         target = ".multilpa_classified_posteriors",
         replace = function(original) function(x) .multilpa_class_posteriors(x),
         suite = "equivalence", filter = "missing-patterns|montecarlo"),
    list(name = "empty_groups_counted",
         description = "Groups with nothing observed counted in group statistics",
         target = ".multilpa_informative_groups",
         replace = function(original) function(x)
           if (is.null(x$group_posteriors)) NULL else rep(TRUE, nrow(x$group_posteriors)),
         suite = "equivalence", filter = "missing-patterns"),
    list(name = "avepp_transposed",
         description = "Average posterior matrix transposed",
         target = ".multilpa_average_posterior_matrix",
         replace = function(original) function(probabilities) t(original(probabilities)),
         suite = "equivalence", filter = "montecarlo|mplus-classification"),
    list(name = "parameter_count",
         description = "One parameter too many",
         target = ".multilpa_count_parameters",
         replace = function(original) function(...) original(...) + 1L,
         suite = "equivalence", filter = "montecarlo|mclust"),
    list(name = "gaussian_density",
         description = "Gaussian log density scaled by 1.001",
         target = ".latents_gaussian_log_density",
         replace = function(original) function(...) original(...) * 1.001,
         suite = "equivalence", filter = "montecarlo|mclust"),
    list(name = "variance_mstep",
         description = "M-step variances inflated by 1%",
         target = ".multilpa_maximization",
         replace = function(original) function(...) {
           parameters <- original(...)
           if (!is.null(parameters$variances)) parameters$variances <- parameters$variances * 1.01
           parameters
         },
         suite = "equivalence", filter = "montecarlo"),
    list(name = "mixing_update",
         description = "Mixing (class share) update shrunk toward equal shares",
         target = ".mixture_update_mixing",
         # The update returns the whole parameter list; its class logits are
         # `gamma` (one level) or `delta` and `class_logits` (two levels).
         replace = function(original) function(spec, params, expectation) {
           updated <- original(spec, params, expectation)
           shrink <- intersect(c("gamma", "delta", "class_logits"), names(updated))
           updated[shrink] <- lapply(updated[shrink], \(logits) logits * 0.99)
           updated
         },
         suite = "equivalence", filter = "mixture-regression-flexmix|growth-mixture-lcmm"),
    list(name = "weights_dropped",
         description = "Sampling weights ignored in the E-step",
         target = ".multilpa_weigh",
         replace = function(original) function(expectation, weights, group_index) expectation,
         suite = "equivalence", filter = "mplus-mixture-weights"),
    list(name = "scaling_denominator",
         description = "Standard deviations for scaling use n instead of n - 1",
         target = ".multilpa_indicator_scales",
         replace = function(original) function(x, group_index, n_groups, scaling) {
           scales <- original(x, group_index, n_groups, scaling)
           if (is.null(scales)) return(scales)
           counts <- if (identical(scaling, "grand")) matrix(colSums(!is.na(x)), 1L) else
             rowsum((!is.na(x)) * 1, group_index, reorder = TRUE)
           scales * sqrt(pmax(counts - 1, 1) / pmax(counts, 1))
         },
         suite = "equivalence", filter = "mplus-standardize|mplus-person-scaling"),
    list(name = "blrt_p_value",
         description = "Bootstrap p-value without the +1 correction",
         target = ".latents_blrt_p_value",
         replace = function(original) function(replicates, observed, iter) {
           tally <- original(replicates, observed, iter)
           tally$p_value <- max(0, (tally$p_value * (iter + 1) - 1) / iter)
           tally
         },
         # Against TECH14 the p-values agree only within Monte Carlo error,
         # which this defect stays inside; the exact formula is a unit test.
         suite = "unit", filter = "^enumeration$"),
    list(name = "lrt_degrees_of_freedom",
         description = "Likelihood-ratio test with one degree of freedom too many",
         target = ".latents_lrt_pair",
         replace = function(original) function(...) {
           test <- original(...)
           test$df <- test$df + 1L
           test
         },
         suite = "equivalence", filter = "^lrt$"),
    list(name = "lrt_covariance_nesting",
         description = "Every covariance structure taken as nested in every other",
         target = ".multilpa_lrt_covariance_nested",
         replace = function(original) function(null, alternative) TRUE,
         suite = "unit", filter = "^compare-lrt$")
  )
}

#' Inject one defect, run its tests, restore the original
run_defect <- function(defect, root) {
  namespace <- asNamespace("latents")
  original <- get(defect$target, envir = namespace)
  replacement <- defect$replace(original)
  # The replacement must see the package's internals and still reach
  # `original`: an environment holding `original`, inside the namespace.
  # (Setting the namespace itself would drop `original`, every wrapped call
  # would error, and the error would pass for a detection.)
  scope <- new.env(parent = namespace)
  scope$original <- original
  environment(replacement) <- scope
  utils::assignInNamespace(defect$target, replacement, ns = "latents")
  on.exit(utils::assignInNamespace(defect$target, original, ns = "latents"), add = TRUE)
  directory <- if (identical(defect$suite, "unit")) file.path(root, "tests", "testthat") else
    file.path(root, "equivalence")
  started <- Sys.time()
  results <- tryCatch(as.data.frame(testthat::test_dir(
    directory, filter = defect$filter, reporter = "silent", stop_on_failure = FALSE,
    package = "latents", load_package = "none")),
    error = function(condition) data.frame(failed = 1L, error = TRUE, nb = 0L,
                                           note = conditionMessage(condition)))
  data.frame(defect = defect$name, description = defect$description,
             target = defect$target, filter = defect$filter,
             expectations = sum(results$nb), failures = sum(results$failed),
             errors = sum(results$error),
             detected = sum(results$failed) + sum(results$error) > 0L,
             minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1))
}

Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
root <- getwd()
requested <- commandArgs(trailingOnly = TRUE)
defects <- sabotage_defects()
if (length(requested)) {
  defects <- Filter(\(defect) grepl(paste(requested, collapse = "|"), defect$name), defects)
}
results <- do.call(rbind, lapply(defects, run_defect, root = root))
utils::write.csv(results, file.path(root, "equivalence", "sabotage", "results.csv"),
                 row.names = FALSE)
print(subset(results, TRUE, c(defect, failures, errors, detected, minutes)))
control <- subset(results, defect == "control_no_defect")
if (nrow(control) && isTRUE(control$detected)) {
  stop("the control run failed without a defect; the detections are not valid", call. = FALSE)
}
