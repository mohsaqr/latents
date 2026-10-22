# Break the VLMR implementation in memory, never on disk, and confirm the
# unit tests catch it. Run from the project root:
#   NOT_CRAN=true Rscript equivalence/features/vlmr/sabotage.R
pkgload::load_all(".", quiet = TRUE)
namespace <- asNamespace("latents")
attached <- as.environment("package:latents")

replace_function <- function(name, replacement) {
  original <- get(name, envir = namespace)
  environment(replacement) <- namespace
  invisible(lapply(list(namespace, attached), function(env) {
    if (bindingIsLocked(name, env)) unlockBinding(name, env)
    assign(name, replacement, envir = env)
  }))
  original
}

run_tests <- function() {
  results <- as.data.frame(testthat::test_file(
    file.path("tests", "testthat", "test-features-vlmr.R"), reporter = "silent",
    stop_on_failure = FALSE))
  c(failed = sum(results$failed), passed = sum(results$passed), errors = sum(results$error))
}

sabotages <- list(
  # The chi-square limit in place of the weighted reference distribution.
  chi_square_tail = list(name = ".latents_imhof_upper", replacement = function(q, lambda) {
    stats::pchisq(q, sum(lambda > 0), lower.tail = FALSE)
  }),
  # Vuong's sign matrix dropped: both models' whitened scores added.
  unsigned_weights = list(name = ".latents_vuong_weights", replacement = local({
    body <- body(.latents_vuong_weights)
    text <- deparse(body, width.cutoff = 500L)
    text <- sub("rep(-1, ncol(null$scores))", "rep(1, ncol(null$scores))", text, fixed = TRUE)
    stopifnot("the sabotage must change the body" =
                !identical(text, deparse(body, width.cutoff = 500L)))
    f <- function(alternative, null) NULL
    body(f) <- parse(text = paste(text, collapse = "\n"))[[1L]]
    f
  })),
  # The adjustment with n in place of log n.
  adjustment = list(name = ".latents_vuong_row", replacement = local({
    text <- deparse(body(.latents_vuong_row), width.cutoff = 500L)
    changed <- sub("added * log(n_units)", "added * n_units", text, fixed = TRUE)
    stopifnot("the sabotage must change the body" = !identical(changed, text))
    f <- function(alternative, null) NULL
    body(f) <- parse(text = paste(changed, collapse = "\n"))[[1L]]
    f
  })))

baseline <- run_tests()
outcome <- lapply(names(sabotages), function(label) {
  sabotage <- sabotages[[label]]
  original <- replace_function(sabotage$name, sabotage$replacement)
  on.exit(replace_function(sabotage$name, original), add = TRUE)
  c(sabotage = label, as.list(run_tests()))
})
print(rbind(data.frame(sabotage = "none", as.list(baseline)),
            do.call(rbind, lapply(outcome, as.data.frame))))
