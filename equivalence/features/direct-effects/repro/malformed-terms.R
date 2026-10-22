# Malformed `direct_effects` values that are not refused as documented.
# From the project root: Rscript equivalence/features/direct-effects/repro/malformed-terms.R
pkgload::load_all(".", quiet = TRUE)
set.seed(1)
data <- data.frame(z = stats::rnorm(200), y1 = stats::rnorm(200), y2 = stats::rnorm(200),
                   y3 = stats::rnorm(200))
outcome <- function(effects) {
  result <- tryCatch(suppressWarnings(suppressMessages(
    lpa(data, c("y1", "y2", "y3"), 2, n_starts = 1, seed = 1, direct_effects = effects))),
    error = identity)
  if (inherits(result, "error")) paste(class(result), collapse = "/") else
    sprintf("fitted, effects: %s", paste(names(result$direct_effects), collapse = ", "))
}
print(data.frame(direct_effects = c('"y3 ~ z +"', "NA_character_", "1", '"y3 ~ + z"'),
                 outcome = c(outcome("y3 ~ z +"), outcome(NA_character_), outcome(1),
                             outcome("y3 ~ + z"))))
