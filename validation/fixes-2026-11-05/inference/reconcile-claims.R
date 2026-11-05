claims_path <- "equivalence/CLAIMS.csv"
before <- utils::read.csv(claims_path, stringsAsFactors = FALSE)
stopifnot(nrow(before) == 44L)
changes <- c(
  `9` = '`"robust"` is the sandwich over the',
  `10` = 'independent units (groups), with the `G / (G - 1)` small-sample factor;',
  `11` = 'reports every other parameter conditional on them, as Mplus does. The',
  `12` = '`scaling_correction`, the MLR factor used by scaled likelihood-ratio',
  `13` = 'default of `sandwich::vcovCL()`. Without it the',
  `18` = '`bic_rows` on the rows)',
  `19` = 'in no sample size and no entropy.',
  `21` = 'combinations are the fourteen mclust models:',
  `22` = 'covariates; the mixing proportions have no prior. As in mclust, the',
  `23` = "The component's density is `1 / V`, with `V`",
  `24` = "mclust's `hypvol()`, and it adds two parameters; a supplied `start`",
  `25` = 'the classification: `subject_profiles` is `0` there, the assignment',
  `26` = 'and classification tables, as Mplus does, keeping their rows with `NA`',
  `27` = "same indicators, each group's indicator means being between-level",
  `28` = 'positive weight as Mplus does, so integer weights act as repeated units.',
  `29` = "by all profiles and one location per profile, Latent GOLD's default",
  `30` = 'variance `mu + alpha * mu^2`) for counts more dispersed than a Poisson',
  `31` = "unit's likelihood is the product of its rows' densities.",
  `32` = 'matrix, `lcmm::hlme(nwg = FALSE)`) or `"proportional"` (a shared matrix',
  `33` = '`TRUE` for uncorrelated random effects (`idiag =',
  `34` = 'autoregressive process (covariance `sd^2 * exp(-rate * |t - s|)`) and',
  `36` = 'The prior follows `mclust::priorControl()`.',
  `44` = "posteriors, as in Latent GOLD's step-three model. The likelihood is"
)
after <- before
after$phrase[as.integer(names(changes))] <- unname(changes)
present <- vapply(seq_len(nrow(after)), function(i) {
  lines <- readLines(after$source[i], warn = FALSE)
  if (startsWith(after$source[i], "R/")) lines <- sub("^#' ?", "", lines)
  any(grepl(after$phrase[i], lines, fixed = TRUE))
}, logical(1))
stopifnot(all(present), nrow(after) == 44L,
          identical(after[setdiff(names(after), "phrase")],
                    before[setdiff(names(before), "phrase")]),
          sum(after$test_file == "NONE") == 1L,
          !anyDuplicated(after[c("source", "phrase", "reference")]))
utils::write.csv(after, claims_path, row.names = FALSE, na = "")
cat("Updated", sum(before$phrase != after$phrase),
    "phrases; retained all", nrow(after), "rows and all reference/test metadata.\n")
print(data.frame(row = as.integer(names(changes)), source = after$source[as.integer(names(changes))],
                 phrase = unname(changes)), row.names = FALSE)
