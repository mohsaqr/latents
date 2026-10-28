output <- file.path("validation", "release-0.96.1-2026-10-20", "claims")
registry <- "equivalence/CLAIMS.csv"
original <- file.path(output, "CLAIMS-before.csv")
before <- utils::read.csv(if (file.exists(original)) original else registry,
                          stringsAsFactors = FALSE)
stopifnot(nrow(before) == 67L, !anyNA(before))
if (!file.exists(original)) stopifnot(file.copy(registry, original))
replacements <- c(
  `23` = "`sandwich::vcovCL()`, applied to every robust standard error in the",
  `24` = "`sqrt(G / (G - 1))` times Mplus's. Groups are the independent units:",
  `25` = "estimator), as `glca` reports; it agrees with",
  `26` = "on them, as Mplus does; the held probabilities are omitted from the table",
  `27` = "scaling correction factor `tr(A^-1 B) / q` (with `B` without the CR1",
  `34` = "of independent units, `bic_rows` by the number of rows, as flexmix",
  `35` = "level exclude the same rows and groups, as in Mplus.",
  `37` = "counts match mclust's for all fourteen.",
  `38` = "mclust's, so a fit from the same start reproduces",
  `39` = "The mixing proportions have no prior. As in mclust, the fit's",
  `41` = "hypervolume), so the criteria equal mclust's. A random start puts a tenth",
  `42` = "classification: `subject_profiles` is `0` for an observation assigned to",
  `43` = "path as mclust's `me()` from the same start.",
  `44` = "prior. With `TRUE` they are excluded, as in Mplus, from the individual and",
  `47` = "GOLD's default ordinal model, against which it is checked. The category",
  `49` = "`flexmix(y ~ x | id)`. A group's likelihood contribution is the product",
  `50` = "estimates the two likelihoods agree",
  `52` = "covariance, as `idiag = TRUE` in lcmm); default `FALSE`.",
  `56` = "does for `mclust::Mclust()`. Pass the result as `multilpa(prior = )`. The",
  `59` = "structure-specific departures from a textbook conjugate update, so the",
  `60` = "(default `0.01`, as in mclust). Zero leaves the means unshrunk.",
  `63` = "others by EM. Agreement with Mplus (User's Guide examples 8.13, 8.14;",
  `64` = "others by EM. Agreement with Mplus (User's Guide examples 8.13, 8.14;",
  `65` = "others by EM. Agreement with Mplus (User's Guide examples 8.13, 8.14;",
  `66` = "otherwise). Results agree with Latent GOLD 6.1's `step3 modal ml` latent",
  `67` = "`\"robust\"` is the sandwich over persons and `\"observed\"` the information"
)
indices <- as.integer(names(replacements))
after <- before
after$phrase[indices] <- unname(replacements)
stopifnot(identical(before[setdiff(names(before), "phrase")],
                    after[setdiff(names(after), "phrase")]),
          sum(before$phrase != after$phrase) == 26L,
          nrow(after) == nrow(before), !anyNA(after),
          !anyDuplicated(after[c("source", "phrase", "reference")]))
present <- vapply(seq_len(nrow(after)), function(index) {
  any(grepl(after$phrase[index], readLines(after$source[index], warn = FALSE),
            fixed = TRUE))
}, logical(1))
stopifnot(all(present))
locations <- vapply(indices, function(index) {
  paste(which(grepl(after$phrase[index], readLines(after$source[index], warn = FALSE),
                    fixed = TRUE)), collapse = ";")
}, character(1))
changes <- data.frame(row = indices, source = after$source[indices],
                      line = locations, reference = after$reference[indices],
                      old_phrase = before$phrase[indices],
                      new_phrase = after$phrase[indices])
utils::write.csv(changes, file.path(output, "phrase-changes.csv"), row.names = FALSE)
utils::write.csv(after, registry, row.names = FALSE)
cat("All 67 claims retained; 26 exact source phrases reconciled.\n")
cat("Sources, references, backing tests, notes and test code unchanged.\n")
