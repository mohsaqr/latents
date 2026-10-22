# Every comparison of the direct-effect kit, written as tidy tables.
#
# Run from the project root after make-kit.R and run.sh:
#   Rscript equivalence/features/direct-effects/compare.R
# Writes comparison.csv (one row per compared quantity, case and program),
# maxima.csv (each program's log likelihood against the best) and
# reference.csv (the independent-implementation battery) next to this file.
# equivalence/test-features-direct-effects.R asserts the same rows.
pkgload::load_all(".", quiet = TRUE)
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
sys.source(file.path(root, "reference.R"), envir = globalenv())
set.seed(20261020)

external <- de_compare_all(root)
utils::write.csv(external$comparison, file.path(root, "comparison.csv"), row.names = FALSE)
utils::write.csv(external$maxima, file.path(root, "maxima.csv"), row.names = FALSE)
if (length(external$failed)) message("references that did not terminate: ",
                                     paste(external$failed, collapse = ", "))
design <- ref_design()
battery <- do.call(rbind, lapply(seq_len(nrow(design$grid)), \(i) {
  ref_check(design$grid[i, ], design$effect_sets)
}))
utils::write.csv(battery, file.path(root, "reference.csv"), row.names = FALSE)

worst <- aggregate(cbind(difference, ratio = difference / pmax(tolerance, 1e-300)) ~
                     quantity + program, data = external$comparison, FUN = max)
print(worst)
print(subset(external$comparison, difference > tolerance))
print(subset(external$maxima, abs(gap_to_best) > 1e-6 * abs(best_reference)))
