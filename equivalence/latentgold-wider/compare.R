# Compares latents with the retained Latent GOLD 6.1 listings in returned/.
# Run from the project root: Rscript equivalence/latentgold-wider/compare.R
# Writes comparison.csv beside this script and exits non-zero on any
# disagreement. Latent GOLD is not needed: its output is kept in returned/.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "latentgold-wider")
source(file.path(directory, "compare-functions.R"))
cases <- readRDS(file.path(directory, "kit", "cases.rds"))
returned <- file.path(directory, "returned")

comparison <- compare_all(cases, returned)

utils::write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
summary_table <- do.call(rbind, lapply(split(comparison, comparison$case), function(rows)
  data.frame(case = rows$case[1L], compared = nrow(rows), agreed = sum(rows$agrees),
             worst_difference = max(rows$difference[rows$tolerance > 0]))))
print(summary_table, row.names = FALSE)
cat(sprintf("\n%d quantities compared, %d agree, %d disagree\n", nrow(comparison),
            sum(comparison$agrees), sum(!comparison$agrees)))
if (!all(comparison$agrees)) {
  print(subset(comparison, !agrees), row.names = FALSE)
  quit(status = 1L, save = "no")
}
