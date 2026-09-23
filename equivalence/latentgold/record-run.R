# Record a Latent GOLD run, so the output is retained as evidence and Latent
# GOLD is not needed again.
#
# Run from the project root after copying Latent GOLD's output into
# equivalence/latentgold/returned/:
#   Rscript equivalence/latentgold/record-run.R
#
# Writes into returned/:
#   INPUTS.md5    the md5 of every .dat and .lgs Latent GOLD was given, which
#                 compare.R checks, so a rebuilt kit cannot be compared against
#                 output produced from a different one;
#   OUTPUTS.md5   the md5 of every retained listing and posterior file;
#   PROVENANCE.md what produced them, on what, and when.
#
# The environment is not read from the machine: pass it, so a record can be
# written for output produced elsewhere.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "latentgold")
source(file.path(directory, "compare-functions.R"))
kit <- file.path(directory, "kit")
returned <- file.path(directory, "returned")

program <- "Latent GOLD 6.1 (Basic+Advanced/Syntax+Choice), free academic single-user licence"
platform <- paste("Wine Devel 11.17 (Gcenx build) on macOS 26.3.1, arm64;",
                  "prefix ~/.wine-latentgold")
command <- "lg61.exe <case>.lgs /b /o <case>.lst"

outputs <- sort(list.files(returned, pattern = "[.](lst|txt)$"))
if (length(outputs) == 0L) {
  stop(errorCondition("No Latent GOLD output in `returned/` to record.",
                      class = "multilpa_missing_latentgold_output", call = NULL))
}

.write_md5 <- function(fingerprint, path) {
  writeLines(sprintf("%s  %s", unname(fingerprint), names(fingerprint)), path)
}

inputs <- lg_fingerprint(kit)
.write_md5(inputs, file.path(returned, "INPUTS.md5"))
.write_md5(stats::setNames(unname(tools::md5sum(file.path(returned, outputs))), outputs),
           file.path(returned, "OUTPUTS.md5"))

writeLines(c(
  "# Retained Latent GOLD output",
  "",
  sprintf("Produced %s by %s,", format(Sys.Date()), program),
  sprintf("running under %s.", platform),
  "",
  sprintf("Each model was estimated in batch as `%s`, from a directory inside the", command),
  "Wine drive. The kit that produced it was built by `make-kit.R` under multilpa",
  sprintf("%s.", utils::packageVersion("latents")),
  "",
  "## Why it is kept",
  "",
  "Latent GOLD is commercial and Windows-only, and its licence is per user. The",
  "listings and posterior files here are the evidence, so `compare.R` runs",
  "offline and nobody needs Latent GOLD again unless the kit's data or syntax",
  "change. `INPUTS.md5` fingerprints exactly those inputs; `compare.R` warns with",
  "`multilpa_stale_latentgold_output` when the kit no longer matches.",
  "",
  "## Files",
  "",
  sprintf("- %d listings and posterior files (`OUTPUTS.md5`)", length(outputs)),
  sprintf("- from %d input files (`INPUTS.md5`)", length(inputs)),
  "",
  "These are Latent GOLD's own output, kept for verification only. They are not",
  "redistributed as part of the package: `equivalence/` is in `.Rbuildignore`."),
  file.path(returned, "PROVENANCE.md"))

message(sprintf("Recorded %d outputs from %d inputs in %s.", length(outputs),
                length(inputs), returned))
