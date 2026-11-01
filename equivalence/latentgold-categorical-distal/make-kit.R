# Categorical distal outcomes in three_step() against Latent GOLD 6.1 Step3.
# Run from the project root:
#   Rscript equivalence/latentgold-categorical-distal/make-kit.R
#   sh equivalence/latentgold-categorical-distal/run-kit.sh
#   Rscript equivalence/latentgold-categorical-distal/compare.R
#
# Each case (common.R) is simulated and fitted with latents; each unit's
# step-one posteriors are written with the distal outcome, and Latent GOLD
# reads them (`posterior = (...)`), assigns modally, builds the
# classification errors and fits the step-three model of the outcome as a
# nominal dependent variable, once per method, with model-based and with
# robust standard errors. returned/INPUTS.md5 fingerprints the kit Latent
# GOLD was given.
pkgload::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-categorical-distal")
source(file.path(base, "common.R"))
kit_dir <- file.path(base, "kit")
dir.create(kit_dir, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(distal_cases, function(case) {
  made <- distal_case(case)
  utils::write.table(distal_data_file(made), file.path(kit_dir, sprintf("%s.dat", case)),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  lapply(names(distal_methods), function(method) {
    lapply(c(FALSE, TRUE), function(robust) {
      file <- sprintf("%s_%s%s.lgs", case, method, if (robust) "_robust" else "")
      writeLines(distal_syntax(case, method, robust, made$fit$n_profiles),
                 file.path(kit_dir, file))
    })
  })
}))
files <- sort(list.files(kit_dir, pattern = "[.](dat|lgs)$", full.names = TRUE), method = "radix")
writeLines(sprintf("%s  %s", unname(tools::md5sum(files)), basename(files)),
           file.path(kit_dir, "MANIFEST.md5"))
cat("kit written:", length(files), "files\n")
