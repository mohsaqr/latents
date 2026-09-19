# Compare Latent GOLD output against this package. Run from the project root
# once the Latent GOLD output files are in this directory:
#   Rscript validation/latentgold/compare.R
#
# Latent GOLD has not been run yet, so this script reports what is missing
# rather than pretending to a comparison it cannot make.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("validation", "latentgold")
targets <- readRDS(file.path(directory, "multilpa-targets.rds"))

output_files <- list.files(directory, pattern = "[.](lst|txt|out)$",
                           full.names = TRUE)
if (length(output_files) == 0L) {
  cat("No Latent GOLD output found in", directory, "\n\n")
  cat("Targets this package produces, awaiting a comparison:\n\n")
  cat("--- model 1, bvr.dat: two-level, 2 x 2 classes, three continuous ---\n")
  cat("log likelihood ", format(targets$bvr$log_likelihood, digits = 10),
      " on ", targets$bvr$n_parameters, " parameters\n", sep = "")
  cat("bivariate residuals, worst three:\n")
  print(head(targets$bvr$residuals[c("profile", "indicator_1", "indicator_2",
                                     "observed", "p_value")], 3),
        digits = 4, row.names = FALSE)
  cat("\n--- model 2, threestep.dat: 2 classes, then Step-3 ---\n")
  cat("log likelihood ", format(targets$three_step$log_likelihood, digits = 10),
      "\n", sep = "")
  cat("classification errors:\n")
  print(targets$three_step$errors, digits = 5, row.names = FALSE)
  cat("distal outcome under BCH:\n")
  print(targets$three_step$distal_bch[c("class", "estimate", "standard_error")],
        digits = 5, row.names = FALSE)
  cat("covariate on membership:\n")
  print(targets$three_step$covariate[c("outcome", "term", "estimate",
                                       "standard_error")],
        digits = 5, row.names = FALSE)
  cat("\nAdd the Latent GOLD output to", directory,
      "and run this again.\n")
  quit(save = "no", status = 0L)
}

cat("Found Latent GOLD output:\n"); print(basename(output_files))
cat("\nThe parser is not written yet, because the exact output layout is not\n")
cat("known until a real file exists. Write it against these files, following\n")
cat("the pattern in validation/mplus/*/compare.R, which parses genuine .out\n")
cat("artifacts the same way.\n")
