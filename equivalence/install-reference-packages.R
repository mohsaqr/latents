# Install the reference implementations the equivalence suite compares against.
# They are not in DESCRIPTION (the folder never ships), so CI installs them here.
# From the project root: Rscript equivalence/install-reference-packages.R
reference_packages <- c(
  "devtools", "pkgload", "testthat", "withr", "numDeriv", "mclust", "MASS", "sandwich",
  "knitr", "rmarkdown", "glca", "poLCA", "multilevLCA", "depmixS4", "flexmix", "lcmm",
  "nnet", "nlme", "lme4", "mvtnorm", "LMest", "VGAM", "mice", "CompQuadForm", "tidySEM",
  "ordinal", "randomLCA", "BayesLCA", "Rmixmod", "EMCluster", "mixtools", "MoEClust",
  "seqHMM", "TraMineR", "mitml", "mitools", "tidyLPA", "tna", "dplyr", "ggplot2")
missing_packages <- setdiff(reference_packages, rownames(utils::installed.packages()))
if (length(missing_packages)) utils::install.packages(missing_packages)
still_missing <- setdiff(reference_packages, rownames(utils::installed.packages()))
if (length(still_missing)) {
  stop("Not installed: ", paste(still_missing, collapse = ", "), call. = FALSE)
}
versions <- vapply(reference_packages, \(p) as.character(utils::packageVersion(p)), character(1))
print(data.frame(package = reference_packages, version = versions, row.names = NULL))
