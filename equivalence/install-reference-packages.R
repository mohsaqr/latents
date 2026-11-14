# Install the reference implementations the equivalence suite compares against.
# They are not in DESCRIPTION (the folder never ships), so CI installs them here.
# From the project root: Rscript equivalence/install-reference-packages.R
reference_packages <- c(
  "devtools", "pkgload", "testthat", "withr", "numDeriv", "mclust", "MASS", "sandwich",
  "knitr", "rmarkdown", "glca", "poLCA", "multilevLCA", "depmixS4", "flexmix", "lcmm",
  "nnet", "nlme", "lme4", "mvtnorm", "LMest", "VGAM", "mice", "CompQuadForm", "tidySEM",
  "ordinal", "randomLCA", "BayesLCA", "Rmixmod", "EMCluster", "mixtools", "MoEClust",
  "seqHMM", "TraMineR", "mitml", "mitools", "tidyLPA", "tna", "dplyr", "ggplot2",
  "nnls", "Rmpfr", "jsonlite")
missing_packages <- setdiff(reference_packages, rownames(utils::installed.packages()))
if (length(missing_packages)) utils::install.packages(missing_packages)
# Pinned reference versions, installed from the CRAN archive when another
# version is present. lcmm 2.2.2 returns no subject-specific predictions for
# the growth mixtures of test-growth-predict-lcmm.R (its blocks then skip and
# the lcmm reference is incomplete); 2.2.1 does, and is what the 2026-11-01
# matrix was built with.
pinned <- c(lcmm = "2.2.1")
outdated <- names(pinned)[vapply(names(pinned), function(p) {
  !p %in% rownames(utils::installed.packages()) ||
    as.character(utils::packageVersion(p)) != pinned[[p]]
}, logical(1))]
invisible(lapply(outdated, function(p) {
  utils::install.packages(
    sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", p, p, pinned[[p]]),
    repos = NULL, type = "source")
  if (as.character(utils::packageVersion(p)) != pinned[[p]]) {
    stop(sprintf("%s %s could not be installed", p, pinned[[p]]), call. = FALSE)
  }
}))
still_missing <- setdiff(reference_packages, rownames(utils::installed.packages()))
if (length(still_missing)) {
  stop("Not installed: ", paste(still_missing, collapse = ", "), call. = FALSE)
}
versions <- vapply(reference_packages, \(p) as.character(utils::packageVersion(p)), character(1))
print(data.frame(package = reference_packages, version = versions, row.names = NULL))
