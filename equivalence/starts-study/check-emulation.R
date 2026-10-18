# Checks the premise of simulate.R: the first s starts of a longer run are the
# starts of an s-start fit with the same seed, so blocks of a 500-start run
# emulate s-start fits. Run from the project root:
#   Rscript equivalence/starts-study/check-emulation.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("equivalence", "starts-study", "simulate.R"))
row <- subset(starts_design(), dataset == "lpa_k3_moderate_n300_r1")
data <- starts_simulate(row)$data
long <- starts_fit(data, "lpa", 4L, 40L, seed = 7)$fit
ten <- starts_fit(data, "lpa", 4L, 10L, seed = 7)$fit
twenty <- starts_fit(data, "lpa", 4L, 20L, seed = 7)$fit
stopifnot(
  "10 starts match the first 10 of 40" =
    isTRUE(all.equal(head(long$starts$log_likelihood, 10L), ten$starts$log_likelihood)),
  "20 starts match the first 20 of 40" =
    isTRUE(all.equal(head(long$starts$log_likelihood, 20L), twenty$starts$log_likelihood))
)
message("emulation premise holds")
