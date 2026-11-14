# Recovery of known parameters: every public function fitted to data that
# simulab drew from known parameters, the estimates compared with those
# parameters over many replicates (bias within 3.5 Monte Carlo standard
# errors, 95% interval coverage within its 99.9% binomial range), and the
# selection and test functions checked for choosing the generating number of
# classes (simulab-recovery/). Runs as its own ci job, equivalence-r.

skip_if_not_installed("simulab")
skip_if_not_installed("mice")

recovery_dir <- file.path("simulab-recovery")
source(file.path(recovery_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

kit_reference_tests(
  "simulab-recovery", compare_all_recovery, n_rows = 279L, groups = "simulation recovery")
