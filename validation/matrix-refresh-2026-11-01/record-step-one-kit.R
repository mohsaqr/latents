# Rerun the step-one correction comparison with Latent GOLD's saved listings
# and bind comparison.csv to the runtime source it was produced with, so the
# verification matrix can count it as current. Run from the project root.
source(file.path("equivalence", "verification-matrix", "provenance.R"))
snapshot <- matrix_capture_source(".")
source(file.path("equivalence", "latentgold-step-one", "compare.R"))
matrix_record_evidence(file.path("equivalence", "latentgold-step-one", "comparison.csv"),
                       snapshot, ".")
