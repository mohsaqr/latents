# Finding 2: a class identified by a single unit gets a standard error of
# exactly 0 (and a zero-width interval) under variance = "equal", silently.
#
# From the project root:
#   Rscript equivalence/features/ml-three-step/repro-zero-se.R
#
# sweep_041 of the kit: four profiles on 50 units; one unit is alone in
# profile 3 with posterior 1 there and 0 elsewhere. With a shared variance
# the likelihood is bounded and the class mean is that unit's outcome; its
# score for the mean is then zero at the estimate, so the sandwich has a zero
# row and column for it. The model-based (inverse information) error is about
# 66. With variance = "varying" the same class's variance collapses and
# three_step() refuses (latents_singular_information) after a boundary
# warning, which is the behaviour one would want here too.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
sys.source(file.path("equivalence", "features", "ml-three-step", "helpers.R"),
           envir = environment())
sys.source(file.path("equivalence", "features", "ml-three-step", "reference.R"),
           envir = environment())
case <- readRDS(file.path(ml3_kit(), "cases.rds"))$sweep_041
fit <- suppressWarnings(ml3_fit(case))
print(tabulate(fit$subject_profiles, 4L))
print(three_step(fit, case$data, "y", method = "ml", variance = "equal",
                 vcov_type = "independent"))
print(three_step(fit, case$data, "y", method = "ml", variance = "equal",
                 vcov_type = "independent", contrast = "pairs"))
reference <- ref_step3(unname(fit$subject_posteriors), case$data$y, "equal")
print(reference$classes)
print(tryCatch(three_step(fit, case$data, "y", method = "ml", variance = "varying"),
               error = \(e) class(e)))
