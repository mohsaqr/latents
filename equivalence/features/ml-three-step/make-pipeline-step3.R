# The second pass of the Latent GOLD pipeline cases: step three read from
# Latent GOLD's own step-one posteriors (its outfile, six decimals).
#
# From the project root, after the step-one cases ran:
#   Rscript equivalence/features/ml-three-step/make-pipeline-step3.R
#   JOBS=2 sh equivalence/features/run.sh ml-three-step latentgold _s3

sys.source(file.path("equivalence", "features", "ml-three-step", "helpers.R"),
           envir = environment())
cases <- readRDS(file.path(ml3_kit(), "cases.rds"))
pipeline <- Filter(\(case) isTRUE(case$pipeline), cases)
invisible(lapply(pipeline, \(case) {
  step1 <- sprintf("pl_%s_s1", case$name)
  step3 <- sprintf("pl_%s_s3", case$name)
  source_file <- file.path(ml3_returned(), paste0(step1, "_post.txt"))
  stopifnot("run the step-one pipeline cases first" = file.exists(source_file))
  file.copy(source_file, file.path(ml3_kit(), paste0(step3, "-post.txt")), overwrite = TRUE)
  models <- data.frame(se = "robust", variance = c("varying", "equal"), tag = c("vr", "er"))
  writeLines(ml3_lg_step3(step3, paste0(step3, "-post.txt"), case$spec$k, models),
             file.path(ml3_kit(), paste0(step3, ".lgs")), sep = "\r\n")
}))
cat(sprintf("%d pipeline step-three files written\n", length(pipeline)))
