# Latent GOLD's step three starts from non-random values (the proportional
# moments), so on a multimodal step-three likelihood it can stop at a lower
# maximum than latents. These inputs refit the cases where that can happen
# with `simultaneous` estimation and 100 random start sets, so the comparison
# is against Latent GOLD's best maximum, not its first.
#
# From the project root, after make-kit.R:
#   Rscript equivalence/features/ml-three-step/make-random-starts.R
#   JOBS=2 sh equivalence/features/run.sh ml-three-step latentgold x_

sys.source(file.path("equivalence", "features", "ml-three-step", "helpers.R"),
           envir = environment())
cases <- readRDS(file.path(ml3_kit(), "cases.rds"))
# Sweep cases on which the reference found a higher maximum than latents
# (finding 1, repro-local-maximum.R): their default
# Latent GOLD run is written here too, so both of its maxima are on record.
local_sweep <- sprintf("sweep_%03d", c(18L, 28L, 30L, 83L, 105L, 108L))
invisible(lapply(cases[local_sweep], \(case) {
  ml3_write_lg_data(file.path(ml3_kit(), paste0(case$name, "-data.txt")), case$data$id,
                    case$data$id, case$data$y, case$posteriors)
  writeLines(ml3_lg_step3(case$name, paste0(case$name, "-data.txt"), case$spec$k,
                          data.frame(se = "robust", variance = c("varying", "equal"),
                                     tag = c("vr", "er"))),
             file.path(ml3_kit(), paste0(case$name, ".lgs")), sep = "\r\n")
}))
multimodal <- c(Filter(\(case) isTRUE(case$lg) && (startsWith(case$name, "e_") ||
                                                     grepl("_lo_", case$name)), cases),
                cases[local_sweep])
invisible(lapply(multimodal, \(case) {
  name <- paste0("x_", case$name)
  file.copy(file.path(ml3_kit(), paste0(case$name, "-data.txt")),
            file.path(ml3_kit(), paste0(name, "-data.txt")), overwrite = TRUE)
  models <- data.frame(se = "robust", variance = c("varying", "equal"), tag = c("vr", "er"))
  lines <- ml3_lg_step3(name, paste0(name, "-data.txt"), case$spec$k, models)
  lines <- sub("step3 modal ml;", "step3 modal ml simultaneous;", lines, fixed = TRUE)
  at <- grep("^   bayes$", lines)
  lines <- append(lines, c("   startvalues",
                           "      seed=20261020 sets=100 tolerance=1e-005 iterations=250;"),
                  after = at[1L] - 1L)
  at <- grep("^   bayes$", lines)
  if (length(at) > 1L) {
    lines <- append(lines, c("   startvalues",
                             "      seed=20261020 sets=100 tolerance=1e-005 iterations=250;"),
                    after = at[2L] - 1L)
  }
  writeLines(lines, file.path(ml3_kit(), paste0(name, ".lgs")), sep = "\r\n")
}))
cat(sprintf("%d random-start inputs written\n", length(multimodal)))
