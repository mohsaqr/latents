# The numbers of README.md: every comparison of the kit, and the largest
# difference per quantity against each reference. Run from the project root:
#   Rscript equivalence/features/vlmr/summarise.R
# VLMR_FIT_CACHE=<dir> reuses fitted grids kept by an earlier run.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
here <- file.path("equivalence", "features", "vlmr")
old <- setwd("equivalence")
on.exit(setwd(old), add = TRUE)
sys.source(file.path("features", "vlmr", "helpers.R"), envir = environment())
sys.source(file.path("features", "vlmr", "reference.R"), envir = environment())
run <- vlmr_compute_all(names(vlmr_cases()), cores = 2L)
setwd(old)
table <- run$table
utils::write.csv(table, file.path(here, "comparisons.csv"), row.names = FALSE)

largest <- function(values) if (all(is.na(values))) NA_real_ else max(abs(values), na.rm = TRUE)
matched <- subset(table, !is.na(lg_lr) & abs(ll0 - lg_ll0) < 1e-4 & abs(ll1 - lg_ll1) < 1e-4)
mplus <- subset(table, !is.na(mp_lr) & !is.na(lr))
reference <- subset(table, !is.na(ref_lr))
summary <- rbind(
  data.frame(reference = "Latent GOLD", pairs = nrow(matched),
             quantity = c("log likelihood (null, alternative)", "LR", "adjusted LR", "VLMR p",
                          "aLMR p", "reference mean", "reference sd"),
             max_difference = c(largest(c(matched$ll0 - matched$lg_ll0,
                                          matched$ll1 - matched$lg_ll1)),
                                largest(matched$lr - matched$lg_lr),
                                largest(matched$adjusted_lr - matched$lg_adjusted_lr),
                                largest(matched$vlmr_p - matched$lg_vlmr_p),
                                largest(matched$almr_p - matched$lg_almr_p),
                                largest(matched$reference_mean - matched$lg_mean),
                                largest(matched$reference_sd - matched$lg_sd))),
  data.frame(reference = "Mplus TECH11", pairs = nrow(mplus),
             quantity = c("H0 log likelihood", "H1 log likelihood", "LR", "adjusted LR",
                          "VLMR p", "aLMR p", "reference mean", "reference sd"),
             max_difference = c(largest(mplus$mp_h0_ll - mplus$ll0),
                                largest(mplus$mp_h1_ll - mplus$ll1),
                                largest(mplus$lr - mplus$mp_lr),
                                largest(mplus$adjusted_lr - mplus$mp_adjusted_lr),
                                largest(mplus$vlmr_p - mplus$mp_vlmr_p),
                                largest(mplus$almr_p - mplus$mp_almr_p),
                                largest(mplus$reference_mean - mplus$mp_mean),
                                largest(mplus$reference_sd - mplus$mp_sd))),
  data.frame(reference = "independent R reference", pairs = nrow(reference),
             quantity = c("LR (relative)", "weights (largest eigenvalue gap)", "reference mean",
                          "reference sd", "VLMR p vs CompQuadForm::imhof",
                          "VLMR p vs CompQuadForm::davies", "aLMR p vs davies",
                          "VLMR p vs Monte Carlo (in binomial SEs)",
                          "summed score at the estimate, per unit"),
             max_difference = c(largest((reference$lr - reference$ref_lr) / reference$lr),
                                largest(reference$lambda_gap),
                                largest(reference$reference_mean - reference$ref_mean),
                                largest(reference$reference_sd - reference$ref_sd),
                                largest(reference$vlmr_p - reference$ref_imhof_p),
                                largest(reference$vlmr_p - reference$ref_davies_p),
                                largest(reference$almr_p - reference$ref_davies_almr_p),
                                largest((reference$vlmr_p - reference$ref_mc_p) /
                                          sqrt(pmax(reference$vlmr_p * (1 - reference$vlmr_p),
                                                    1e-6) / reference$ref_draws)),
                                largest(reference$score_sum / reference$n_units))))
stress <- vlmr_imhof_stress()
stress_summary <- do.call(rbind, lapply(split(stress, stress$family), function(part) {
  data.frame(family = part$family[1L], cases = nrow(part), reference = part$reference_kind[1L],
             max_abs_error = max(part$abs_error),
             max_rel_error = largest(part$rel_error))
}))
utils::write.csv(summary, file.path(here, "summary.csv"), row.names = FALSE)
utils::write.csv(stress_summary, file.path(here, "imhof-stress.csv"), row.names = FALSE)
print(summary, digits = 3)
print(stress_summary, digits = 3)
print(run$errors)
