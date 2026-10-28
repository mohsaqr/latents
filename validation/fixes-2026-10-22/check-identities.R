Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1")
pkgload::load_all(quiet = TRUE)
before <- readRDS("/private/tmp/lta-thin-before.rds")
after <- readRDS("/private/tmp/lta-thin-after.rds")
same_lta <- unlist(Map(function(left, right) {
  Map(function(a, b) identical(a$fit, b$fit), left, right)
}, before, after))
stopifnot(length(same_lta) == 4L, all(same_lta))
cat("LTA returned fits identical before/after:", paste(same_lta, collapse = ", "), "\n")

# Compare the previous bootstrap implementation against the corrected path
# on an unaffected two-start fit, including the caller's RNG state.
old_code <- system2("git", c("show", "HEAD:R/bootstrap-inference.R"), stdout = TRUE)
previous <- new.env(parent = asNamespace("latents"))
eval(parse(text = old_code), envir = previous)
quiet <- function(expr) withCallingHandlers(expr, warning = function(condition) {
  if (inherits(condition, c("latents_unconverged", "latents_boundary",
                            "latents_small_classes", "latents_thin_profiles",
                            "latents_failed_starts"))) invokeRestart("muffleWarning")
})
set.seed(11)
data <- data.frame(g = rep(seq_len(15L), each = 8L),
                   a = rnorm(120, rep(c(-2, 2), each = 60L)),
                   b = rnorm(120, rep(c(-2, 2), each = 60L)))
fit <- quiet(fit_staged(data, c("a", "b"), "g", 2L, 2L,
                         n_starts = 2L, n_screen = 8L, search = "emEM",
                         max_iter = 300L, seed = 1L))
set.seed(2)
old <- quiet(previous$.multilpa_bootstrap_replicates(fit, data, 2L, 2L, 300L, 1e-8))
old_rng <- .Random.seed
set.seed(2)
new <- quiet(.multilpa_bootstrap_replicates(fit, data, 2L, 2L, 300L, 1e-8))
new_rng <- .Random.seed
stopifnot(all(complete.cases(old$estimates)), identical(old, new),
          identical(old_rng, new_rng))
cat("Unaffected staged-bootstrap estimates/messages identical:", identical(old, new), "\n")
cat("Unaffected staged-bootstrap RNG state identical:", identical(old_rng, new_rng), "\n")
saveRDS(list(lta_identical = same_lta, old_bootstrap = old, new_bootstrap = new,
             rng_identical = identical(old_rng, new_rng)),
        "/private/tmp/latents-fixes-20261022/identities.rds")
