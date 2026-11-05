suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
ns <- asNamespace("latents"); tests <- "/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"
mutate <- function(name, from, to) { f <- get(name, ns); txt <- deparse(f, width.cutoff = 500L); hit <- grepl(from, txt, fixed = TRUE)
  stopifnot("mutation target not found" = any(hit)); txt <- gsub(from, to, txt, fixed = TRUE); g <- eval(parse(text = txt)); environment(g) <- ns
  unlockBinding(name, ns); assign(name, g, ns); lockBinding(name, ns); f }
restore <- function(name, f) { unlockBinding(name, ns); assign(name, f, ns); lockBinding(name, ns) }
run <- function(files) { res <- vapply(files, function(file) { r <- as.data.frame(testthat::test_file(file.path(tests, file), reporter = "silent", package = "latents", load_package = "none", stop_on_failure = FALSE)); c(sum(r$failed), sum(r$error), sum(r$nb)) }, numeric(3)); rownames(res) <- c("failed", "error", "n"); res }
trial <- function(label, name, from, to, files) { old <- mutate(name, from, to); on.exit(restore(name, old)); t0 <- Sys.time(); out <- run(files); cat("\n==", label, "(", round(difftime(Sys.time(), t0, units = "secs")), "s)\n"); print(out) }
cl_files <- c("test-cross-level.R", "test-review-1015-families.R", "test-review-followup.R", "test-review-families.R")
trial("cross-level count: full counts prevalences once", ".cross_level_count", "n_group_classes * (n_profiles - 1L)", "(n_profiles - 1L)", cl_files)
trial("restricted composition replaced by a pooled share", ".cross_level_result", "crossprod(weight, subject_posteriors)/colSums(weight)", "matrix(colSums(subject_posteriors)/nrow(subject_posteriors), n_group_classes, n_profiles, byrow = TRUE)", cl_files)
trial("full prevalence: unweighted class average", ".cross_level_result", "drop(crossprod(p$group_probabilities, p$profile_probabilities))", "colMeans(p$profile_probabilities)", cl_files)
ad_files <- c("test-additive-compare.R", "test-audit-1101-additive-compare.R", "test-review-1026-additive.R")
trial("family bootstrap simulate: sd = variance (missing sqrt)", ".additive_simulate", "sqrt(x$between_variances[classes, , drop = FALSE])", "x$between_variances[classes, , drop = FALSE]", ad_files)
trial("additive log-interval doubled", ".multilpa_wald_bounds", "spread <- critical * standard_error[positive]/estimate[positive]", "spread <- 2 * critical * standard_error[positive]/estimate[positive]", c("test-additive-methods.R", "test-review-1015-families.R"))
trial("logit interval replaced by Wald spread", ".multilpa_wald_bounds", "low[probability] <- stats::plogis(centre - spread)", "low[probability] <- estimate[probability] - critical * standard_error[probability]", c("test-additive-methods.R", "test-review-1015-families.R"))
