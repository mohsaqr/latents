# Calibration of the two transition features added after 0.9.8:
#  (a) percentile-interval coverage of parameter_inference(method = "bootstrap")
#      on a full-covariance (VVV) occasion-transition model, where Wald
#      inference is unavailable;
#  (b) size of bootstrap_lrt() with a full-covariance null (first vs second
#      order), which needs full-covariance simulation.
# Truth is a fitted model; datasets are simulated from it on a fixed design.
# Run from the repository root:
#   Rscript validation/review-2026-10-02-followup/transition-bootstrap-validation.R
devtools::load_all(".", quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(20261002)
out_dir <- file.path("validation", "review-2026-10-02-followup")
n_datasets <- as.integer(Sys.getenv("TB_N", "200"))
# Sensitivity runs: TB_ITER replicates per interval; TB_SIZE = 0 skips (b).
boot_iter <- as.integer(Sys.getenv("TB_ITER", "99"))
run_size <- !identical(Sys.getenv("TB_SIZE", "1"), "0")
suffix <- if (identical(boot_iter, 99L)) "" else sprintf("-iter%d", boot_iter)
cores <- max(1L, parallel::detectCores() - 2L)

# Design: 240 persons, 4 occasions, two correlated indicators.
template <- function(n, occasions = 4L) {
  data.frame(id = rep(seq_len(n), each = occasions), time = rep(seq_len(occasions), n),
             y1 = 0, y2 = 0)
}
design <- template(240L)
seed_data <- local({
  set.seed(41)
  states <- matrix(NA_integer_, 240L, 4L)
  states[, 1L] <- sample.int(2L, 240L, replace = TRUE)
  invisible(lapply(2:4, function(t) {
    states[, t] <<- vapply(states[, t - 1L], function(from)
      sample.int(2L, 1L, prob = if (from == 1L) c(0.85, 0.15) else c(0.15, 0.85)),
      integer(1))
  }))
  s <- as.vector(t(states))
  first <- c(0, 3)[s] + rnorm(length(s))
  data.frame(id = design$id, time = design$time, y1 = first,
             y2 = c(0.5, 3.5)[s] + 0.6 * first + rnorm(length(s)))
})
truth_fit <- lta(seed_data, c("y1", "y2"), "id", 2, time = "time",
                 transitions = "occasion", model = "VVV", n_starts = 5, seed = 1)
truth <- .lta_bootstrap_estimates(truth_fit)

coverage_one <- function(i) {
  data <- .lta_simulate(truth_fit, design)
  fit <- tryCatch(suppressWarnings(lta(data, c("y1", "y2"), "id", 2, time = "time",
                                       transitions = "occasion", model = "VVV",
                                       n_starts = 3, seed = i)),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  boot <- tryCatch(suppressWarnings(parameter_inference(
    fit, data = data, method = "bootstrap", iter = boot_iter, n_starts = 2, seed = i)),
    error = function(e) NULL)
  if (is.null(boot)) return(NULL)
  # The fit's labels are its own (truth profile j is fit profile order[j]).
  # Each replicate is relabelled into the truth's labels, rebasing included,
  # and the percentile interval is read there.
  order <- .multilpa_match_order(.lta_profile_signature(truth_fit, truth_fit),
                                 .lta_profile_signature(fit, truth_fit))
  replicates <- attr(boot, "replicates")
  in_truth <- t(apply(replicates, 1L, function(row) {
    .lta_relabel_profiles(stats::setNames(row, colnames(replicates)), order)[names(truth)]
  }))
  bounds <- apply(in_truth, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  estimate <- .lta_align_estimates(fit, truth_fit)[names(truth)]
  data.frame(dataset = i, parameter = names(truth), truth = unname(truth),
             estimate = unname(estimate),
             covered = unname(truth >= bounds[1L, ] & truth <= bounds[2L, ]),
             relabelled = !identical(order, seq_along(order)))
}

size_one <- function(i) {
  data <- .lta_simulate(truth_fit, design)
  fits <- lapply(1:2, function(order) tryCatch(suppressWarnings(lta(
    data, c("y1", "y2"), "id", 2, time = "time", transitions = "occasion",
    order = order, model = "VVV", n_starts = 3, seed = i)), error = function(e) NULL))
  if (any(vapply(fits, is.null, logical(1)))) return(data.frame(dataset = i, p_value = NA))
  test <- tryCatch(suppressWarnings(bootstrap_lrt(fits[[1]], fits[[2]], data = data,
                                                  iter = 39, n_starts = 1, seed = i)),
                   error = function(e) NULL)
  data.frame(dataset = i, p_value = if (is.null(test)) NA_real_ else
    get_results(test, "test")$p_value)
}

started <- Sys.time()
coverage <- do.call(rbind, parallel::mclapply(seq_len(n_datasets), coverage_one,
                                              mc.cores = cores, mc.set.seed = TRUE))
size <- if (run_size) do.call(rbind, parallel::mclapply(seq_len(n_datasets), size_one,
                                          mc.cores = cores, mc.set.seed = TRUE)) else NULL
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
utils::write.csv(coverage, file.path(out_dir, sprintf(
  "transition-bootstrap-coverage%s.csv", suffix)), row.names = FALSE)
if (run_size) utils::write.csv(size, file.path(out_dir, "transition-bootstrap-lrt-size.csv"),
                               row.names = FALSE)

usable <- subset(coverage, !is.na(covered))
by_parameter <- aggregate(covered ~ parameter, data = usable, FUN = mean)
by_parameter$mcse <- sqrt(by_parameter$covered * (1 - by_parameter$covered) /
                            as.vector(table(usable$parameter)[by_parameter$parameter]))
p <- if (run_size) size$p_value[is.finite(size$p_value)] else numeric()
lines <- c(
  sprintf("Transition bootstrap validation (%s; %.1f min; %d cores)", format(Sys.Date()),
          elapsed, cores),
  sprintf("Truth: VVV occasion-transition lta on 240 persons x 4 occasions; %d datasets.",
          n_datasets),
  "",
  sprintf("(a) 95%% percentile-interval coverage, iter = %d per dataset:", boot_iter),
  sprintf("  usable datasets %d of %d (of which label-switched and matched: %d)",
          length(unique(usable$dataset)), n_datasets,
          length(unique(subset(usable, relabelled)$dataset))),
  utils::capture.output(print(transform(by_parameter, covered = round(covered, 3),
                                        mcse = round(mcse, 3)), row.names = FALSE)),
  sprintf("  pooled coverage %.3f", mean(usable$covered)),
  sprintf("  by label switching: matched-as-is %.3f, relabelled %.3f",
          mean(subset(usable, !relabelled)$covered), mean(subset(usable, relabelled)$covered)),
  "",
  if (run_size) c(
  "(b) bootstrap_lrt size, VVV first vs second order, null true, iter = 39:",
  sprintf("  usable %d of %d; rejection rate at 0.05: %.3f (MCSE %.3f); at 0.10: %.3f",
          length(p), n_datasets, mean(p <= 0.05),
          sqrt(mean(p <= 0.05) * (1 - mean(p <= 0.05)) / length(p)), mean(p <= 0.10)),
  sprintf("  KS test of p-value uniformity: p = %.3f",
          suppressWarnings(stats::ks.test(p, "punif")$p.value))))
writeLines(lines, file.path(out_dir, sprintf("transition-bootstrap-validation%s.txt", suffix)))
writeLines(lines)
