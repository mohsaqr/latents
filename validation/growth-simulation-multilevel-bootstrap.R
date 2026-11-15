# Predeclared simulation of the multilevel growth mixture's cluster bootstrap
# (roadmap item 5, "multilevel growth bootstrap"): coverage of the percentile
# intervals of parameter_inference(method = "bootstrap"), their agreement
# with the observed-information (Hessian) Wald errors, and the consistency of
# the class and group-class matching across replicates. It adds its own
# registry rows (growth-simulation-multilevel-bootstrap-registry.csv) beside
# those of growth-simulation.R, which it does not change. The registry is
# written BEFORE any replicate runs and the script refuses to run if it
# changes.
#
# Design: clusters (schools) in two group classes of equal probability; each
# person's trajectory class depends on the school's group class
# (P(improving) 0.85 or 0.2). Improving: 50 + 3 t; flat: 55. Random
# intercept and slope, SD 2 and 0.4, uncorrelated, shared by both classes;
# residual SD 1.5 in each class. Fitted as generated: two classes, two group
# classes, random = ~ 1 + wave, random_covariance = "equal".
#
# Criteria, fixed before the first replicate (verification article, Layer 4):
#   availability    >= 98% of datasets give a converged, non-degenerate fit
#                   with Wald inference and at least 90% usable resamples;
#   bias            |mean estimate - truth| <= 3.5 Monte Carlo SE;
#   coverage        share of 95% intervals containing the truth within the
#                   central 99.9% binomial range for the available datasets,
#                   for the bootstrap percentile and the Wald intervals;
#   agreement       median over datasets of bootstrap SE / Wald SE in
#                   [0.8, 1.25] for every parameter;
#   label matching  in >= 99% of usable resamples the matched replicate keeps
#                   the original's order of the classes (by slope) and of the
#                   group classes (by the improving class's logit).
#
# Run from the project root:
#   Rscript validation/growth-simulation-multilevel-bootstrap.R
# Environment: GROWTH_ML_BOOT_N (datasets, default the registry's 200),
# GROWTH_ML_BOOT_CORES (PSOCK workers, default 8), GROWTH_ML_BOOT_OUT
# (output folder, default validation), GROWTH_ML_BOOT_ITER (resamples; fewer
# than declared only for a pilot written outside validation/). A run with fewer datasets than the
# registry declares says so in its outputs.
# Outputs: growth-simulation-multilevel-bootstrap-{replicates,summary}.csv
# and growth-simulation-multilevel-bootstrap.txt.

pkgload::load_all(".", quiet = TRUE)
out_dir <- Sys.getenv("GROWTH_ML_BOOT_OUT", "validation")
n_cores <- as.integer(Sys.getenv("GROWTH_ML_BOOT_CORES", "8"))

registry <- data.frame(
  cell = "multilevel_separated", clusters = 50L, persons_per_cluster = 8L,
  times = "0,1,2,3", group_share = "0.5,0.5", improving_by_group_class = "0.85,0.2",
  intercept = "50,55", slope = "3,0", random_sd = "2,0.4", random_correlation = 0,
  sigma = "1.5,1.5", n_starts = 3L, bootstrap_iter = 99L, bootstrap_n_starts = 1L,
  n_datasets = 200L, base_seed = 26111500L, stringsAsFactors = FALSE)
registry_file <- file.path(out_dir, "growth-simulation-multilevel-bootstrap-registry.csv")
if (file.exists(registry_file)) {
  stored <- utils::read.csv(registry_file, stringsAsFactors = FALSE)
  current <- utils::read.csv(text = paste(utils::capture.output(
    utils::write.csv(registry, row.names = FALSE)), collapse = "\n"),
    stringsAsFactors = FALSE)
  stopifnot("registry changed after it was declared" = identical(stored, current))
} else {
  utils::write.csv(registry, registry_file, row.names = FALSE)
}
n_datasets <- as.integer(Sys.getenv("GROWTH_ML_BOOT_N", registry$n_datasets))
# A pilot may use fewer resamples, but only outside the retained results.
iter <- as.integer(Sys.getenv("GROWTH_ML_BOOT_ITER", registry$bootstrap_iter))
stopifnot("fewer resamples than declared are for pilots outside validation/" =
            iter == registry$bootstrap_iter ||
            normalizePath(out_dir) != normalizePath("validation"))

# The generating values, read back from the registry.
design <- local({
  numbers <- function(text) as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]])
  list(J = registry$clusters, n = registry$persons_per_cluster,
       times = numbers(registry$times), group_share = numbers(registry$group_share),
       improving = numbers(registry$improving_by_group_class),
       intercept = numbers(registry$intercept), slope = numbers(registry$slope),
       random_sd = numbers(registry$random_sd), sigma = numbers(registry$sigma))
})

draw <- function(design, seed) {
  set.seed(seed)
  group_class <- 1L + (stats::runif(design$J) >= design$group_share[1L])
  persons <- design$J * design$n
  school <- rep(seq_len(design$J), each = design$n)
  class <- ifelse(stats::runif(persons) < design$improving[group_class[school]], 1L, 2L)
  effects <- cbind(stats::rnorm(persons, 0, design$random_sd[1L]),
                   stats::rnorm(persons, 0, design$random_sd[2L]))
  n_times <- length(design$times)
  frame <- data.frame(school = rep(school, each = n_times),
                      student = rep(seq_len(persons), each = n_times),
                      wave = rep(design$times, persons))
  k <- class[frame$student]
  frame$score <- design$intercept[k] + design$slope[k] * frame$wave +
    effects[frame$student, 1L] + effects[frame$student, 2L] * frame$wave +
    stats::rnorm(nrow(frame), 0, design$sigma[k])
  frame
}

# The truth of every natural parameter in the true labels: class 1 improves,
# group class 1 holds more of it. Membership logits are those of class 2
# against class 1 within each group class, and of group class 2 against 1.
truth_rows <- function(design) {
  within <- rbind(design$improving, 1 - design$improving)  # classes x group classes
  rbind(
    data.frame(block = "coefficient", class = rep(c("class_1", "class_2"), each = 2L),
               term = rep(c("(Intercept)", "wave"), 2L),
               truth = as.vector(rbind(design$intercept, design$slope))),
    data.frame(block = "group_membership", class = "group_class_2", term = "(Intercept)",
               truth = log(design$group_share[2L]) - log(design$group_share[1L])),
    data.frame(block = "membership", class = "class_2",
               term = c("(Intercept):group_class_1", "(Intercept):group_class_2"),
               truth = log(within[2L, ]) - log(within[1L, ])),
    data.frame(block = "residual_sd", class = c("class_1", "class_2"), term = "residual",
               truth = design$sigma),
    data.frame(block = c("random_variance", "random_variance", "random_sd", "random_sd",
                         "random_correlation"),
               class = "all",
               term = c("(Intercept)", "wave", "(Intercept)", "wave", "(Intercept):wave"),
               truth = c(design$random_sd^2, design$random_sd, 0)),
    data.frame(block = "share", class = c("class_1", "class_2"), term = "share",
               truth = as.vector(within %*% design$group_share)),
    data.frame(block = "group_share", class = c("group_class_1", "group_class_2"),
               term = "share", truth = design$group_share))
}

# A fitted table in the true labels: fitted class j is true class
# class_map[j] and fitted group class h true group class group_map[h]. With
# two classes, a logit whose reference class is relabelled changes sign (and
# its interval's ends swap); every other row is only renamed.
relabel <- function(table, class_map, group_map) {
  rename <- function(text, prefix, map) {
    pattern <- paste0(prefix, "([0-9]+)$")
    at <- grepl(pattern, text)
    number <- as.integer(sub(paste0(".*", pattern), "\\1", text[at]))
    text[at] <- paste0(sub(pattern, "", text[at]), prefix, map[number])
    text
  }
  sign <- ifelse(table$block == "membership" & class_map[1L] == 2L |
                   table$block == "group_membership" & group_map[1L] == 2L, -1, 1)
  is_class <- table$block %in% c("coefficient", "residual_sd", "share")
  table$class[is_class] <- rename(table$class[is_class], "class_", class_map)
  is_group <- table$block == "group_share"
  table$class[is_group] <- rename(table$class[is_group], "group_class_", group_map)
  is_logit <- table$block == "membership"
  table$term[is_logit] <- rename(table$term[is_logit], "group_class_", group_map)
  low <- table$conf_low * sign
  high <- table$conf_high * sign
  table$estimate <- table$estimate * sign
  table$conf_low <- pmin(low, high)
  table$conf_high <- pmax(low, high)
  table
}

one_replicate <- function(i) {
  seed <- registry$base_seed + i
  data <- draw(design, seed)
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, c("latents_failed_starts", "latents_degenerate_start",
                      "latents_bootstrap_dropped"))) invokeRestart("muffleWarning")
  })
  failure <- function(reason) {
    data.frame(replicate = i, block = NA_character_, class = NA_character_,
               term = NA_character_, truth = NA_real_, estimate = NA_real_,
               se_wald = NA_real_, low_wald = NA_real_, high_wald = NA_real_,
               se_boot = NA_real_, low_boot = NA_real_, high_boot = NA_real_,
               n_valid = NA_integer_, class_kept = NA_real_, group_kept = NA_real_,
               available = FALSE, reason = reason, stringsAsFactors = FALSE)
  }
  fit <- tryCatch(quiet(mixture_regression(
    score ~ wave, data, n_classes = 2, id = "student", class_level = "group",
    random = ~ 1 + wave, random_covariance = "equal", n_group_classes = 2,
    cluster = "school", n_starts = registry$n_starts, seed = i)),
    error = function(e) e)
  if (inherits(fit, "error")) return(failure(conditionMessage(fit)))
  if (!isTRUE(fit$converged)) return(failure("unconverged"))
  wald <- tryCatch(parameter_inference(fit, vcov_type = "observed"),
                   error = function(e) e)
  if (inherits(wald, "error")) return(failure(conditionMessage(wald)))
  boot <- tryCatch(quiet(parameter_inference(
    fit, method = "bootstrap", iter = iter,
    n_starts = registry$bootstrap_n_starts, seed = i)), error = function(e) e)
  if (inherits(boot, "error")) return(failure(conditionMessage(boot)))
  if (attr(boot, "n_valid") < 0.9 * iter) {
    return(failure(sprintf("only %d usable resamples", attr(boot, "n_valid"))))
  }
  # Name the fitted classes by the true ones: the improving class has the
  # larger slope, the improving group class the larger share of it.
  slope <- wald$estimate[wald$block == "coefficient" & wald$term == "wave"]
  class_map <- if (slope[1L] >= slope[2L]) 1:2 else 2:1
  logits <- wald$estimate[wald$block == "membership"]  # log P(class_2)/P(class_1)
  improving_logit <- if (class_map[1L] == 1L) -logits else logits
  group_map <- if (improving_logit[1L] >= improving_logit[2L]) 1:2 else 2:1
  truth <- truth_rows(design)
  key <- function(table) paste(table$block, table$class, table$term)
  stopifnot(identical(key(wald), key(boot)))
  wald_true <- relabel(wald, class_map, group_map)
  boot_true <- relabel(boot, class_map, group_map)
  at <- match(key(truth), key(wald_true))
  stopifnot(!anyNA(at), identical(key(wald_true), key(boot_true)))
  # Label matching: each matched resample keeps the original's order of the
  # classes (by slope) and of the group classes (by the improving logit).
  replicates <- attr(boot, "replicates")
  columns <- paste(wald$block, wald$class, wald$term, sep = ".")
  colnames(replicates) <- columns
  slopes <- replicates[, c("coefficient.class_1.wave", "coefficient.class_2.wave"),
                       drop = FALSE]
  class_kept <- mean((slopes[, 1L] >= slopes[, 2L]) == (slope[1L] >= slope[2L]))
  group_logits <- replicates[, c("membership.class_2.(Intercept):group_class_1",
                                 "membership.class_2.(Intercept):group_class_2"),
                             drop = FALSE]
  group_kept <- mean((group_logits[, 1L] >= group_logits[, 2L]) ==
                       (logits[1L] >= logits[2L]))
  data.frame(replicate = i, truth[c("block", "class", "term", "truth")],
             estimate = wald_true$estimate[at], se_wald = wald_true$standard_error[at],
             low_wald = wald_true$conf_low[at], high_wald = wald_true$conf_high[at],
             se_boot = boot_true$standard_error[at], low_boot = boot_true$conf_low[at],
             high_boot = boot_true$conf_high[at], n_valid = attr(boot, "n_valid"),
             class_kept = class_kept, group_kept = group_kept,
             available = TRUE, reason = NA_character_, stringsAsFactors = FALSE)
}

started <- Sys.time()
workers <- parallel::makePSOCKcluster(n_cores)
results <- tryCatch({
  parallel::clusterCall(workers, function() {
    pkgload::load_all(".", quiet = TRUE)
    NULL
  })
  parallel::clusterExport(workers, c("registry", "design", "iter", "draw", "truth_rows",
                                     "relabel", "one_replicate"))
  # Every dataset and every bootstrap draws from its own seed, so the
  # results do not depend on the number of workers or their schedule.
  parallel::parLapplyLB(workers, seq_len(n_datasets), function(i) {
    tryCatch(one_replicate(i), error = function(e) {
      stop(sprintf("dataset %d: %s", i, conditionMessage(e)))
    })
  })
}, finally = parallel::stopCluster(workers))
replicates <- do.call(rbind, results)
utils::write.csv(replicates,
                 file.path(out_dir, "growth-simulation-multilevel-bootstrap-replicates.csv"),
                 row.names = FALSE)

# The central 99.9% binomial range of the coverage of n intervals at 0.95.
coverage_range <- function(n) stats::qbinom(c(0.0005, 0.9995), n, 0.95) / n
summarise <- function(replicates, n_datasets) {
  ok <- replicates[replicates$available, ]
  keys <- unique(ok[c("block", "class", "term")])
  available <- length(unique(ok$replicate)) / n_datasets
  do.call(rbind, lapply(seq_len(nrow(keys)), function(r) {
    rows <- ok[ok$block == keys$block[r] & ok$class == keys$class[r] &
                 ok$term == keys$term[r], ]
    n <- nrow(rows)
    range <- coverage_range(n)
    error <- rows$estimate - rows$truth
    mcse <- stats::sd(rows$estimate) / sqrt(n)
    covered <- function(low, high) mean(low <= rows$truth & rows$truth <= high)
    data.frame(keys[r, ], truth = rows$truth[1L], n = n, availability = available,
               bias = mean(error), mcse = mcse, empirical_sd = stats::sd(rows$estimate),
               mean_se_wald = mean(rows$se_wald), mean_se_boot = mean(rows$se_boot),
               median_ratio = stats::median(rows$se_boot / rows$se_wald),
               coverage_wald = covered(rows$low_wald, rows$high_wald),
               coverage_boot = covered(rows$low_boot, rows$high_boot),
               coverage_low = range[1L], coverage_high = range[2L],
               stringsAsFactors = FALSE)
  }))
}
summary <- summarise(replicates, n_datasets)
summary$pass_bias <- abs(summary$bias) <= 3.5 * summary$mcse
summary$pass_coverage_boot <- summary$coverage_boot >= summary$coverage_low &
  summary$coverage_boot <= summary$coverage_high
summary$pass_coverage_wald <- summary$coverage_wald >= summary$coverage_low &
  summary$coverage_wald <= summary$coverage_high
summary$pass_agreement <- summary$median_ratio >= 0.8 & summary$median_ratio <= 1.25
rownames(summary) <- NULL
utils::write.csv(summary,
                 file.path(out_dir, "growth-simulation-multilevel-bootstrap-summary.csv"),
                 row.names = FALSE)

ok <- replicates[replicates$available & !duplicated(replicates$replicate), ]
reasons <- table(replicates$reason[!replicates$available], useNA = "no")
sink(file.path(out_dir, "growth-simulation-multilevel-bootstrap.txt"))
cat(sprintf(paste0("Multilevel growth bootstrap simulation: %d datasets (registry %d), ",
                   "%d resamples each, %s elapsed on %d workers\n\n"),
            n_datasets, registry$n_datasets, iter,
            format(round(Sys.time() - started, 1)), n_cores))
print(summary, digits = 3, row.names = FALSE)
cat(sprintf("\nAvailability: %.3f\n", length(unique(ok$replicate)) / n_datasets))
cat(sprintf("Label matching, share of usable resamples keeping the order: classes %.4f (min %.3f), group classes %.4f (min %.3f)\n",
            mean(ok$class_kept), min(ok$class_kept), mean(ok$group_kept),
            min(ok$group_kept)))
cat(sprintf("Usable resamples per dataset: mean %.1f, min %d\n", mean(ok$n_valid),
            min(ok$n_valid)))
cat("\nUnavailable datasets by reason\n")
print(reasons)
cat("\n")
print(utils::sessionInfo())
sink()
cat(readLines(file.path(out_dir, "growth-simulation-multilevel-bootstrap.txt")), sep = "\n")
