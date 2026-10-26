# The ML three-step kit: simulated datasets (simulab), their latents
# measurement fits, and the Latent GOLD / Mplus inputs that compare the
# step-three model on the same posteriors.
#
# From the project root:
#   Rscript equivalence/features/ml-three-step/make-kit.R
#   JOBS=2 sh equivalence/features/run.sh ml-three-step all
#   Rscript equivalence/features/ml-three-step/make-pipeline-step3.R
#   JOBS=2 sh equivalence/features/run.sh ml-three-step latentgold _s3
#
# Every case is stored in kit/cases.rds with its data, spec, truth, seed and
# the posteriors written for Latent GOLD (latents' own, full precision).
# Rerunning rewrites the same files.

source(file.path("validation", "fork-runtime.R"))
invisible(.validation_fork_cores(2L))
stopifnot("simulab is needed to build the kit" = requireNamespace("simulab", quietly = TRUE))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
sys.source(file.path("equivalence", "features", "ml-three-step", "helpers.R"),
           envir = environment())
kit_dir <- ml3_kit()
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(kit_dir, full.names = TRUE))

patterns <- rbind(c(0, 0, 0, 0), c(1, 1, -1, 0.5), c(-1, 1, 1, -0.5), c(1, -1, 0.5, 1),
                  c(-0.5, -1, -1, 1))
lca_patterns <- rbind(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1), c(1, 1, 0, 1, 0, 0),
                      c(0, 1, 1, 0, 1, 1))
proportions <- list(c(0.6, 0.4), c(0.45, 0.35, 0.2), c(0.3, 0.3, 0.2, 0.2),
                    c(0.25, 0.2, 0.2, 0.2, 0.15))
effects <- c(0, 1, 2, -1, 0.5)
varying_sds <- c(1, 1.6, 0.7, 1.3, 0.9)

#' The measurement part: indicators and the true class
measurement <- function(family, k, n, separation, seed) {
  if (identical(family, "lpa")) {
    d <- as.data.frame(simulab::simulate_lpa(
      n = n, means = separation * patterns[seq_len(k), , drop = FALSE], sds = 1,
      proportions = proportions[[k - 1L]], seed = seed))
    names(d) <- sub("^indicator_", "y", names(d))
    return(list(data = data.frame(id = d$id, d[paste0("y", 1:4)]),
                truth_class = as.integer(sub("\\D+", "", d$profile)),
                vars = paste0("y", 1:4),
                parameters = list(means = separation * patterns[seq_len(k), , drop = FALSE],
                                  sds = 1, proportions = proportions[[k - 1L]])))
  }
  endorse <- ifelse(lca_patterns[seq_len(k), , drop = FALSE] == 1, separation, 1 - separation)
  d <- as.data.frame(simulab::simulate_lca(
    n = n, probabilities = array(c(endorse, 1 - endorse), dim = c(k, 6L, 2L)),
    proportions = proportions[[k - 1L]], seed = seed))
  items <- grep("^item_", names(d), value = TRUE)
  out <- data.frame(id = d$id, stats::setNames(lapply(d[items], as.integer),
                                               sub("item_", "u", items)))
  list(data = out, truth_class = as.integer(sub("\\D+", "", d$latent_class)),
       vars = paste0("u", 1:6),
       parameters = list(endorse = endorse, proportions = proportions[[k - 1L]]))
}

#' A distal outcome from the true class
outcome <- function(truth_class, k, kind, seed) {
  set.seed(seed)
  n <- length(truth_class)
  means <- 500 + 60 * effects[seq_len(k)]
  sds <- switch(kind, varying = 50 * varying_sds[seq_len(k)], 50 + 0 * seq_len(k))
  y <- switch(kind,
    varying = , equal = stats::rnorm(n, means[truth_class], sds[truth_class]),
    skewed = 300 + 60 * exp(stats::rnorm(n, effects[seq_len(k)][truth_class] * 0.6, 0.7)),
    heavy = means[truth_class] + 30 * stats::rt(n, df = 3),
    null = stats::rnorm(n, 500, 50))
  list(y = y, truth = list(kind = kind, means = means, sds = sds))
}

new_case <- function(name, config, replicate, seed, family, k, n, separation, kind,
                     level = "individuals", lg = TRUE, extra = list()) {
  m <- measurement(family, k, n, separation, seed)
  o <- outcome(m$truth_class, k, kind, seed + 7L)
  data <- m$data
  data$y <- o$y
  data$truth <- m$truth_class
  utils::modifyList(list(name = name, config = config, replicate = replicate, seed = seed,
                         level = level, outcome = "y", expect = "ok", lg = lg, mplus = FALSE,
                         pipeline = FALSE, psu = FALSE, data = data,
                         spec = list(family = family, k = k, vars = m$vars, model = "VVI",
                                     separation = separation, n = n),
                         truth = c(m$parameters, list(outcome = o$truth))), extra)
}

cases <- list()
add <- function(case) cases[[case$name]] <<- case

# -- regular configurations, five replicates each ------------------------------------

regular <- list(
  lpa_k2_hi = list("lpa", 2L, 500L, 2.2, "varying"),
  lpa_k2_lo = list("lpa", 2L, 500L, 1.2, "equal"),
  lpa_k3_hi = list("lpa", 3L, 500L, 2.0, "varying"),
  lpa_k3_lo = list("lpa", 3L, 600L, 1.15, "varying"),
  lpa_k4_hi = list("lpa", 4L, 800L, 2.2, "equal"),
  lpa_k5_hi = list("lpa", 5L, 1500L, 2.3, "varying"),
  lca_k2_hi = list("lca", 2L, 500L, 0.85, "varying"),
  lca_k3_lo = list("lca", 3L, 800L, 0.70, "equal"),
  small_k2 = list("lpa", 2L, 60L, 2.2, "varying"),
  small_k3 = list("lpa", 3L, 100L, 2.4, "equal"),
  large_lca_k3 = list("lca", 3L, 5000L, 0.85, "varying"))
mplus_subset <- c("lpa_k2_hi_01", "lpa_k2_hi_02", "lpa_k2_lo_01", "lpa_k2_lo_02",
                  "lpa_k3_hi_01", "lpa_k3_hi_02", "lca_k2_hi_01", "lca_k2_hi_02",
                  "lca_k3_lo_01", "lca_k3_lo_02")
pipeline_subset <- c("lpa_k2_hi_01", "lpa_k3_lo_01", "lpa_k4_hi_01", "lca_k2_hi_01",
                     "lca_k3_lo_01", "small_k2_01")
invisible(lapply(seq_along(regular), \(index) {
  r <- regular[[index]]
  lapply(seq_len(5L), \(replicate) {
    name <- sprintf("%s_%02d", names(regular)[index], replicate)
    add(new_case(name, names(regular)[index], replicate,
                 20261020L + 1000L * index + replicate, r[[1L]], r[[2L]], r[[3L]],
                 r[[4L]], r[[5L]],
                 extra = list(mplus = name %in% mplus_subset,
                              pipeline = name %in% pipeline_subset)))
  })
}))
add(new_case("huge_lpa_k2_01", "huge_lpa_k2", 1L, 20262001L, "lpa", 2L, 20000L, 1.3,
             "varying"))

# -- edge cases --------------------------------------------------------------------------

base <- cases$lpa_k2_hi_01
base$mplus <- FALSE
base$pipeline <- FALSE
offset <- base
offset$name <- "e_offset_01"; offset$config <- "e_offset"
offset$data$y <- base$data$y + 1e6; offset$truth$outcome$transform <- "y + 1e6"
add(offset)
units <- base
units$name <- "e_units_01"; units$config <- "e_units"
units$data$y <- base$data$y * 1e-4; units$truth$outcome$transform <- "y * 1e-4"
add(units)
add(new_case("e_skewed_01", "e_skewed", 1L, 20263001L, "lpa", 3L, 500L, 2.0, "skewed"))
add(new_case("e_heavy_01", "e_heavy", 1L, 20263002L, "lpa", 3L, 500L, 2.0, "heavy"))
add(new_case("e_null_01", "e_null", 1L, 20263003L, "lpa", 3L, 600L, 2.0, "null"))
add(new_case("e_perfect_01", "e_perfect", 1L, 20263004L, "lpa", 3L, 500L, 8, "varying"))
constant <- base
constant$name <- "e_constant_01"; constant$config <- "e_constant"
constant$data$y <- 500; constant$expect <- "latents_bad_outcome"; constant$lg <- FALSE
add(constant)
missing_case <- base
missing_case$name <- "e_na_01"; missing_case$config <- "e_na"
missing_case$data$y[c(3L, 50L)] <- NA; missing_case$expect <- "outcome must not be missing"
missing_case$lg <- FALSE
add(missing_case)
weighted <- base
weighted$name <- "e_weighted_01"; weighted$config <- "e_weighted"
weighted$data$w <- rep(c(1, 2), length.out = nrow(base$data))
weighted$spec$extra <- list(weights = "w")
weighted$expect <- "latents_unsupported_weights"; weighted$lg <- FALSE
add(weighted)

# A class to which (almost) no unit is modally assigned: two true classes,
# fitted with three whose means and variances are held at known values, the
# middle one between the others.
middle_case <- function(name, seed, middle_variance) {
  d <- as.data.frame(simulab::simulate_lpa(n = 400, means = rbind(c(0, 0), c(3, 3)), sds = 1,
                                           seed = seed))
  names(d) <- sub("^indicator_", "y", names(d))
  truth_class <- as.integer(sub("\\D+", "", d$profile))
  o <- outcome(truth_class, 2L, "varying", seed + 7L)
  means <- rbind(c(0, 0), c(1.5, 1.5), c(3, 3))
  start <- list(means = means, variances = rbind(c(1, 1), rep(middle_variance, 2), c(1, 1)),
                profile_probabilities = matrix(c(0.45, 0.1, 0.45), 1L),
                group_probabilities = 1)
  list(name = name, config = sub("_[0-9]+$", "", name), replicate = 1L, seed = seed,
       level = "individuals", outcome = "y", expect = "ok", lg = TRUE, mplus = FALSE,
       pipeline = FALSE, psu = FALSE,
       data = data.frame(id = d$id, d[c("y1", "y2")], y = o$y, truth = truth_class),
       spec = list(family = "lpa", k = 3L, vars = c("y1", "y2"), model = "VVI",
                   n_starts = 1L, n = 400L,
                   extra = list(start = start, fixed = "measurement")),
       truth = list(means = means, middle_variance = middle_variance, outcome = o$truth))
}
modal_counts <- function(case) tabulate(ml3_fit(case)$subject_profiles, case$spec$k)
add(middle_case("e_empty_01", 20264001L, 1))
stopifnot("the empty-column case has a modal unit in every class" =
            min(modal_counts(cases$e_empty_01)) == 0L)
few <- Find(\(case) min(modal_counts(case)) %in% 1:3,
            lapply(20264100L + seq_len(40L), \(seed) middle_case("e_few_01", seed, 0.2)))
stopifnot("no seed gives a class with 1 to 3 modal units" = !is.null(few))
add(few)

# Two-level fits: an individual outcome with a group effect (cluster against
# independent variances) and a group outcome by group class (level = "groups").
two_level <- function(replicate) {
  seed <- 20265000L + replicate
  d <- as.data.frame(simulab::simulate_ml_lpa(
    clusters = 100L, cluster_size = 10L, means = rbind(c(0, 0, 0), c(2, 1.8, -1.6)),
    profile_probabilities = rbind(c(0.8, 0.2), c(0.2, 0.8)), sds = 1,
    cluster_class_proportions = c(0.55, 0.45), seed = seed))
  names(d) <- sub("^indicator_", "y", names(d))
  profile <- as.integer(sub("\\D+", "", d$profile))
  group_class <- as.integer(sub("\\D+", "", d$cluster_class))
  set.seed(seed + 7L)
  group_effect <- stats::rnorm(100L, 0, 30)[d$cluster]
  y <- 500 + c(0, 60)[profile] + group_effect + stats::rnorm(nrow(d), 0, c(50, 70)[profile])
  y_group <- stats::rnorm(100L, c(480, 540)[tapply(group_class, d$cluster, `[`, 1L)],
                          c(40, 25)[tapply(group_class, d$cluster, `[`, 1L)])[d$cluster]
  data <- data.frame(id = d$id, g = d$cluster, d[c("y1", "y2", "y3")], y = y,
                     y_group = y_group, truth = profile, truth_group = group_class)
  spec <- list(family = "ml_lpa", k = 2L, h = 2L, vars = c("y1", "y2", "y3"),
               model = "VVI", n = nrow(data))
  common <- list(replicate = replicate, seed = seed, expect = "ok", lg = TRUE,
                 mplus = FALSE, pipeline = FALSE, data = data, spec = spec,
                 truth = list(individual = list(means = c(500, 560), sds = c(50, 70),
                                                group_sd = 30),
                              group = list(means = c(480, 540), sds = c(40, 25))))
  add(c(list(name = sprintf("ml_ind_%02d", replicate), config = "ml_ind",
             level = "individuals", outcome = "y", psu = TRUE), common))
  add(c(list(name = sprintf("ml_grp_%02d", replicate), config = "ml_grp",
             level = "groups", outcome = "y_group", psu = FALSE), common))
}
invisible(lapply(seq_len(3L), two_level))

# -- fit every case, write the inputs ----------------------------------------------------

fitted <- parallel::mclapply(cases, \(case) {
  fit <- suppressWarnings(ml3_fit(case))
  inputs <- ml3_level_inputs(fit, case)
  list(posteriors = inputs$posteriors, modal = inputs$modal, values = inputs$values,
       ids = inputs$ids, clusters = inputs$clusters,
       step1_log_likelihood = fit$log_likelihood, entropy = ml3_entropy(inputs$posteriors),
       modal_counts = tabulate(inputs$modal, ncol(inputs$posteriors)))
}, mc.cores = .validation_fork_cores(2L), mc.preschedule = FALSE)
lost <- vapply(fitted, \(f) !is.list(f), logical(1))
stopifnot("every case must fit" = !any(lost))
cases <- Map(\(case, f) c(case, f), cases, fitted)

invisible(lapply(cases, \(case) {
  if (!isTRUE(case$lg)) return(NULL)
  ml3_write_lg_data(file.path(kit_dir, paste0(case$name, "-data.txt")), case$ids,
                    case$clusters, case$values, case$posteriors)
  writeLines(ml3_lg_step3(case$name, paste0(case$name, "-data.txt"), case$spec$k,
                          ml3_lg_models(case)),
             file.path(kit_dir, paste0(case$name, ".lgs")), sep = "\r\n")
  if (isTRUE(case$mplus)) {
    data <- case$data[c("id", case$spec$vars)]
    data$ydu <- case$data$y
    data$yde <- case$data$y
    utils::write.table(format(data, digits = 15), file.path(kit_dir, paste0(case$name, ".dat")),
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(ml3_mplus_input(case$name, case$spec, names(data)),
               file.path(kit_dir, paste0(case$name, ".inp")))
  }
  if (isTRUE(case$pipeline)) {
    step1 <- sprintf("pl_%s_s1", case$name)
    data <- case$data[c("id", case$spec$vars, "y")]
    lines <- c(paste(names(data), collapse = "\t"),
               apply(format(data, digits = 17, trim = TRUE), 1L, paste, collapse = "\t"))
    writeLines(lines, file.path(kit_dir, paste0(step1, "-data.txt")), sep = "\r\n")
    writeLines(ml3_lg_step1(step1, paste0(step1, "-data.txt"), case$spec),
               file.path(kit_dir, paste0(step1, ".lgs")), sep = "\r\n")
  }
  NULL
}))

# -- the reference sweep: cheap datasets compared with reference.R only -------------------

sweep_grid <- expand.grid(family = c("lpa", "lca"), k = 2:4, n = c(50L, 200L, 1000L),
                          separation = c("hi", "lo"), stringsAsFactors = FALSE)
sweep_grid <- rbind(sweep_grid, sweep_grid, sweep_grid,
                    expand.grid(family = c("lpa", "lca"), k = 2:4, n = 3000L,
                                separation = c("hi", "lo"), stringsAsFactors = FALSE))
sweep_grid$kind <- rep(c("varying", "equal", "skewed", "null", "heavy"),
                       length.out = nrow(sweep_grid))
separation_of <- function(family, separation) {
  switch(paste(family, separation), "lpa hi" = 2.2, "lpa lo" = 1.2, "lca hi" = 0.85,
         "lca lo" = 0.72)
}
sweep <- lapply(seq_len(nrow(sweep_grid)), \(i) {
  row <- sweep_grid[i, ]
  case <- new_case(sprintf("sweep_%03d", i), "sweep", i, 20266000L + i, row$family, row$k,
                   row$n, separation_of(row$family, row$separation), row$kind, lg = FALSE)
  case$spec$n_starts <- 5L
  case
})
names(sweep) <- vapply(sweep, `[[`, character(1), "name")
sweep_fitted <- parallel::mclapply(sweep, \(case) {
  fit <- suppressWarnings(ml3_fit(case))
  list(posteriors = unname(fit$subject_posteriors), modal = fit$subject_profiles,
       step1_log_likelihood = fit$log_likelihood,
       entropy = ml3_entropy(fit$subject_posteriors),
       modal_counts = tabulate(fit$subject_profiles, case$spec$k))
}, mc.cores = .validation_fork_cores(2L), mc.preschedule = FALSE)
stopifnot("every sweep case must fit" = all(vapply(sweep_fitted, is.list, logical(1))))
sweep <- Map(\(case, f) c(case, f), sweep, sweep_fitted)

saveRDS(c(cases, sweep), file.path(kit_dir, "cases.rds"), compress = "xz")
cat(sprintf("%d kit cases (%d for Latent GOLD, %d for Mplus, %d pipelines) and %d sweep cases\n",
            length(cases), sum(vapply(cases, \(c) isTRUE(c$lg), logical(1))),
            sum(vapply(cases, \(c) isTRUE(c$mplus), logical(1))),
            sum(vapply(cases, \(c) isTRUE(c$pipeline), logical(1))), length(sweep)))
summary_table <- data.frame(
  name = names(cases), level = vapply(cases, `[[`, character(1), "level"),
  k = vapply(cases, \(c) c$spec$k, integer(1)), n = vapply(cases, \(c) nrow(c$posteriors), integer(1)),
  entropy = round(vapply(cases, `[[`, numeric(1), "entropy"), 3),
  min_modal = vapply(cases, \(c) min(c$modal_counts), integer(1)), row.names = NULL)
print(summary_table)
