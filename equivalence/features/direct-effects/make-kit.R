# Direct effects of covariates on continuous indicators: the external kit.
#
# Run from the project root (deterministic, rerunnable):
#   Rscript equivalence/features/direct-effects/make-kit.R
#   JOBS=2 sh equivalence/features/run.sh direct-effects all
# The outputs land in equivalence/features/direct-effects/returned/, which
# equivalence/test-features-direct-effects.R reads offline.
#
# Profile data are drawn with simulab::simulate_lpa(); simulab has no
# direct-effect generator, so each covariate is drawn here from its own seed
# and b z is added to the indicator. Every case stores its data, its
# specification, the truth (simulab's parameter table plus the effects and the
# covariates' distribution) and its seed in kit/cases.rds. Every model is
# within the Mplus 9 Demo limits (4 dependent, at most 2 independent).

stopifnot("simulab is needed to build the kit" =
            requireNamespace("simulab", quietly = TRUE))
feature_dir <- file.path("equivalence", "features", "direct-effects")
kit_dir <- file.path(feature_dir, "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(kit_dir, full.names = TRUE))

vars <- paste0("y", 1:4)
means_2 <- rbind(c(0, 0, 0, 0), c(1.8, 1.5, 1.6, 1.2))
means_2_wide <- rbind(c(0, 0, 0, 0), c(2.6, 2.2, 2.4, 2.0))
means_3 <- rbind(c(0, 0, 0, 0), c(1.8, 1.4, -1.0, 0.8), c(-1.4, 1.6, 1.4, -1.0))
sds_2 <- rbind(c(1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1))
sds_3 <- rbind(c(1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1), c(1.2, 0.9, 1, 0.8))
sds_equal <- c(1, 1.2, 0.8, 1)

#' A configuration: how its profiles, covariates and effects are drawn
configuration <- function(k, model, n, means, sds, proportions, effects,
                          covariates = list(z = c(mean = 0.5, sd = 1.5)),
                          kind = "normal", noise = "normal", replicates = 5L,
                          mplus_scale = 1) {
  list(k = k, model = model, n = n, means = means, sds = sds,
       proportions = proportions, effects = effects, covariates = covariates,
       kind = kind, noise = noise, replicates = replicates, mplus_scale = mplus_scale)
}
one_effect <- data.frame(indicator = "y3", covariate = "z", b = 0.8)
two_covariates <- data.frame(indicator = "y3", covariate = c("z1", "z2"), b = c(0.8, -0.5))
two_indicators <- data.frame(indicator = c("y3", "y1"), covariate = "z", b = c(0.8, -0.6))
pair <- list(z1 = c(mean = 0.5, sd = 1.5), z2 = c(mean = -1, sd = 1))

configurations <- list(
  # Regular configurations, five replicates each.
  vvi2 = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), one_effect),
  eei2 = configuration(2L, "EEI", 500L, means_2, sds_equal, c(0.55, 0.45), one_effect),
  vvi3 = configuration(3L, "VVI", 500L, means_3, sds_3, c(0.45, 0.35, 0.2), one_effect),
  eei3 = configuration(3L, "EEI", 500L, means_3, sds_equal, c(0.45, 0.35, 0.2), one_effect),
  vvi2cov2 = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), two_covariates,
                           covariates = pair),
  vvi2ind2 = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), two_indicators),
  vvi2small = configuration(2L, "VVI", c(40L, 50L, 60L, 70L, 80L), means_2_wide,
                            sds_2, c(0.55, 0.45), one_effect),
  eei2small = configuration(2L, "EEI", c(40L, 50L, 60L, 70L, 80L), means_2_wide,
                            sds_equal, c(0.55, 0.45), one_effect),
  vvi2large = configuration(2L, "VVI", 5000L, means_2, sds_2, c(0.55, 0.45), one_effect),
  eei3large = configuration(3L, "EEI", 5000L, means_3, sds_equal, c(0.45, 0.35, 0.2),
                            two_covariates, covariates = pair),
  # One very large case.
  vvi2huge = configuration(2L, "VVI", 20000L, means_2, sds_2, c(0.55, 0.45), one_effect,
                           replicates = 1L),
  # Degenerate and edge configurations.
  farmean = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), one_effect,
                          covariates = list(z = c(mean = 1000, sd = 1.5)), replicates = 2L),
  unitsbig = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45),
                           transform(one_effect, b = 0.8 / 1e4),
                           covariates = list(z = c(mean = 0.5e4, sd = 1.5e4)), replicates = 1L,
                           mplus_scale = 1e-4),
  unitssmall = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45),
                             transform(one_effect, b = 0.8 / 1e-4),
                             covariates = list(z = c(mean = 0.5e-4, sd = 1.5e-4)),
                             replicates = 1L, mplus_scale = 1e4),
  zero = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45),
                       transform(one_effect, b = 0), replicates = 2L),
  collinear = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), one_effect,
                            kind = "collinear", replicates = 2L),
  binary = configuration(2L, "VVI", 500L, means_2, sds_2, c(0.55, 0.45), one_effect,
                         kind = "binary", replicates = 2L),
  heavy = configuration(2L, "VVI", 500L, means_2_wide, sds_2, c(0.55, 0.45), one_effect,
                        noise = "t3", replicates = 2L),
  # One profile: the model is a linear regression of each indicator.
  k1 = configuration(1L, "VVI", 300L, rbind(c(1, 2, 0.5, -1), c(1, 2, 0.5, -1)),
                     sds_equal, c(0.5, 0.5), two_indicators, replicates = 1L),
  k1cov2 = configuration(1L, "EEI", 300L, rbind(c(1, 2, 0.5, -1), c(1, 2, 0.5, -1)),
                         sds_equal, c(0.5, 0.5), two_covariates, covariates = pair,
                         replicates = 1L),
  # A profile of very few rows (about 3% of 400).
  fewrows = configuration(2L, "EEI", 400L, rbind(c(0, 0, 0, 0), c(4, 3.5, 4, 3)), sds_equal,
                          c(0.97, 0.03), one_effect, replicates = 2L))

#' Draw one case of a configuration
draw_case <- function(name, config, replicate, seed) {
  n <- if (length(config$n) > 1L) config$n[[replicate]] else config$n
  simulated <- simulab::simulate_lpa(n = n, means = config$means, sds = config$sds,
                                     proportions = config$proportions, seed = seed)
  frame <- as.data.frame(simulated)
  profile <- match(frame$profile, sort(unique(frame$profile)))
  parameters <- attr(simulated, "simulab_tables")$parameters
  y <- as.matrix(frame[paste0("indicator_", 1:4)])
  colnames(y) <- vars
  set.seed(seed + 100000L)
  if (identical(config$noise, "t3")) {
    # Heavy tails: keep simulab's profiles and means, replace the normal noise
    # by t(3) noise with the same standard deviation.
    sds <- if (is.matrix(config$sds)) config$sds[profile, ] else
      matrix(config$sds, n, 4L, byrow = TRUE)
    y <- config$means[profile, , drop = FALSE] +
      sds * matrix(stats::rt(n * 4L, df = 3) / sqrt(3), n, 4L)
    colnames(y) <- vars
  }
  covariates <- vapply(names(config$covariates), function(covariate) {
    spread <- config$covariates[[covariate]]
    switch(config$kind,
           normal = stats::rnorm(n, spread[["mean"]], spread[["sd"]]),
           binary = as.numeric(stats::runif(n) < 0.4),
           # Nearly collinear with membership: the covariate is mostly the
           # profile, with a little noise.
           collinear = 2 * (profile == 2L) + stats::rnorm(n, 0, 0.3))
  }, numeric(n))
  covariates <- matrix(covariates, n, length(config$covariates),
                       dimnames = list(NULL, names(config$covariates)))
  shift <- vapply(vars, function(indicator) {
    rows <- config$effects$indicator == indicator
    as.vector(covariates[, config$effects$covariate[rows], drop = FALSE] %*%
                config$effects$b[rows])
  }, numeric(n))
  y <- y + shift
  data <- data.frame(id = seq_len(n), covariates, y)
  # The files carry ten decimals; the stored data are exactly what they hold.
  data[-1L] <- lapply(data[-1L], \(column) round(column, 10))
  effects <- config$effects
  list(name = name, seed = seed, data = data,
       spec = list(config = sub("_[0-9]+$", "", name), k = config$k, model = config$model,
                   vars = vars, covariates = names(config$covariates),
                   # Mplus refuses a variable whose variance is above 1e6 or
                   # near zero; its DEFINE multiplies the covariate by this.
                   mplus_scale = config$mplus_scale,
                   effects = effects[c("indicator", "covariate")],
                   direct_effects = vapply(unique(effects$indicator), \(indicator)
                     sprintf("%s ~ %s", indicator, paste(effects$covariate[
                       effects$indicator == indicator], collapse = " + ")),
                     character(1), USE.NAMES = FALSE)),
       truth = list(parameters = parameters, effects = effects,
                    covariates = list(kind = config$kind, distribution = config$covariates),
                    noise = config$noise, membership = profile))
}

seeds <- seq(20261020L, by = 7919L, length.out = 200L)
counter <- 0L
cases <- unlist(lapply(names(configurations), function(config_name) {
  config <- configurations[[config_name]]
  lapply(seq_len(config$replicates), function(replicate) {
    counter <<- counter + 1L
    name <- sprintf("%s_%02d", config_name, replicate)
    draw_case(name, config, replicate, seeds[[counter]])
  })
}), recursive = FALSE)
names(cases) <- vapply(cases, `[[`, character(1), "name")

# -- Latent GOLD -------------------------------------------------------------

latentgold_syntax <- function(case) {
  spec <- case$spec
  equations <- vapply(spec$vars, function(indicator) {
    covariates <- spec$effects$covariate[spec$effects$indicator == indicator]
    sprintf("   %s <- 1 | Cluster%s;", indicator,
            paste0(" + ", covariates, collapse = "", recycle0 = TRUE))
  }, character(1))
  variances <- if (identical(spec$model, "VVI")) sprintf("   %s | Cluster;", spec$vars) else
    character()
  lines <- c(
    "//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", case$name), "",
    "model", sprintf("title '%s';", case$name), "options", "   maxthreads=2;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261020 sets=100 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;", "   output",
    sprintf("      parameters=first standarderrors profile write='%s_write.txt';", case$name),
    sprintf("   outfile '%s_posteriors.txt' classification keep id;", case$name),
    "variables",
    sprintf("   independent %s;", paste(sprintf("%s numeric", spec$covariates), collapse = ", ")),
    sprintf("   dependent %s;", paste(sprintf("%s continuous", spec$vars), collapse = ", ")),
    sprintf("   latent Cluster nominal %d;", spec$k),
    "equations", "   Cluster <- 1;", equations, variances, "end model")
  writeLines(lines, file.path(kit_dir, paste0(case$name, ".lgs")), sep = "\r\n")
}

# -- Mplus -------------------------------------------------------------------

mplus_syntax <- function(case) {
  spec <- case$spec
  y <- paste(spec$vars, collapse = " ")
  regressions <- vapply(unique(spec$effects$indicator), function(indicator) {
    sprintf(" %s ON %s;", indicator,
            paste(spec$effects$covariate[spec$effects$indicator == indicator], collapse = " "))
  }, character(1), USE.NAMES = FALSE)
  # Mplus frees the residual covariance of two dependent variables that are
  # both regressed on covariates; the profile model holds them independent.
  regressed <- toupper(unique(spec$effects$indicator))
  independence <- if (length(regressed) > 1L) {
    pairs <- utils::combn(regressed, 2L)
    sprintf(" %s WITH %s@0;", pairs[1L, ], pairs[2L, ])
  }
  class_blocks <- unlist(lapply(seq_len(spec$k), function(j) {
    c(sprintf(" %%c#%d%%", j), sprintf(" [%s];", y),
      if (identical(spec$model, "VVI")) sprintf(" %s;", y),
      if (identical(spec$model, "VVI")) independence)
  }))
  lines <- c(
    sprintf("TITLE: latents direct effects %s;", case$name),
    sprintf("DATA: FILE = %s.dat;", case$name),
    "VARIABLE:",
    sprintf(" NAMES = id %s %s;", paste(spec$covariates, collapse = " "), y),
    sprintf(" USEVARIABLES = %s %s;", y, paste(spec$covariates, collapse = " ")),
    " IDVARIABLE = id;", sprintf(" CLASSES = c(%d);", spec$k),
    if (spec$mplus_scale != 1) {
      c("DEFINE:", sprintf(" %s = %s * %s;", spec$covariates, spec$covariates,
                           format(spec$mplus_scale, scientific = FALSE)))
    },
    "ANALYSIS:", " TYPE = MIXTURE;", " ESTIMATOR = ML;", " INFORMATION = OBSERVED;",
    " STARTS = 400 80;", " STITERATIONS = 40;", " STSEED = 20261020;",
    " MITERATIONS = 10000;", " LOGCRITERION = 0.00000001;",
    " RLOGCRITERION = 0.000000001;", " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", regressions, independence,
    # EEI: means by class, variances shared (Mplus's default).
    if (identical(spec$model, "EEI")) sprintf(" %s;", y),
    class_blocks,
    "OUTPUT: TECH1;", "SAVEDATA:",
    sprintf(" FILE = %s-posteriors.dat;", case$name), " SAVE = CPROBABILITIES;",
    " FORMAT = F20.12;", sprintf(" RESULTS = %s-results.dat;", case$name))
  writeLines(lines, file.path(kit_dir, paste0(case$name, ".inp")))
}

write_data <- function(case) {
  data <- case$data
  text <- data
  text[-1L] <- lapply(data[-1L], \(column) sprintf("%.10f", column))
  utils::write.table(text, file.path(kit_dir, paste0(case$name, "-lg.dat")), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  utils::write.table(text, file.path(kit_dir, paste0(case$name, ".dat")), sep = " ",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  # The stored data must be the values the programs read.
  back <- utils::read.delim(file.path(kit_dir, paste0(case$name, "-lg.dat")))
  stopifnot("written data differ from the stored data" =
              isTRUE(all.equal(back, data, tolerance = 0, check.attributes = FALSE)))
}

invisible(lapply(cases, \(case) {
  write_data(case)
  latentgold_syntax(case)
  mplus_syntax(case)
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
message(sprintf("%d cases written to %s", length(cases), kit_dir))
