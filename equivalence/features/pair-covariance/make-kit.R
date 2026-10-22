# Covariances for chosen pairs of continuous indicators: the external kit.
#
# Run from the project root:
#   Rscript equivalence/features/pair-covariance/make-kit.R
#   JOBS=2 sh equivalence/features/run.sh pair-covariance all
# Every dataset is drawn with simulab from known parameters; the
# per-profile correlation matrices are nonzero on the chosen pairs only
# (misspecified cases add one correlation outside them). simulab's truth
# table has no correlations, so each case stores its own truth. Each case is
# written for Latent GOLD 6.1 (`y1 <-> y2 | Cluster`, or without `| Cluster`
# for one shared covariance) and the Mplus 9 Demo (`y1 WITH y2;` in each
# class section, or in %OVERALL% for one shared covariance). Rerunning
# rewrites the same files.

stopifnot("simulab is needed to build the kit" =
            requireNamespace("simulab", quietly = TRUE))
options(scipen = 100)
root <- file.path("equivalence", "features", "pair-covariance")
kit_dir <- file.path(root, "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)

correlation_matrix <- function(d, pairs, values) {
  m <- diag(d)
  invisible(Map(function(pair, value) {
    at <- as.integer(sub("y", "", strsplit(pair, ":", fixed = TRUE)[[1L]]))
    m[at[1L], at[2L]] <<- value
    m[at[2L], at[1L]] <<- value
  }, pairs, values))
  m
}

# One configuration: generating parameters and the model the programs fit.
# `r` holds one row per profile and one column per correlated pair (the
# fitted pairs, then any `extra` pairs of a misspecified case).
configuration <- function(n, means, sds, proportions, pairs, r, variance = "varying",
                          extra = character(), replicates = 5L, min_variance = NULL,
                          floor = NULL, missing = FALSE, starts = c(200L, 40L),
                          sets = 100L, compare_full = FALSE,
                          programs = c("latentgold", "mplus")) {
  k <- nrow(means)
  r <- matrix(r, k, length(pairs) + length(extra), byrow = TRUE)
  list(n = n, means = means, sds = sds, proportions = proportions, pairs = pairs,
       extra = extra, r = r, variance = variance, replicates = replicates,
       min_variance = min_variance, floor = floor, missing = missing,
       starts = starts, sets = sets, compare_full = compare_full, programs = programs)
}

two <- function(a, b) rbind(a, b, deparse.level = 0)

configurations <- list(
  # -- regular ----------------------------------------------------------------
  vvi_one_k2_m = configuration(
    500, two(c(0, 0, 0, 0), c(2, 1.5, 1.2, -1)),
    two(c(1, 1, 1, 1), c(0.8, 1.3, 0.9, 1.1)), c(0.55, 0.45), "y1:y2", c(0.6, -0.4)),
  vvi_two_k3_m = configuration(
    600, rbind(c(0, 0, 0, 0, 0), c(2, 1.5, -1, 1.2, 0.5), c(-1.5, 2, 1.5, -1, 1.5)),
    rbind(c(1, 1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1, 1), c(1.1, 0.9, 1.2, 0.8, 0.9)),
    c(0.4, 0.35, 0.25), c("y1:y2", "y3:y4"), c(0.5, -0.5, -0.3, 0.4, 0.7, 0.2)),
  eei_one_k2_m = configuration(
    500, two(c(0, 0, 0, 0), c(2, 1.6, 1.2, -1.2)),
    two(c(1, 1.2, 0.9, 1.1), c(1, 1.2, 0.9, 1.1)), c(0.5, 0.5), "y1:y2", c(-0.5, -0.5),
    variance = "equal"),
  eei_two_k3_m = configuration(
    600, rbind(c(0, 0, 0, 0, 0), c(2, 1.5, -1, 1.5, 0.8), c(-1.5, 2, 1.5, -1, 1.6)),
    matrix(c(1, 1.1, 0.9, 1, 1.2), 3, 5, byrow = TRUE), c(0.4, 0.35, 0.25),
    c("y1:y2", "y4:y5"), c(0.6, -0.4, 0.6, -0.4, 0.6, -0.4), variance = "equal"),
  vvi_two_k4_m = configuration(
    1000, rbind(c(0, 0, 0, 0, 0), c(3, 2.5, 0, 2, -1), c(-2.5, 2.5, 2.5, -1, 1.5),
                c(1, -2.5, -2, 2.5, 2.5)),
    rbind(c(1, 1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1, 1), c(1.1, 0.9, 1.2, 0.8, 0.9),
          c(0.9, 1, 0.8, 1.2, 1.1)),
    c(0.3, 0.25, 0.25, 0.2), c("y1:y3", "y2:y5"),
    c(0.5, -0.4, -0.5, 0.6, 0.3, 0.3, 0.7, -0.6), starts = c(400L, 80L), sets = 200L),
  vvi_one_k2_s = configuration(
    60, two(c(0, 0, 0), c(3, 3, 2.5)), two(c(1, 1, 1), c(0.8, 1.2, 1)), c(0.5, 0.5),
    "y2:y3", c(0.6, -0.5), sets = 200L),
  eei_two_k2_s = configuration(
    80, two(c(0, 0, 0, 0), c(3, 2.5, 2.5, -2.5)),
    two(c(1, 1, 1, 1), c(1, 1, 1, 1)), c(0.5, 0.5), c("y1:y2", "y3:y4"),
    c(0.5, -0.5, 0.5, -0.5), variance = "equal", sets = 200L),
  vvi_one_k3_l = configuration(
    5000, rbind(c(0, 0, 0, 0), c(2, 1.5, -1.2, 1), c(-1.5, 2, 1.2, -1.2)),
    rbind(c(1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1), c(1.1, 0.9, 1.2, 0.8)),
    c(0.45, 0.35, 0.2), "y1:y2", c(-0.6, -0.3, -0.7)),
  # -- edges ------------------------------------------------------------------
  near_one = configuration(
    500, two(c(0, 0, 0), c(2, 1.5, -1.5)), two(c(1, 1, 1), c(0.8, 1.2, 1)), c(0.5, 0.5),
    "y1:y2", c(0.98, -0.98), replicates = 2L),
  zero_cov = configuration(
    500, two(c(0, 0, 0, 0), c(2, 1.5, 1.2, -1)), two(c(1, 1, 1, 1), c(0.8, 1.3, 0.9, 1.1)),
    c(0.5, 0.5), "y1:y2", c(0, 0), replicates = 2L),
  # y4 has a standard deviation of 0.01 in the first profile: a floor of 0.01
  # on the variance (1e-4) binds.
  floor_unpaired = configuration(
    400, two(c(0, 0, 0, 0), c(2, 1.5, 1.2, 1)), two(c(1, 1, 1, 0.01), c(0.8, 1.3, 0.9, 1)),
    c(0.5, 0.5), "y1:y2", c(0.6, -0.4), replicates = 2L, min_variance = 0.01,
    floor = "variance"),
  # The pair's smaller eigenvalue in the first profile is 0.05^2 * 0.01 =
  # 2.5e-5; a floor of 1e-3 binds on it while both variances stay above it.
  floor_paired = configuration(
    400, two(c(0, 0, 0), c(1, 1, 2)), two(c(0.05, 0.05, 1), c(0.8, 1.2, 1)),
    c(0.5, 0.5), "y1:y2", c(0.99, -0.3), replicates = 2L, min_variance = 1e-3,
    floor = "eigenvalue"),
  small_class = configuration(
    400, rbind(c(0, 0, 0, 0), c(3, 3, -3, 3), c(-3, 3, 3, -3)),
    rbind(c(1, 1, 1, 1), c(0.8, 1.2, 0.9, 1.1), c(1, 1, 1, 1)),
    c(0.88, 0.09, 0.03), "y1:y2", c(0.5, -0.5, 0.5), replicates = 2L, sets = 200L,
    starts = c(400L, 80L)),
  # y1 in units of 1e-4 (variances near 1e-8) paired with y2 in units of 1e4.
  # The default min_variance (1e-6) is an absolute bound that would bind on
  # y1, so these are fitted with min_variance = 1e-14, which does not. Mplus
  # refuses the data ("One or more variables have a variance of zero"): its
  # input is kept, the comparison is with Latent GOLD only.
  scales = configuration(
    500, two(c(0.001, 0, 0, 0), c(0.0013, 15000, 1.5, -1.2)),
    two(c(1e-4, 1e4, 1, 1), c(0.8e-4, 1.3e4, 0.9, 1.1)), c(0.5, 0.5), "y1:y2",
    c(0.5, -0.5), replicates = 2L, min_variance = 1e-14),
  two_vvi = configuration(
    500, two(c(0, 0), c(2.5, 1.5)), two(c(1, 1), c(0.8, 1.3)), c(0.55, 0.45), "y1:y2",
    c(0.6, -0.4), replicates = 2L, compare_full = TRUE),
  two_eei = configuration(
    500, two(c(0, 0), c(2.5, 1.5)), two(c(1, 1.2), c(1, 1.2)), c(0.55, 0.45), "y1:y2",
    c(0.5, 0.5), variance = "equal", replicates = 2L, compare_full = TRUE),
  missing = configuration(
    500, two(c(0, 0, 0, 0), c(2, 1.5, 1.2, -1)), two(c(1, 1, 1, 1), c(0.8, 1.3, 0.9, 1.1)),
    c(0.55, 0.45), "y1:y2", c(0.6, -0.4), replicates = 3L, missing = TRUE),
  # y3 and y4 are correlated in the data; only y1:y2 is fitted.
  misspecified = configuration(
    500, two(c(0, 0, 0, 0), c(2, 1.5, 1.2, -1)), two(c(1, 1, 1, 1), c(0.8, 1.3, 0.9, 1.1)),
    c(0.55, 0.45), "y1:y2", c(0.6, 0.5, -0.4, 0.5), extra = "y3:y4", replicates = 2L),
  huge = configuration(
    20000, two(c(0, 0, 0, 0), c(1.5, 1.2, 1, -0.8)), two(c(1, 1, 1, 1), c(0.8, 1.3, 0.9, 1.1)),
    c(0.55, 0.45), "y1:y2", c(0.6, -0.4), replicates = 1L))

# -- data -------------------------------------------------------------------------

simulate_case <- function(config, seed) {
  d <- ncol(config$means)
  all_pairs <- c(config$pairs, config$extra)
  correlations <- lapply(seq_len(nrow(config$means)), function(profile) {
    correlation_matrix(d, all_pairs, config$r[profile, ])
  })
  simulated <- simulab::simulate_lpa(n = config$n, means = config$means, sds = config$sds,
                                     proportions = config$proportions,
                                     correlations = correlations, seed = seed)
  frame <- as.data.frame(simulated)
  vars <- paste0("y", seq_len(d))
  data <- stats::setNames(frame[c("id", paste0("indicator_", seq_len(d)))], c("id", vars))
  if (config$missing) {
    data[vars] <- simulab::inject_missingness(data[vars], mechanism = "MAR",
                                              proportion = 0.15,
                                              variables = vars[-1L], predictor = "y1",
                                              seed = seed + 1L)
  }
  covariances <- lapply(seq_along(correlations), function(profile) {
    s <- config$sds[profile, ]
    correlations[[profile]] * tcrossprod(s)
  })
  list(data = data,
       truth = list(means = config$means, sds = config$sds, correlations = correlations,
                    covariances = covariances, proportions = config$proportions,
                    pairs = config$pairs, extra = config$extra,
                    simulab = attr(simulated, "simulab_tables")$parameters),
       membership = as.integer(sub("Profile ", "", frame$profile)))
}

# -- Mplus --------------------------------------------------------------------------

mplus_with <- function(pairs) {
  vapply(strsplit(pairs, ":", fixed = TRUE),
         \(p) sprintf(" %s WITH %s;", p[1L], p[2L]), character(1))
}

# `fixed = TRUE` writes the floor case's companion run: the binding bound
# held as an equality in the first class (y4@floor, or the pair's smaller
# eigenvalue equal to the floor), every other parameter free. Where only that
# bound binds, its maximum is latents' constrained maximum.
mplus_input <- function(name, spec, columns, fixed = FALSE) {
  k <- spec$k
  vars <- paste(spec$vars, collapse = " ")
  varying <- identical(spec$variance, "varying")
  # Labels for the floor constraints: a variance per class, or the pair's two
  # variances and covariance per class.
  floor_labels <- function(j) {
    if (identical(spec$floor, "variance")) return(sprintf(" y4 (f%d);", j))
    if (identical(spec$floor, "eigenvalue")) {
      return(c(sprintf(" y1 (a%d);", j), sprintf(" y2 (c%d);", j),
               sprintf(" y1 WITH y2 (b%d);", j)))
    }
    character()
  }
  class_block <- unlist(lapply(seq_len(k), function(j) {
    if (fixed && j == 1L && identical(spec$floor, "variance")) {
      return(c(" %c#1%", sprintf(" [%s];", vars),
               sprintf(" %s;", paste(setdiff(spec$vars, "y4"), collapse = " ")),
               sprintf(" y4@%s;", format(spec$min_variance)), mplus_with(spec$pairs)))
    }
    c(sprintf(" %%c#%d%%", j), sprintf(" [%s];", vars),
      if (varying) c(sprintf(" %s;", vars), mplus_with(spec$pairs)),
      if (varying && (!fixed || identical(spec$floor, "eigenvalue"))) floor_labels(j))
  }))
  constraint <- if (fixed && identical(spec$floor, "eigenvalue")) {
    c("MODEL CONSTRAINT:", sprintf(
      " 0 = (a1 + c1)/2 - SQRT(((a1 - c1)/2)**2 + b1**2) - %s;", format(spec$min_variance)))
  } else if (fixed) {
    character()
  } else if (identical(spec$floor, "variance")) {
    c("MODEL CONSTRAINT:", sprintf(" f%d > %s;", seq_len(k), format(spec$min_variance)))
  } else if (identical(spec$floor, "eigenvalue")) {
    c("MODEL CONSTRAINT:", sprintf(
      " 0 < (a%d + c%d)/2 - SQRT(((a%d - c%d)/2)**2 + b%d**2) - %s;",
      seq_len(k), seq_len(k), seq_len(k), seq_len(k), seq_len(k),
      format(spec$min_variance)))
  } else character()
  c(sprintf("TITLE: latents pair covariance %s;", name),
    sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", vars), " IDVARIABLE = id;", " MISSING = ALL(-999);",
    sprintf(" CLASSES = c(%d);", k),
    "ANALYSIS:", " TYPE = MIXTURE;", " ESTIMATOR = ML;",
    sprintf(" STARTS = %d %d;", spec$starts[1L], spec$starts[2L]),
    " STITERATIONS = 40;", " STSEED = 20261020;", " MITERATIONS = 10000;",
    " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
    " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", if (!varying) mplus_with(spec$pairs), class_block,
    constraint, "OUTPUT: TECH1;", "SAVEDATA:",
    sprintf(" FILE = %s-posteriors.dat;", name), " SAVE = CPROBABILITIES;",
    " FORMAT = F20.12;", sprintf(" RESULTS = %s-results.dat;", name))
}

# -- Latent GOLD ------------------------------------------------------------------

latentgold_input <- function(name, spec) {
  varying <- identical(spec$variance, "varying")
  given <- if (varying) " | Cluster" else ""
  covariances <- vapply(strsplit(spec$pairs, ":", fixed = TRUE),
                        \(p) sprintf("   %s <-> %s%s;", p[1L], p[2L], given), character(1))
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=2;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues",
    sprintf("      seed=20261020 sets=%d tolerance=1e-005 iterations=250;", spec$sets),
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    sprintf("   missing %s;", if (isTRUE(spec$missing)) "includeall" else "excludeall"),
    "   output",
    sprintf("      parameters=first standarderrors profile write='%s_write.txt';", name),
    sprintf("   outfile '%s_posteriors.txt' classification keep id;", name),
    "variables",
    sprintf("   dependent %s;", paste(sprintf("%s continuous", spec$vars), collapse = ", ")),
    "   latent", sprintf("      Cluster nominal %d;", spec$k),
    "equations", "   Cluster <- 1;",
    sprintf("   %s <- 1 + Cluster;", spec$vars),
    sprintf("   %s%s;", spec$vars, given), covariances, "end model")
}

# -- write --------------------------------------------------------------------------

# The floor cases constrain variances or eigenvalues, which Latent GOLD has no
# inequality for; they go to Mplus (MODEL CONSTRAINT) and the R reference.
cases <- list()
invisible(lapply(seq_along(configurations), function(index) {
  config_name <- names(configurations)[index]
  config <- configurations[[config_name]]
  lapply(seq_len(config$replicates), function(replicate) {
    seed <- 20261020L + 1000L * index + replicate
    generated <- simulate_case(config, seed)
    name <- sprintf("%s_%02d", config_name, replicate)
    spec <- list(k = nrow(config$means), vars = paste0("y", seq_len(ncol(config$means))),
                 pairs = config$pairs, variance = config$variance,
                 missing = config$missing, min_variance = config$min_variance,
                 floor = config$floor, starts = config$starts, sets = config$sets,
                 compare_full = config$compare_full,
                 programs = if (!is.null(config$floor)) "mplus" else
                   if (identical(config_name, "scales")) "latentgold" else config$programs)
    data <- generated$data
    utils::write.table(data, file.path(kit_dir, paste0(name, ".dat")), sep = " ",
                       quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
    writeLines(mplus_input(name, spec, names(data)), file.path(kit_dir, paste0(name, ".inp")))
    if (!is.null(spec$floor)) {
      companion <- paste0(name, "_fixed")
      utils::write.table(data, file.path(kit_dir, paste0(companion, ".dat")), sep = " ",
                         quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
      writeLines(mplus_input(companion, spec, names(data), fixed = TRUE),
                 file.path(kit_dir, paste0(companion, ".inp")))
    }
    if ("latentgold" %in% spec$programs) {
      utils::write.table(data, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                         quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
      writeLines(latentgold_input(name, spec), file.path(kit_dir, paste0(name, ".lgs")),
                 sep = "\r\n")
    }
    cases[[name]] <<- c(list(name = name, config = config_name, replicate = replicate,
                             seed = seed, spec = spec), generated)
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("%d cases in %d configurations written to %s\n", length(cases),
            length(configurations), kit_dir))
