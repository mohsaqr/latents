# Extended Monte Carlo equivalence: model types the first kit does not cover.
#
# - ml_lpa_eee, ml_lpa_vvv: two-level profiles with a full within-profile
#   covariance, invariant (EEE) or profile-specific (VVV); Mplus and Latent GOLD.
# - lta_fiml: a latent transition model (two states, four occasions,
#   homogeneous transitions) with missing indicators, some occasions entirely
#   unobserved; Latent GOLD's Markov module.
# - lta_equal: the same model with one variance per indicator shared by the
#   states (`variance_model = "equal"`); Latent GOLD.
#
# Run from the project root:
#   Rscript equivalence/montecarlo/extended/make-kit.R
#   sh equivalence/montecarlo/extended/run.sh

stopifnot(requireNamespace("simulab", quietly = TRUE))
directory <- file.path("equivalence", "montecarlo", "extended")
kit_dir <- file.path(directory, "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
n_replicates <- 10L

two_level <- function(seed, correlations, sds) {
  d <- simulab::simulate_ml_lpa(
    clusters = 60, cluster_size = 15, means = rbind(c(0, 0, 0), c(1.5, 1.2, -1)),
    profile_probabilities = rbind(c(0.75, 0.25), c(0.25, 0.75)),
    sds = sds, cluster_class_proportions = c(0.55, 0.45),
    correlations = correlations, seed = seed)
  data <- as.data.frame(d)
  out <- data.frame(id = data$id, g = data$cluster, data[paste0("indicator_", 1:3)])
  names(out) <- c("id", "g", "y1", "y2", "y3")
  list(data = out, truth = attr(d, "simulab_tables"), membership = data$profile)
}

transitions <- function(seed, equal) {
  set.seed(seed)
  n_sequences <- 150L
  n_occasions <- 4L
  first <- sample(1:2, n_sequences, replace = TRUE, prob = c(0.6, 0.4))
  states <- Reduce(function(previous, occasion) {
    stay <- ifelse(previous == 1L, 0.8, 0.7)
    ifelse(stats::runif(n_sequences) < stay, previous, 3L - previous)
  }, seq_len(n_occasions - 1L), accumulate = TRUE, init = first)
  state <- as.vector(t(do.call(cbind, states)))
  n <- n_sequences * n_occasions
  sd2 <- if (equal) c(1, 1) else c(1, 0.7)
  data <- data.frame(id = seq_len(n), g = rep(seq_len(n_sequences), each = n_occasions),
                     t = rep(seq_len(n_occasions), n_sequences),
                     y1 = stats::rnorm(n, ifelse(state == 1L, 0, 1.8)),
                     y2 = stats::rnorm(n, ifelse(state == 1L, 0, 1.4), sd2[state]))
  if (!equal) {
    # Missing at random on y1, and 4% of occasions entirely unobserved.
    data$y2[stats::runif(n) < stats::plogis(-2 + 0.8 * data$y1)] <- NA
    data[stats::runif(n) < 0.04, c("y1", "y2")] <- NA
  }
  list(data = data, membership = state)
}

generators <- list(
  ml_lpa_eee = function(seed) {
    r <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3)
    two_level(seed, list(r, r), sds = rbind(c(1, 1.1, 0.9), c(1, 1.1, 0.9)))
  },
  ml_lpa_vvv = function(seed) {
    two_level(seed, list(matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3),
                         matrix(c(1, -0.3, 0.1, -0.3, 1, 0.4, 0.1, 0.4, 1), 3)),
              sds = rbind(c(1, 1, 1), c(0.8, 1.2, 0.9)))
  },
  lta_fiml = function(seed) transitions(seed, equal = FALSE),
  lta_equal = function(seed) transitions(seed, equal = TRUE))

specification <- list(
  ml_lpa_eee = list(level = "two", k = 2L, h = 2L, vars = paste0("y", 1:3), model = "EEE",
                    programs = c("mplus", "latentgold")),
  ml_lpa_vvv = list(level = "two", k = 2L, h = 2L, vars = paste0("y", 1:3), model = "VVV",
                    programs = c("mplus", "latentgold")),
  lta_fiml = list(level = "transition", k = 2L, vars = c("y1", "y2"), variance = "varying",
                  missing = TRUE, programs = "latentgold"),
  lta_equal = list(level = "transition", k = 2L, vars = c("y1", "y2"), variance = "equal",
                   programs = "latentgold"))

mplus_input <- function(name, spec, columns) {
  vars <- paste(spec$vars, collapse = " ")
  covariances <- unlist(lapply(seq_len(length(spec$vars) - 1L), function(i)
    sprintf(" %s WITH %s;", spec$vars[i],
            paste(spec$vars[(i + 1L):length(spec$vars)], collapse = " "))))
  class_specific <- identical(spec$model, "VVV")
  class_block <- unlist(lapply(seq_len(spec$k), function(j)
    c(sprintf(" %%cw#%d%%", j), sprintf(" [%s];", vars),
      if (class_specific) c(sprintf(" %s;", vars), covariances))))
  c(sprintf("TITLE: latents extended Monte Carlo %s;", name),
    sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", vars), " IDVARIABLE = id;", " MISSING = ALL(-999);",
    " CLUSTER = g;", sprintf(" WITHIN = %s;", vars),
    sprintf(" CLASSES = cb(%d) cw(%d);", spec$h, spec$k), " BETWEEN = cb;",
    "ANALYSIS:", " TYPE = TWOLEVEL MIXTURE;", " ESTIMATOR = ML;", " STARTS = 200 40;",
    " STITERATIONS = 40;", " STSEED = 20261018;", " MITERATIONS = 10000;",
    " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
    " MCONVERGENCE = 0.000001;",
    "MODEL:", " %WITHIN%", " %OVERALL%", if (!class_specific) covariances,
    " %BETWEEN%", " %OVERALL%", " cw#1 ON cb;",
    "MODEL cw:", " %WITHIN%", class_block,
    "SAVEDATA:", sprintf(" FILE = %s-posteriors.dat;", name),
    " SAVE = CPROBABILITIES;", " FORMAT = F20.12;", sprintf(" RESULTS = %s-results.dat;", name))
}

latentgold_input <- function(name, spec) {
  options <- c("options", "   maxthreads=all;", "   algorithm",
               "      tolerance=1e-010 emtolerance=1e-008 emiterations=5000 nriterations=500;",
               "   startvalues", "      seed=20261018 sets=100 tolerance=1e-005 iterations=250;",
               "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
               "   missing includeall;",
               "   output", "      parameters=first profile probmeans=posterior;")
  dependents <- sprintf("   dependent %s;",
                        paste(sprintf("%s continuous", spec$vars), collapse = ", "))
  if (identical(spec$level, "transition")) {
    variances <- if (identical(spec$variance, "equal")) character() else
      sprintf("   %s | State;", spec$vars)
    return(c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
             "model", sprintf("title '%s';", name), options,
             sprintf("   outfile '%s_posteriors.txt' classification keep id t;", name),
             "variables", "   caseid g;", dependents,
             "   latent", "      State nominal dynamic 2;",
             "equations", "   State[=0] <- 1;", "   State <- 1 | State[-1];",
             sprintf("   %s <- 1 + State;", spec$vars), variances, "end model"))
  }
  covariances <- unlist(lapply(seq_len(length(spec$vars) - 1L), function(i)
    sprintf("   %s <-> %s%s;", spec$vars[i], spec$vars[(i + 1L):length(spec$vars)],
            if (identical(spec$model, "VVV")) " | Cluster" else "")))
  variances <- if (identical(spec$model, "VVV")) sprintf("   %s | Cluster;", spec$vars) else
    character()
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), options,
    sprintf("   outfile '%s_posteriors.txt' classification keep id g;", name),
    "variables", "   groupid g;", dependents, "   latent",
    sprintf("      GClass group nominal %d,", spec$h),
    sprintf("      Cluster nominal %d;", spec$k),
    "equations", "   GClass <- 1;", "   Cluster <- 1 + GClass;",
    sprintf("   %s <- 1 + Cluster;", spec$vars), variances, covariances, "end model")
}

cases <- list()
invisible(lapply(seq_along(generators), function(index) {
  config <- names(generators)[index]
  spec <- specification[[config]]
  lapply(seq_len(n_replicates), function(replicate) {
    seed <- 20261018L + 1000L * index + replicate
    generated <- generators[[config]](seed)
    name <- sprintf("%s_%02d", config, replicate)
    data <- generated$data
    utils::write.table(data, file.path(kit_dir, paste0(name, ".dat")), sep = " ",
                       quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
    utils::write.table(data, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
    if ("mplus" %in% spec$programs) {
      writeLines(mplus_input(name, spec, names(data)), file.path(kit_dir, paste0(name, ".inp")))
    }
    writeLines(latentgold_input(name, spec), file.path(kit_dir, paste0(name, ".lgs")),
               sep = "\r\n")
    cases[[name]] <<- c(list(config = config, replicate = replicate, seed = seed,
                             spec = spec), generated)
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("%d cases written to %s\n", length(cases), kit_dir))
