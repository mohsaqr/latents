# The Vuong-Lo-Mendell-Rubin kit: simulated enumerations written for Latent
# GOLD 6.1 (one listing per case, `latent Cluster nominal 1:4`, which prints
# the VLMR block for every neighbouring pair) and the Mplus 9 Demo (one input
# per number of classes K, TECH11 testing K - 1 against K).
#
# Run from the project root (deterministic, rerunnable):
#   Rscript equivalence/features/vlmr/make-kit.R
#   JOBS=2 sh equivalence/features/run.sh vlmr all
# Every dataset comes from simulab, which keeps the generating parameters.

stopifnot("simulab is needed to build the kit" =
            requireNamespace("simulab", quietly = TRUE))
kit_dir <- file.path("equivalence", "features", "vlmr", "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(kit_dir, full.names = TRUE))

renamed <- function(data, prefix = "y") {
  data <- as.data.frame(data)
  names(data) <- sub("^(indicator_|V)", prefix, names(data))
  data
}

lpa_case <- function(n, means, sds, proportions, correlations = NULL) {
  function(seed) {
    d <- simulab::simulate_lpa(n = n, means = means, sds = sds, proportions = proportions,
                               correlations = correlations, seed = seed)
    vars <- paste0("indicator_", seq_len(ncol(means)))
    list(data = renamed(as.data.frame(d)[c("id", vars)]),
         truth = attr(d, "simulab_tables")$parameters,
         membership = as.data.frame(d)$profile)
  }
}

lca_case <- function(n, endorse, proportions) {
  function(seed) {
    probabilities <- array(c(endorse, 1 - endorse), dim = c(nrow(endorse), ncol(endorse), 2))
    d <- simulab::simulate_lca(n = n, probabilities = probabilities,
                               proportions = proportions, seed = seed)
    data <- as.data.frame(d)
    items <- grep("^item_", names(data), value = TRUE)
    list(data = data.frame(id = data$id, stats::setNames(lapply(data[items], as.integer),
                                                          sub("item_", "u", items))),
         truth = attr(d, "simulab_tables")$parameters, membership = data$latent_class)
  }
}

two_means <- rbind(c(0, 0, 0, 0), c(1.6, 1.4, 1.2, 1.0))
poor_means <- rbind(c(0, 0, 0, 0), c(1.0, 0.8, -0.6, 0.5), c(-0.8, 0.9, 0.8, -0.6))
three_means <- rbind(c(0, 0, 0), c(2, 1.5, -1.5), c(-1.5, 2, 1.5))
shared_correlation <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3)
endorse_two <- rbind(c(0.85, 0.80, 0.75, 0.25, 0.20), c(0.20, 0.30, 0.25, 0.80, 0.85))
endorse_three <- rbind(c(0.85, 0.80, 0.75, 0.20, 0.25, 0.30),
                       c(0.20, 0.25, 0.30, 0.80, 0.75, 0.85),
                       c(0.80, 0.75, 0.85, 0.80, 0.80, 0.75))

# One entry per configuration: its generator, the fitted structure, the
# enumeration range, the true number of classes and the replicates.
configurations <- list(
  eei_k2 = list(make = lpa_case(500, two_means, 1, c(0.55, 0.45)), model = "EEI",
                range = 1:4, true_k = 2L, replicates = 5L,
                condition = "null true for 2 vs 3, well separated, n 500"),
  eei_k3poor = list(make = lpa_case(500, poor_means, 1, c(0.4, 0.35, 0.25)), model = "EEI",
                    range = 1:4, true_k = 3L, replicates = 5L,
                    condition = "alternative true, poorly separated, n 500"),
  eei_small = list(make = lpa_case(80, two_means[, 1:3] * 1.5, 1, c(0.5, 0.5)), model = "EEI",
                   range = 1:4, true_k = 2L, replicates = 5L,
                   condition = "small n 80, null true for 2 vs 3"),
  eei_large = list(make = lpa_case(5000, three_means * 0.6, 1, c(0.45, 0.35, 0.2)),
                   model = "EEI", range = 2:4, true_k = 3L, replicates = 5L,
                   condition = "large n 5000, alternative true for 2 vs 3"),
  eei_k1 = list(make = function(seed) {
    d <- simulab::simulate_correlated(n = 300, means = c(0, 1, -1, 0.5), sds = c(1, 2, 0.5, 1),
                                      seed = seed)
    list(data = renamed(d), truth = attr(d, "simulab_tables")$parameters,
         membership = rep("Profile 1", 300))
  }, model = "EEI", range = 1:4, true_k = 1L, replicates = 5L,
  condition = "one true class: 1 vs 2 under the null"),
  vvi_k2 = list(make = lpa_case(500, two_means, rbind(c(1, 1, 1, 1), c(0.7, 1.3, 0.9, 1.1)),
                                c(0.55, 0.45)), model = "VVI", range = 1:4, true_k = 2L,
                replicates = 5L, condition = "VVI, null true for 2 vs 3"),
  vvi_k3 = list(make = lpa_case(500, three_means, rbind(c(1, 1, 1), c(0.7, 1.2, 0.9),
                                                        c(1.2, 0.8, 1)), c(0.4, 0.35, 0.25)),
                model = "VVI", range = 1:4, true_k = 3L, replicates = 5L,
                condition = "VVI, alternative true, well separated"),
  eee_k2 = list(make = lpa_case(500, two_means[, 1:3] * 1.2, c(1, 1.2, 0.8), c(0.6, 0.4),
                                list(shared_correlation, shared_correlation)),
                model = "EEE", range = 1:4, true_k = 2L, replicates = 5L,
                condition = "EEE shared covariance, null true for 2 vs 3"),
  lca_k2 = list(make = lca_case(500, endorse_two, c(0.55, 0.45)), model = "LCA",
                range = 1:4, true_k = 2L, replicates = 5L,
                condition = "LCA 5 binary items, null true for 2 vs 3"),
  lca_k3 = list(make = lca_case(1000, endorse_three, c(0.4, 0.35, 0.25)), model = "LCA",
                range = 1:4, true_k = 3L, replicates = 5L,
                condition = "LCA 6 binary items, alternative true"),
  weights = list(make = function(seed) {
    generated <- lpa_case(400, two_means, 1, c(0.55, 0.45))(seed)
    # Weights 0.5 and 1.5 in equal numbers sum to n, so frequency and
    # normalised sampling weights give the same log likelihood.
    set.seed(seed)
    generated$data$w <- sample(rep(c(0.5, 1.5), 200))
    generated
  }, model = "EEI", range = 1:4, true_k = 2L, replicates = 3L, weights = "w",
  condition = "sampling weights summing to n"),
  missing = list(make = function(seed) {
    generated <- lpa_case(500, two_means, 1, c(0.55, 0.45))(seed)
    vars <- paste0("y", 1:4)
    generated$data[vars] <- simulab::inject_missingness(generated$data[vars], mechanism = "MCAR",
                                                        proportion = 0.1, seed = seed + 1L)
    # The third replicate also has 20 rows with nothing observed.
    if (seed %% 10L == 3L) generated$data[seq_len(20L), vars] <- NA
    generated
  }, model = "EEI", range = 1:4, true_k = 2L, replicates = 3L, missing = TRUE,
  condition = "MCAR 10% (FIML); replicate 3 with 20 empty rows"),
  vvi_small = list(make = lpa_case(60, two_means[, 1:3] * 1.5,
                                   rbind(c(1, 1, 1), c(0.6, 1.4, 0.8)), c(0.5, 0.5)),
                   model = "VVI", range = 1:4, true_k = 2L, replicates = 3L,
                   condition = "VVI n 60: over-extracted fits near the boundary"),
  eei_tiny = list(make = lpa_case(500, rbind(c(0, 0, 0), c(1.8, 1.5, 1.2), c(5, -4, 5)), 1,
                                  c(0.6, 0.37, 0.03)),
                  model = "EEI", range = 1:4, true_k = 3L, replicates = 3L,
                  condition = "a 3% class far from the rest"))

# Reference-only datasets: fitted by latents and checked against the
# independent reference (reference.R), not run through the two programs.
# They widen the range of sizes (50 to 4000) and structures.
internal_configurations <- list(
  ref_eei_n50 = list(make = lpa_case(50, two_means * 1.5, 1, c(0.5, 0.5)), model = "EEI",
                     replicates = 4L, true_k = 2L),
  ref_eei_n200 = list(make = lpa_case(200, poor_means, 1, c(0.4, 0.35, 0.25)), model = "EEI",
                      replicates = 4L, true_k = 3L),
  ref_eei_n1000 = list(make = lpa_case(1000, two_means, c(1, 2, 0.5, 1), c(0.7, 0.3)),
                       model = "EEI", replicates = 4L, true_k = 2L),
  ref_eei_n4000 = list(make = lpa_case(4000, three_means * 0.5, 1, c(0.4, 0.35, 0.25)),
                       model = "EEI", replicates = 3L, true_k = 3L),
  ref_vvi_n100 = list(make = lpa_case(100, three_means, 1, c(0.4, 0.35, 0.25)), model = "VVI",
                      replicates = 4L, true_k = 3L),
  ref_vvi_n600 = list(make = lpa_case(600, two_means, rbind(c(1, 1, 1, 1), c(0.6, 1.5, 1, 0.8)),
                                      c(0.5, 0.5)), model = "VVI", replicates = 4L, true_k = 2L),
  ref_vvi_n2000 = list(make = lpa_case(2000, three_means * 0.7, rbind(c(1, 1, 1), c(0.7, 1.2, 0.9),
                                                                      c(1.2, 0.8, 1)),
                                       c(0.4, 0.35, 0.25)), model = "VVI", replicates = 3L,
                       true_k = 3L),
  ref_lca_n150 = list(make = lca_case(150, endorse_two, c(0.5, 0.5)), model = "LCA",
                      replicates = 4L, true_k = 2L),
  ref_lca_n3000 = list(make = lca_case(3000, endorse_three, c(0.4, 0.35, 0.25)), model = "LCA",
                       replicates = 3L, true_k = 3L),
  ref_eee_n300 = list(make = lpa_case(300, two_means[, 1:3], c(1, 1.2, 0.8), c(0.6, 0.4),
                                      list(shared_correlation, shared_correlation)),
                      model = "EEE", replicates = 4L, true_k = 2L),
  ref_k1_n2000 = list(make = function(seed) {
    d <- simulab::simulate_correlated(n = 2000, means = c(0, 1, -1), sds = c(1, 2, 0.5),
                                      seed = seed)
    list(data = renamed(d), truth = attr(d, "simulab_tables")$parameters,
         membership = rep("Profile 1", 2000))
  }, model = "EEI", replicates = 4L, true_k = 1L),
  ref_weights_vvi = list(make = function(seed) {
    generated <- lpa_case(300, two_means, rbind(c(1, 1, 1, 1), c(0.7, 1.3, 0.9, 1.1)),
                          c(0.55, 0.45))(seed)
    set.seed(seed)
    raw <- stats::runif(300, 0.2, 3)
    generated$data$w <- raw * 300 / sum(raw)
    generated
  }, model = "VVI", replicates = 4L, true_k = 2L, weights = "w"))

# -- Latent GOLD ---------------------------------------------------------------

latentgold_input <- function(name, spec, seed) {
  categorical <- identical(spec$model, "LCA")
  type <- if (categorical) "nominal" else "continuous"
  vars <- spec$vars
  variances <- if (identical(spec$model, "VVI")) sprintf("   %s | Cluster;", vars) else character()
  covariances <- if (identical(spec$model, "EEE"))
    unlist(lapply(seq_len(length(vars) - 1L), function(i)
      sprintf("   %s <-> %s;", vars[i], vars[(i + 1L):length(vars)]))) else character()
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=2;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", sprintf("      seed=%d sets=200 tolerance=1e-005 iterations=250;", seed),
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    if (isTRUE(spec$missing)) "   missing includeall;" else "   missing excludeall;",
    "   output", "      parameters=first standarderrors;",
    "variables",
    if (!is.null(spec$weights)) sprintf("   caseweight %s;", spec$weights),
    sprintf("   dependent %s;", paste(sprintf("%s %s", vars, type), collapse = ", ")),
    sprintf("   latent Cluster nominal %d:%d;", min(spec$range), max(spec$range)),
    "equations", "   Cluster <- 1;",
    sprintf("   %s <- 1 + Cluster;", vars), variances, covariances, "end model")
}

# -- Mplus ---------------------------------------------------------------------

mplus_input <- function(name, spec, k, columns, seed) {
  vars <- spec$vars
  class_blocks <- if (identical(spec$model, "VVI")) unlist(lapply(seq_len(k), function(j)
    c(sprintf(" %%c#%d%%", j), sprintf(" %s;", paste(vars, collapse = " "))))) else character()
  covariances <- if (identical(spec$model, "EEE"))
    unlist(lapply(seq_len(length(vars) - 1L), function(i)
      sprintf(" %s WITH %s;", vars[i], paste(vars[(i + 1L):length(vars)], collapse = " ")))) else
    character()
  c(sprintf("TITLE: latents VLMR %s, %d against %d classes;", name, k - 1L, k),
    sprintf("DATA: FILE = %s_k%d.dat;", name, k), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", paste(vars, collapse = " ")),
    " IDVARIABLE = id;", " MISSING = ALL(-999);",
    if (identical(spec$model, "LCA")) sprintf(" CATEGORICAL = %s;", paste(vars, collapse = " ")),
    if (!is.null(spec$weights)) sprintf(" WEIGHT = %s;", spec$weights),
    sprintf(" CLASSES = c(%d);", k),
    "ANALYSIS:", " TYPE = MIXTURE;", " ESTIMATOR = MLR;", " STARTS = 200 40;",
    " STITERATIONS = 40;",
    sprintf(" STSEED = %d;", seed),
    " MITERATIONS = 10000;", " LOGCRITERION = 0.00000001;",
    " RLOGCRITERION = 0.000000001;", " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", covariances, class_blocks, "OUTPUT: TECH11;")
}

# -- write ---------------------------------------------------------------------

cases <- list()
invisible(lapply(seq_along(configurations), function(index) {
  config <- names(configurations)[index]
  configuration <- configurations[[config]]
  lapply(seq_len(configuration$replicates), function(replicate) {
    # The last digit of the seed is the replicate (the missing generator reads it).
    seed <- 20261020L + 1000L * index + replicate * 11L
    generated <- configuration$make(seed)
    data <- generated$data
    name <- sprintf("%s_%02d", config, replicate)
    vars <- setdiff(names(data), c("id", configuration$weights))
    spec <- list(model = configuration$model, range = configuration$range,
                 true_k = configuration$true_k, vars = vars, weights = configuration$weights,
                 missing = isTRUE(configuration$missing), n = nrow(data),
                 condition = configuration$condition)
    utils::write.table(data, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
    writeLines(latentgold_input(name, spec, seed), file.path(kit_dir, paste0(name, ".lgs")),
               sep = "\r\n")
    invisible(lapply(setdiff(spec$range, 1L), function(k) {
      stem <- sprintf("%s_k%d", name, k)
      utils::write.table(data, file.path(kit_dir, paste0(stem, ".dat")), sep = " ",
                         quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
      writeLines(mplus_input(name, spec, k, names(data), seed),
                 file.path(kit_dir, paste0(stem, ".inp")))
    }))
    cases[[name]] <<- list(config = config, replicate = replicate, seed = seed, spec = spec,
                           external = TRUE, data = data, truth = generated$truth,
                           membership = generated$membership)
  })
}))
invisible(lapply(seq_along(internal_configurations), function(index) {
  config <- names(internal_configurations)[index]
  configuration <- internal_configurations[[config]]
  lapply(seq_len(configuration$replicates), function(replicate) {
    seed <- 20271020L + 1000L * index + replicate * 11L
    generated <- configuration$make(seed)
    data <- generated$data
    vars <- setdiff(names(data), c("id", configuration$weights))
    spec <- list(model = configuration$model, range = 1:3, true_k = configuration$true_k,
                 vars = vars, weights = configuration$weights, missing = FALSE,
                 n = nrow(data), condition = "reference only")
    cases[[sprintf("%s_%02d", config, replicate)]] <<- list(
      config = config, replicate = replicate, seed = seed, spec = spec, external = FALSE,
      data = data, truth = generated$truth, membership = generated$membership)
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("%d cases (%d run through the programs) in %d configurations written to %s\n",
            length(cases), sum(vapply(cases, `[[`, logical(1), "external")),
            length(configurations) + length(internal_configurations), kit_dir))
