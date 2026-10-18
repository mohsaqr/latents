# Monte Carlo equivalence: simulated datasets with known parameters, written
# for Mplus and Latent GOLD.
#
# Run from the project root:
#   Rscript equivalence/montecarlo/make-kit.R
#   sh equivalence/montecarlo/run.sh
# The outputs land in equivalence/montecarlo/returned/, which
# equivalence/test-montecarlo.R reads offline.
#
# Every configuration is replicated with its own seed. Data are generated
# with simulab (single-level and two-level profiles, latent classes,
# missingness), which stores the generating parameters with each dataset;
# the configurations simulab has no generator for (two-level latent classes,
# class predictors) are generated here, with their parameters stored the same
# way. Separation is moderate, so that entropy and average posterior
# probabilities are informative rather than near one, and every model is
# within the Mplus 9 Demo limits (6 dependent, 2 independent variables).

stopifnot("simulab is needed to build the kit" =
            requireNamespace("simulab", quietly = TRUE))
kit_dir <- file.path("equivalence", "montecarlo", "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
n_replicates <- 15L

renamed <- function(data, prefix = "y") {
  names(data) <- sub("^indicator_", prefix, names(data))
  data
}

three_means <- rbind(c(0, 0, 0, 0), c(1.6, 1.2, -0.8, 0.6), c(-1.2, 1.4, 1.2, -0.8))
three_sds <- rbind(c(1, 1, 1, 1), c(0.8, 1.1, 0.9, 1.2), c(1.2, 0.9, 1, 0.8))

generators <- list(
  lpa_eei = function(seed) {
    d <- simulab::simulate_lpa(n = 500, means = three_means, sds = c(1, 1.2, 0.8, 1),
                               proportions = c(0.45, 0.35, 0.20), seed = seed)
    list(data = renamed(as.data.frame(d)[c("id", paste0("indicator_", 1:4))]),
         truth = attr(d, "simulab_tables")$parameters,
         membership = as.data.frame(d)$profile)
  },
  lpa_vvi = function(seed) {
    d <- simulab::simulate_lpa(n = 500, means = three_means, sds = three_sds,
                               proportions = c(0.45, 0.35, 0.20), seed = seed)
    list(data = renamed(as.data.frame(d)[c("id", paste0("indicator_", 1:4))]),
         truth = attr(d, "simulab_tables")$parameters,
         membership = as.data.frame(d)$profile)
  },
  lpa_vvv = function(seed) {
    correlations <- list(matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3),
                         matrix(c(1, -0.3, 0.1, -0.3, 1, 0.4, 0.1, 0.4, 1), 3))
    d <- simulab::simulate_lpa(n = 500, means = rbind(c(0, 0, 0), c(1.8, 1.2, -1)),
                               sds = rbind(c(1, 1, 1), c(0.8, 1.2, 0.9)),
                               proportions = c(0.6, 0.4), correlations = correlations,
                               seed = seed)
    list(data = renamed(as.data.frame(d)[c("id", paste0("indicator_", 1:3))]),
         truth = attr(d, "simulab_tables")$parameters,
         membership = as.data.frame(d)$profile)
  },
  lca_binary = function(seed) {
    endorse <- rbind(c(0.85, 0.80, 0.75, 0.20, 0.25, 0.30),
                     c(0.20, 0.25, 0.30, 0.80, 0.75, 0.85),
                     c(0.80, 0.75, 0.85, 0.80, 0.80, 0.75))
    probabilities <- array(c(endorse, 1 - endorse), dim = c(3, 6, 2))
    d <- simulab::simulate_lca(n = 800, probabilities = probabilities,
                               proportions = c(0.4, 0.35, 0.25), seed = seed)
    data <- as.data.frame(d)
    items <- grep("^item_", names(data), value = TRUE)
    out <- data.frame(id = data$id,
                      stats::setNames(lapply(data[items], as.integer), sub("item_", "u", items)))
    list(data = out, truth = attr(d, "simulab_tables")$parameters,
         membership = data$latent_class)
  },
  ml_lpa = function(seed) {
    d <- simulab::simulate_ml_lpa(
      clusters = 60, cluster_size = 15, means = rbind(c(0, 0, 0), c(1.5, 1.3, -1.1)),
      profile_probabilities = rbind(c(0.75, 0.25), c(0.25, 0.75)),
      sds = rbind(c(1, 1, 1), c(0.9, 1.1, 1)), cluster_class_proportions = c(0.55, 0.45),
      seed = seed)
    data <- as.data.frame(d)
    list(data = data.frame(id = data$id, g = data$cluster,
                           renamed(data[paste0("indicator_", 1:3)])),
         truth = attr(d, "simulab_tables"), membership = data$profile,
         group_membership = data$cluster_class)
  },
  ml_lca = function(seed) {
    set.seed(seed)
    n_groups <- 60L; group_size <- 15L
    group_class <- sample.int(2L, n_groups, replace = TRUE, prob = c(0.55, 0.45))
    g <- rep(seq_len(n_groups), each = group_size)
    share <- rbind(c(0.75, 0.25), c(0.25, 0.75))
    profile <- vapply(group_class[g], \(h) sample.int(2L, 1L, prob = share[h, ]), integer(1))
    endorse <- rbind(c(0.85, 0.80, 0.75, 0.25, 0.20), c(0.20, 0.30, 0.25, 0.80, 0.85))
    items <- vapply(seq_len(5L), \(j)
      ifelse(stats::runif(length(g)) < endorse[profile, j], 1L, 2L), integer(length(g)))
    colnames(items) <- paste0("u", 1:5)
    truth <- data.frame(class = rep(1:2, 5), item = rep(paste0("u", 1:5), each = 2),
                        probability_first = as.vector(endorse))
    list(data = data.frame(id = seq_along(g), g = g, items), truth = truth,
         membership = profile, group_membership = group_class[g])
  },
  lpa_missing = function(seed) {
    d <- simulab::simulate_lpa(n = 500, means = three_means, sds = three_sds,
                               proportions = c(0.45, 0.35, 0.20), seed = seed)
    data <- renamed(as.data.frame(d)[c("id", paste0("indicator_", 1:4))])
    data[c("y1", "y2", "y3", "y4")] <- simulab::inject_missingness(
      data[c("y1", "y2", "y3", "y4")], mechanism = "MAR", proportion = 0.2,
      variables = c("y2", "y3", "y4"), predictor = "y1", seed = seed + 1L)
    set.seed(seed + 2L)
    data[stats::runif(nrow(data)) < 0.05, c("y1", "y2", "y3", "y4")] <- NA
    list(data = data, truth = attr(d, "simulab_tables")$parameters,
         membership = as.data.frame(d)$profile)
  },
  lpa_covariates = function(seed) {
    set.seed(seed)
    n <- 500L
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    profile <- 1L + (stats::runif(n) < stats::plogis(-0.3 + 0.9 * x1 - 0.6 * x2))
    means <- rbind(c(0, 0, 0), c(1.6, 1.2, -1))
    sds <- rbind(c(1, 1, 1), c(0.8, 1.2, 0.9))
    y <- vapply(1:3, \(j) stats::rnorm(n, means[profile, j], sds[profile, j]), numeric(n))
    colnames(y) <- paste0("y", 1:3)
    list(data = data.frame(id = seq_len(n), x1 = x1, x2 = x2, y),
         truth = list(means = means, sds = sds,
                      logit = c(intercept = -0.3, x1 = 0.9, x2 = -0.6)),
         membership = profile)
  })

specification <- list(
  lpa_eei = list(level = "single", k = 3L, vars = paste0("y", 1:4), model = "EEI"),
  lpa_vvi = list(level = "single", k = 3L, vars = paste0("y", 1:4), model = "VVI"),
  lpa_vvv = list(level = "single", k = 2L, vars = paste0("y", 1:3), model = "VVV"),
  lca_binary = list(level = "single", k = 3L, vars = paste0("u", 1:6), categorical = TRUE),
  ml_lpa = list(level = "two", k = 2L, h = 2L, vars = paste0("y", 1:3), model = "VVI"),
  ml_lca = list(level = "two", k = 2L, h = 2L, vars = paste0("u", 1:5), categorical = TRUE),
  lpa_missing = list(level = "single", k = 3L, vars = paste0("y", 1:4), model = "VVI",
                     missing = TRUE),
  lpa_covariates = list(level = "single", k = 2L, vars = paste0("y", 1:3), model = "VVI",
                        covariates = c("x1", "x2")))

# -- Mplus ----------------------------------------------------------------------

mplus_class_block <- function(spec, class_variable) {
  vars <- paste(spec$vars, collapse = " ")
  unlist(lapply(seq_len(spec$k), function(j) {
    head <- sprintf(" %%%s#%d%%", class_variable, j)
    if (isTRUE(spec$categorical)) {
      return(c(head, sprintf(" [%s];", paste0(spec$vars, "$1", collapse = " "))))
    }
    variances <- if (identical(spec$model, "EEI")) character() else sprintf(" %s;", vars)
    covariances <- if (identical(spec$model, "VVV"))
      unlist(lapply(seq_len(length(spec$vars) - 1L), function(i)
        sprintf(" %s WITH %s;", spec$vars[i],
                paste(spec$vars[(i + 1L):length(spec$vars)], collapse = " ")))) else character()
    c(head, sprintf(" [%s];", vars), variances, covariances)
  }))
}

mplus_input <- function(name, spec, columns) {
  analysis <- c(" ESTIMATOR = ML;", " STARTS = 200 40;", " STITERATIONS = 40;",
                " STSEED = 20261017;", " MITERATIONS = 10000;",
                " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
                " MCONVERGENCE = 0.000001;")
  usevariables <- c(spec$vars, spec$covariates)
  common <- c(sprintf("TITLE: latents Monte Carlo equivalence %s;", name),
              sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
              sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
              sprintf(" USEVARIABLES = %s;", paste(usevariables, collapse = " ")),
              " IDVARIABLE = id;", " MISSING = ALL(-999);",
              if (isTRUE(spec$categorical))
                sprintf(" CATEGORICAL = %s;", paste(spec$vars, collapse = " ")))
  savedata <- c("SAVEDATA:", sprintf(" FILE = %s-posteriors.dat;", name),
                " SAVE = CPROBABILITIES;", " FORMAT = F20.12;",
                sprintf(" RESULTS = %s-results.dat;", name))
  if (identical(spec$level, "single")) {
    return(c(common, sprintf(" CLASSES = c(%d);", spec$k),
             "ANALYSIS:", " TYPE = MIXTURE;", analysis,
             "MODEL:", " %OVERALL%",
             if (length(spec$covariates))
               sprintf(" c ON %s;", paste(spec$covariates, collapse = " ")),
             mplus_class_block(spec, "c"), "OUTPUT: TECH1;", savedata))
  }
  c(common, " CLUSTER = g;", sprintf(" WITHIN = %s;", paste(spec$vars, collapse = " ")),
    sprintf(" CLASSES = cb(%d) cw(%d);", spec$h, spec$k), " BETWEEN = cb;",
    "ANALYSIS:", " TYPE = TWOLEVEL MIXTURE;", analysis,
    "MODEL:", " %WITHIN%", " %OVERALL%", " %BETWEEN%", " %OVERALL%",
    sprintf(" cw#1%s ON cb;", if (spec$k > 2L) sprintf("-cw#%d", spec$k - 1L) else ""),
    "MODEL cw:", " %WITHIN%", mplus_class_block(spec, "cw"), "OUTPUT: TECH1;", savedata)
}

# -- Latent GOLD -----------------------------------------------------------------

latentgold_input <- function(name, spec) {
  two <- identical(spec$level, "two")
  type <- if (isTRUE(spec$categorical)) "nominal" else "continuous"
  dependents <- paste(sprintf("%s %s", spec$vars, type), collapse = ", ")
  variances <- if (isTRUE(spec$categorical) || identical(spec$model, "EEI")) character() else
    sprintf("   %s | Cluster;", spec$vars)
  covariances <- if (identical(spec$model, "VVV"))
    unlist(lapply(seq_len(length(spec$vars) - 1L), function(i)
      sprintf("   %s <-> %s | Cluster;", spec$vars[i], spec$vars[(i + 1L):length(spec$vars)]))) else
    character()
  cluster_equation <- if (two) "   Cluster <- 1 + GClass;" else if (length(spec$covariates))
    sprintf("   Cluster <- 1 + %s;", paste(spec$covariates, collapse = " + ")) else
    "   Cluster <- 1;"
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-010 emtolerance=1e-008 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261017 sets=100 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing includeall;",
    "   output", "      parameters=first profile probmeans=posterior;",
    sprintf("   outfile '%s_posteriors.txt' classification keep id%s;", name,
            if (two) " g" else ""),
    "variables",
    if (two) "   groupid g;",
    if (length(spec$covariates))
      sprintf("   independent %s;", paste(spec$covariates, collapse = ", ")),
    sprintf("   dependent %s;", dependents),
    "   latent",
    if (two) sprintf("      GClass group nominal %d,", spec$h),
    sprintf("      Cluster nominal %d;", spec$k),
    "equations",
    if (two) "   GClass <- 1;",
    cluster_equation,
    sprintf("   %s <- 1 + Cluster;", spec$vars), variances, covariances,
    "end model")
}

# -- write ----------------------------------------------------------------------

cases <- list()
invisible(lapply(seq_along(generators), function(index) {
  config <- names(generators)[index]
  spec <- specification[[config]]
  lapply(seq_len(n_replicates), function(replicate) {
    seed <- 20261017L + 1000L * index + replicate
    generated <- generators[[config]](seed)
    name <- sprintf("%s_%02d", config, replicate)
    data <- generated$data
    columns <- names(data)
    utils::write.table(data, file.path(kit_dir, paste0(name, ".dat")), sep = " ",
                       quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
    utils::write.table(data, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
    writeLines(mplus_input(name, spec, columns), file.path(kit_dir, paste0(name, ".inp")))
    writeLines(latentgold_input(name, spec), file.path(kit_dir, paste0(name, ".lgs")),
               sep = "\r\n")
    cases[[name]] <<- c(list(config = config, replicate = replicate, seed = seed,
                             spec = spec), generated)
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("%d cases in %d configurations written to %s\n", length(cases),
            length(generators), kit_dir))
