# Datasets with missing indicators, written for Mplus and Latent GOLD.
#
# Run from the project root:
#   Rscript equivalence/missing-patterns/make-kit.R
# then run both programs with equivalence/missing-patterns/run.sh. Their
# outputs land in returned/, which test-missing-patterns.R reads.
#
# Each case is one missingness pattern on simulated data with moderately
# separated profiles, so entropy and AvePP are informative rather than near 1.
# Missing values are written as -999 for Mplus and "." for Latent GOLD.

kit_dir <- file.path("equivalence", "missing-patterns", "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
vars <- c("y1", "y2", "y3")

# -- data ---------------------------------------------------------------------

#' Single-level data: three profiles, moderately separated
single_level_data <- function(seed, n = 600L) {
  set.seed(seed)
  means <- rbind(c(0, 0, 0), c(1.6, 1.2, -0.8), c(-1.2, 1.6, 1.4))
  sds <- rbind(c(1, 1, 1), c(0.8, 1.1, 0.9), c(1.2, 0.9, 1))
  profile <- sample.int(3L, n, replace = TRUE, prob = c(0.45, 0.35, 0.20))
  values <- vapply(seq_along(vars), \(j)
    stats::rnorm(n, means[profile, j], sds[profile, j]), numeric(n))
  colnames(values) <- vars
  data.frame(id = seq_len(n), profile = profile, values)
}

#' Two-level data: two profiles whose shares differ by group class
two_level_data <- function(seed, n_groups = 60L, group_size = 15L) {
  set.seed(seed)
  group_class <- sample.int(2L, n_groups, replace = TRUE, prob = c(0.55, 0.45))
  g <- rep(seq_len(n_groups), each = group_size)
  share <- c(0.75, 0.25)[group_class[g]]
  profile <- ifelse(stats::runif(length(g)) < share, 1L, 2L)
  means <- rbind(c(0, 0, 0), c(1.5, 1.3, -1.1))
  values <- vapply(seq_along(vars), \(j)
    stats::rnorm(length(g), means[profile, j], 1), numeric(length(g)))
  colnames(values) <- vars
  data.frame(id = seq_along(g), g = g, profile = profile,
             group_class = group_class[g], values)
}

#' Blank a share of indicator cells completely at random
#'
#' A row the draw empties keeps its `y1`, so cell blanking alone never makes a
#' row with nothing observed; those are made deliberately by `blank_rows()`.
blank_cells <- function(data, share) {
  original <- as.matrix(data[vars])
  cells <- original
  cells[stats::runif(length(cells)) < share] <- NA
  emptied <- rowSums(!is.na(cells)) == 0L
  cells[emptied, "y1"] <- original[emptied, "y1"]
  data[vars] <- as.data.frame(cells)
  data
}

#' Blank every indicator of a share of rows
blank_rows <- function(data, share) {
  rows <- stats::runif(nrow(data)) < share
  data[rows, vars] <- NA
  data
}

cases <- list(
  m1 = list(
    description = "Single level, 3 profiles; 15% of cells missing completely at random.",
    level = "single", n_profiles = 3L,
    data = local({ d <- single_level_data(101); set.seed(1101); blank_cells(d, 0.15) })),
  m2 = list(
    description = "m1 plus 10% of rows with every indicator missing.",
    level = "single", n_profiles = 3L,
    data = local({ d <- single_level_data(101); set.seed(1101); d <- blank_cells(d, 0.15)
                   set.seed(2101); blank_rows(d, 0.10) })),
  m3 = list(
    description = "Single level, 3 profiles; y2 and y3 missing more often as y1 rises (MAR on y1).",
    level = "single", n_profiles = 3L,
    data = local({
      d <- single_level_data(103); set.seed(1103)
      p <- stats::plogis(-1.4 + 1.1 * d$y1)
      d$y2[stats::runif(nrow(d)) < p] <- NA
      d$y3[stats::runif(nrow(d)) < p] <- NA
      d })),
  m4 = list(
    description = "Single level, 3 profiles; y3 missing for 70% of rows.",
    level = "single", n_profiles = 3L,
    data = local({ d <- single_level_data(104); set.seed(1104)
                   d$y3[stats::runif(nrow(d)) < 0.70] <- NA; d })),
  m5 = list(
    description = paste("Two level, 2 profiles by 2 group classes; 10% of cells MCAR,",
                        "8% of rows empty, and groups 7 and 31 entirely empty."),
    level = "two", n_profiles = 2L, n_group_classes = 2L,
    data = local({ d <- two_level_data(105); set.seed(1105); d <- blank_cells(d, 0.10)
                   set.seed(2105); d <- blank_rows(d, 0.08)
                   d[d$g %in% c(7L, 31L), vars] <- NA; d }))
)

# -- Mplus --------------------------------------------------------------------

mplus_class_block <- function(class_name, k) {
  unlist(lapply(seq_len(k), \(j) c(sprintf(" %%%s#%d%%", class_name, j),
                                   " [y1 y2 y3];", " y1 y2 y3;")))
}

mplus_input <- function(name, case) {
  common_analysis <- c(" ESTIMATOR = ML;", " STARTS = 400 80;", " STITERATIONS = 50;",
                       " STSEED = 20261017;", " MITERATIONS = 10000;",
                       " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
                       " MCONVERGENCE = 0.000001;")
  savedata <- c("SAVEDATA:", sprintf(" FILE = %s-posteriors.dat;", name),
                " SAVE = CPROBABILITIES;", " FORMAT = F20.12;",
                sprintf(" RESULTS = %s-results.dat;", name))
  if (identical(case$level, "single")) {
    return(c(sprintf("TITLE: latents missing-pattern case %s;", name),
             sprintf("DATA: FILE = %s.dat;", name),
             "VARIABLE:", " NAMES = id y1 y2 y3;", " USEVARIABLES = y1 y2 y3;",
             " IDVARIABLE = id;", " MISSING = ALL(-999);",
             sprintf(" CLASSES = c(%d);", case$n_profiles),
             "ANALYSIS:", " TYPE = MIXTURE;", common_analysis,
             "MODEL:", " %OVERALL%", mplus_class_block("c", case$n_profiles),
             "OUTPUT: TECH1;", savedata))
  }
  k <- case$n_profiles
  c(sprintf("TITLE: latents missing-pattern case %s;", name),
    sprintf("DATA: FILE = %s.dat;", name),
    "VARIABLE:", " NAMES = id g y1 y2 y3;", " USEVARIABLES = y1 y2 y3;",
    " CLUSTER = g;", " IDVARIABLE = id;", " MISSING = ALL(-999);",
    " WITHIN = y1 y2 y3;",
    sprintf(" CLASSES = cb(%d) cw(%d);", case$n_group_classes, k), " BETWEEN = cb;",
    "ANALYSIS:", " TYPE = TWOLEVEL MIXTURE;", common_analysis,
    "MODEL:", " %WITHIN%", " %OVERALL%", " %BETWEEN%", " %OVERALL%",
    sprintf(" cw#1%s ON cb;", if (k > 2L) sprintf("-cw#%d", k - 1L) else ""),
    "MODEL cw:", " %WITHIN%", mplus_class_block("cw", k),
    "OUTPUT: TECH1;", savedata)
}

# -- Latent GOLD --------------------------------------------------------------

latentgold_input <- function(name, case) {
  two <- identical(case$level, "two")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-010 emtolerance=1e-008 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261017 sets=200 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing includeall;",
    "   output", "      parameters=first profile probmeans=posterior;",
    sprintf("   outfile '%s_posteriors.txt' classification keep id%s;", name,
            if (two) " g" else ""),
    "variables",
    if (two) "   groupid g;",
    "   dependent y1 continuous, y2 continuous, y3 continuous;",
    "   latent",
    if (two) sprintf("      GClass group nominal %d,", case$n_group_classes),
    sprintf("      Cluster nominal %d;", case$n_profiles),
    "equations",
    if (two) c("   GClass <- 1;", "   Cluster <- 1 + GClass;") else "   Cluster <- 1;",
    "   y1 <- 1 + Cluster;", "   y2 <- 1 + Cluster;", "   y3 <- 1 + Cluster;",
    "   y1 | Cluster;", "   y2 | Cluster;", "   y3 | Cluster;",
    "end model")
}

# -- write ----------------------------------------------------------------------

invisible(lapply(names(cases), function(name) {
  case <- cases[[name]]
  columns <- if (identical(case$level, "two")) c("id", "g", vars) else c("id", vars)
  data <- case$data[columns]
  utils::write.table(data, file.path(kit_dir, paste0(name, ".dat")), sep = " ",
                     quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
  utils::write.table(data, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                     quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
  writeLines(mplus_input(name, case), file.path(kit_dir, paste0(name, ".inp")))
  writeLines(latentgold_input(name, case), file.path(kit_dir, paste0(name, ".lgs")),
             sep = "\r\n")
}))
saveRDS(lapply(cases, \(case) case[setdiff(names(case), character())]),
        file.path(kit_dir, "cases.rds"))

summary_table <- do.call(rbind, lapply(names(cases), function(name) {
  d <- cases[[name]]$data
  observed <- rowSums(!is.na(d[vars]))
  data.frame(case = name, rows = nrow(d), cells_missing = sum(is.na(d[vars])),
             rows_some_missing = sum(observed < length(vars) & observed > 0),
             rows_all_missing = sum(observed == 0),
             groups = if (is.null(d$g)) NA_integer_ else length(unique(d$g)),
             groups_all_missing = if (is.null(d$g)) NA_integer_ else
               sum(tapply(observed, d$g, sum) == 0))
}))
print(summary_table)
