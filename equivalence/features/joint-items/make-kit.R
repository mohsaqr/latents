# Joint items (local dependence between categorical items): simulated datasets
# with known parameters, written for Latent GOLD and the Mplus 9 Demo.
#
# Run from the project root:
#   Rscript equivalence/features/joint-items/make-kit.R
#   JOBS=2 sh equivalence/features/run.sh joint-items all
# The outputs land in equivalence/features/joint-items/returned/, which
# equivalence/test-features-joint-items.R reads offline.
#
# The items that are independent given the class are drawn by simulab. It has
# no generator for local dependence, so the second item of each pair is then
# redrawn from a stored class-specific conditional table P(second | first,
# class); the case's truth holds that table and the joint table it implies.
# Latent GOLD fits the pair with `first <-> second | Cluster` (nominal items).
# Mplus cannot fit the association directly: it is given the joined item as one
# NOMINAL variable, which is the same model.
#
# Design (generators.R; 54 cases, 62 Latent GOLD runs, 51 Mplus runs):
#   config       cases k items pair                      n      references
#   bin_small      5   2   5  2 x 2                        80   LG, Mplus
#   bin_mod        5   2   5  2 x 2                       800   LG, Mplus
#   cat3_mod       5   3   6  3 x 3                       800   LG, Mplus
#   bin4_mod       5   2   6  2 x 4                       800   LG, Mplus
#   twopair_mod    5   3   6  two 2 x 2 pairs             800   LG, Mplus
#   k4_large       5   4   7  2 x 2                      5000   LG, Mplus
#   huge           1   3   8  3 x 3                     20000   LG
#   sparse         4   2   5  3 x 3, (1,3), (3,1) never   600   LG, Mplus
#   perfect        3   2   5  p2 == p1, 3 categories      600   LG, Mplus
#   boundary       2   2   5  2 x 3, p2 = 2 never in      700   LG, Mplus
#                             class 1
#   partial        3   2   5  2 x 3; 8% rows miss p2,     800   LG raw (includeall),
#                             5% miss p1, 5% miss u1            LG masked, Mplus
#   indep          3   2   5  2 x 2, independent         1000   LG joint and plain, Mplus
#   strings        2   2   5  3 x 2, string labels        600   LG, Mplus
#   numeric        2   2   5  codes {1, 11} x {0, 10, 1}  600   LG, Mplus
#   weights        2   2   5  2 x 2, weights 1..3         400   LG caseweight and
#                                                               duplicated rows, Mplus
#   multilevel     2   2   5  2 x 2, 2 group classes   60 x 15  LG (two-level)
# Every single-level unweighted case but `huge` is also fitted by poLCA; the
# independent reference (reference.R) polishes every single-level fit and runs
# 140 more datasets (compare.R, `ji_study_design()`).
#
# Rerun: this script, then `JOBS=2 sh equivalence/features/run.sh joint-items`,
# then `Rscript equivalence/run.R features-joint-items` (writes report.csv).

stopifnot("simulab is needed to build the kit" =
            requireNamespace("simulab", quietly = TRUE))
kit_dir <- file.path("equivalence", "features", "joint-items", "kit")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(kit_dir, full.names = TRUE))

source(file.path("equivalence", "features", "joint-items", "generators.R"))

# -- writers ---------------------------------------------------------------------

#' Latent GOLD syntax: nominal items, one `a <-> b | Cluster` per pair
latentgold_input <- function(name, items, k, pairs, h = NULL, weight = NULL,
                             associate = TRUE) {
  two <- !is.null(h)
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s-lg.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=1;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261020 sets=150 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing includeall;",
    "   output", "      parameters=first standarderrors profile estimatedvalues=model;",
    sprintf("   outfile '%s_posteriors.txt' classification keep id%s;", name,
            if (two) " g" else ""),
    "variables",
    if (two) "   groupid g;",
    if (!is.null(weight)) sprintf("   caseweight %s;", weight),
    sprintf("   dependent %s;", paste(sprintf("%s nominal", items), collapse = ", ")),
    "   latent",
    if (two) sprintf("      GClass group nominal %d,", h),
    sprintf("      Cluster nominal %d;", k),
    "equations",
    if (two) "   GClass <- 1;",
    if (two) "   Cluster <- 1 + GClass;" else "   Cluster <- 1;",
    sprintf("   %s <- 1 + Cluster;", items),
    if (associate) vapply(strsplit(pairs, ":", fixed = TRUE),
                          \(p) sprintf("   %s <-> %s | Cluster;", p[1], p[2]), ""),
    "end model")
}

#' Mplus input: the joined items as NOMINAL, the others as binary CATEGORICAL
mplus_input <- function(name, columns, nominal, binary, k, weight = NULL) {
  vars <- c(nominal, binary)
  class_block <- unlist(lapply(seq_len(k), function(j) {
    c(sprintf(" %%c#%d%%", j),
      if (length(binary)) sprintf(" [%s];", paste(sprintf("%s$1", names(binary)), collapse = " ")),
      vapply(names(nominal), \(v) sprintf(" [%s];", paste(
        sprintf("%s#%d", v, seq_len(nominal[[v]] - 1L)), collapse = " ")), ""))
  }))
  c(sprintf("TITLE: latents joint items %s;", name),
    sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", paste(names(vars), collapse = " ")),
    " IDVARIABLE = id;", " MISSING = ALL(-999);",
    if (length(binary)) sprintf(" CATEGORICAL = %s;", paste(names(binary), collapse = " ")),
    sprintf(" NOMINAL = %s;", paste(names(nominal), collapse = " ")),
    if (!is.null(weight)) sprintf(" WEIGHT = %s;", weight),
    sprintf(" CLASSES = c(%d);", k),
    "ANALYSIS:", " TYPE = MIXTURE;",
    if (is.null(weight)) " ESTIMATOR = ML;" else " ESTIMATOR = MLR;",
    " STARTS = 400 80;", " STITERATIONS = 40;",
    " STSEED = 20261020;", " MITERATIONS = 10000;",
    " LOGCRITERION = 0.00000001;", " RLOGCRITERION = 0.000000001;",
    " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", class_block, "OUTPUT: TECH1;",
    "SAVEDATA:", sprintf(" FILE = %s-posteriors.dat;", name),
    " SAVE = CPROBABILITIES;", " FORMAT = F20.12;",
    sprintf(" RESULTS = %s-results.dat;", name))
}

write_lg_data <- function(frame, name) {
  utils::write.table(frame, file.path(kit_dir, paste0(name, "-lg.dat")), sep = "\t",
                     quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
}

write_lg <- function(name, ...) {
  writeLines(latentgold_input(name, ...), file.path(kit_dir, paste0(name, ".lgs")),
             sep = "\r\n")
}

#' Integer codes of a joined pair, 1..K over the observed combinations,
#' missing when either item is
joint_codes <- function(first, second) {
  label <- ifelse(is.na(first) | is.na(second), NA_character_, paste(first, second, sep = ":"))
  levels <- unique(label[order(first, second)])
  match(label, levels[!is.na(levels)])
}

# -- write ------------------------------------------------------------------------

cases <- list()
invisible(lapply(seq_along(generate), function(index) {
  config <- names(generate)[index]
  setting <- specification[[config]]
  lapply(seq_len(setting$replicates), function(replicate) {
    seed <- 20261020L + 1000L * index + replicate
    name <- sprintf("%s_%02d", config, replicate)
    drawn <- generate[[config]](seed)
    items <- drawn$items
    # Codes become their labels where a case has labels.
    if (!is.null(drawn$labels)) {
      items[names(drawn$labels)] <- lapply(names(drawn$labels),
                                           \(v) drawn$labels[[v]][items[[v]]])
    }
    data <- data.frame(Filter(Negate(is.null), list(id = seq_len(nrow(items)),
                                                     g = drawn$group, w = drawn$weights)),
                       items)
    vars <- names(items)
    pairs <- setting$pairs %||% "p1:p2"
    split_pairs <- strsplit(pairs, ":", fixed = TRUE)
    h <- setting$h
    weight <- if (!is.null(drawn$weights)) "w"
    lg_runs <- name
    write_lg_data(data, name)
    write_lg(name, items = vars, k = setting$k, pairs = pairs, h = h, weight = weight)
    if (identical(config, "partial")) {
      # The pair made missing wherever one item is: the convention of latents.
      masked <- data
      either <- is.na(masked$p1) | is.na(masked$p2)
      masked$p1[either] <- NA
      masked$p2[either] <- NA
      variant <- sub("partial", "partialmask", name)
      write_lg_data(masked, variant)
      write_lg(variant, items = vars, k = setting$k, pairs = pairs)
      lg_runs <- c(lg_runs, mask = variant)
    }
    if (identical(config, "indep")) {
      variant <- sub("indep", "indepplain", name)
      write_lg_data(data, variant)
      write_lg(variant, items = vars, k = setting$k, pairs = pairs, associate = FALSE)
      lg_runs <- c(lg_runs, plain = variant)
    }
    if (identical(config, "weights")) {
      variant <- sub("weights", "weightsdup", name)
      duplicated_rows <- data[rep(seq_len(nrow(data)), data$w), setdiff(names(data), "w")]
      duplicated_rows$id <- seq_len(nrow(duplicated_rows))
      write_lg_data(duplicated_rows, variant)
      write_lg(variant, items = vars, k = setting$k, pairs = pairs)
      lg_runs <- c(lg_runs, duplicated = variant)
    }
    names(lg_runs)[1L] <- "main"
    mplus_coding <- NULL
    if (isTRUE(setting$mplus)) {
      paired <- unlist(split_pairs)
      nominal_frame <- stats::setNames(lapply(split_pairs, \(p) joint_codes(
        drawn$items[[p[1]]], drawn$items[[p[2]]])), vapply(split_pairs, paste, "",
                                                            collapse = ""))
      binary_frame <- drawn$items[setdiff(vars, paired)]
      mplus_frame <- data.frame(Filter(Negate(is.null), list(id = data$id, w = drawn$weights)),
                                nominal_frame, binary_frame)
      nominal <- vapply(nominal_frame, \(v) max(v, na.rm = TRUE), integer(1))
      binary <- stats::setNames(rep(2L, ncol(binary_frame)), names(binary_frame))
      utils::write.table(mplus_frame, file.path(kit_dir, paste0(name, ".dat")), sep = " ",
                         quote = FALSE, row.names = FALSE, col.names = FALSE, na = "-999")
      writeLines(mplus_input(name, names(mplus_frame), nominal, binary, setting$k, weight),
                 file.path(kit_dir, paste0(name, ".inp")))
      # The joint labels each Mplus code stands for, as latents writes them.
      mplus_coding <- lapply(stats::setNames(seq_along(split_pairs), names(nominal_frame)),
                             function(i) {
        p <- split_pairs[[i]]
        code <- nominal_frame[[i]]
        label <- paste(data[[p[1]]], data[[p[2]]], sep = ":")
        unique(data.frame(code = code, label = label)[!is.na(code), ])
      })
    }
    cases[[name]] <<- list(
      name = name, config = config, replicate = replicate, seed = seed, data = data,
      spec = list(vars = vars, pairs = pairs, k = setting$k, h = h,
                  id = if (!is.null(h)) "g", weights = weight,
                  missing = isTRUE(setting$missing), latentgold = lg_runs,
                  mplus = isTRUE(setting$mplus), mplus_coding = mplus_coding),
      truth = drawn$truth, membership = drawn$membership)
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("%d cases (%d Latent GOLD runs, %d Mplus runs) written to %s\n", length(cases),
            length(list.files(kit_dir, "[.]lgs$")), length(list.files(kit_dir, "[.]inp$")),
            kit_dir))
