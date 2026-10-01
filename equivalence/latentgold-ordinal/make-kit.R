# Latent GOLD 6.1 kit for ordinal (adjacent-category logit, uniform scores)
# and Poisson count indicators in a latent class model. Run from the package
# root: Rscript equivalence/latentgold-ordinal/make-kit.R
kit_dir <- file.path("equivalence", "latentgold-ordinal", "kit")
wine_dir <- file.path(Sys.getenv("HOME"), ".wine-latentgold", "drive_c", "lgordinal")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wine_dir, showWarnings = FALSE, recursive = TRUE)

# Each case: classes, class sizes, and per indicator its type and per-class
# truth (ordinal: intercepts a_2..a_K and locations; count: means; continuous:
# means and sd).
cases <- list(
  o01 = list(seed = 201L, n = 800L, sizes = c(0.6, 0.4), indicators = list(
    y = list(type = "continuous", mean = c(0, 1.5), sd = c(1, 0.8)),
    o = list(type = "ordinal", intercepts = c(0.5, 0, -1), locations = c(1.2, 0)),
    k = list(type = "count", mean = c(1, 5)))),
  o02 = list(seed = 202L, n = 900L, sizes = c(0.4, 0.35, 0.25), indicators = list(
    o1 = list(type = "ordinal", intercepts = c(0.3, -0.4), locations = c(1.5, 0.6, 0)),
    o2 = list(type = "ordinal", intercepts = c(0.8, 0.2, -0.3, -1.2),
              locations = c(-1, 0.8, 0)),
    k = list(type = "count", mean = c(0.5, 2, 6)))),
  o03 = list(seed = 203L, n = 600L, sizes = c(0.5, 0.5), indicators = list(
    k1 = list(type = "count", mean = c(0.8, 4)),
    k2 = list(type = "count", mean = c(3, 1)))),
  # Ordinal only. (Two ordinal indicators with three classes would be
  # saturated, df = 0, with no unique maximum; this case is identified.)
  o04 = list(seed = 204L, n = 700L, sizes = c(0.55, 0.45), indicators = list(
    o1 = list(type = "ordinal", intercepts = c(0.2, 0.1, -0.6), locations = c(1.4, 0)),
    o2 = list(type = "ordinal", intercepts = c(0.5, -0.5), locations = c(-1.2, 0)),
    o3 = list(type = "ordinal", intercepts = c(0.6, 0.3, -0.2, -0.9),
              locations = c(0.9, 0)))),
  # Two-level: 120 groups of 6; the group class shifts the profile shares.
  o05 = list(seed = 205L, n = 720L, sizes = c(0.5, 0.5), group_size = 6L,
             group_classes = c(0.5, 0.5), shares = rbind(c(0.8, 0.2), c(0.25, 0.75)),
             indicators = list(
    y = list(type = "continuous", mean = c(0, 1.5), sd = c(1, 1)),
    o = list(type = "ordinal", intercepts = c(0.4, 0, -0.8), locations = c(1.3, 0)),
    k = list(type = "count", mean = c(1.5, 4)))),
  # Negative binomial (Latent GOLD's `poisson overdispersed`): dispersion by
  # class, then shared and mixed with an ordinal indicator.
  o06 = list(seed = 206L, n = 1200L, sizes = c(0.6, 0.4), indicators = list(
    y = list(type = "continuous", mean = c(0, 1.5), sd = c(1, 1)),
    k = list(type = "negbin", mean = c(2, 8), dispersion = c(0.5, 0.2)))),
  o07 = list(seed = 207L, n = 1000L, sizes = c(0.5, 0.5), dispersion = "equal",
             indicators = list(
    o = list(type = "ordinal", intercepts = c(0.3, -0.2, -0.9), locations = c(1.1, 0)),
    k = list(type = "negbin", mean = c(1.5, 5), dispersion = c(0.4, 0.4)))))

ordinal_probabilities <- function(intercepts, locations) {
  logits <- outer(locations, seq_along(c(0, intercepts)) - 1) +
    matrix(c(0, intercepts), length(locations), length(intercepts) + 1L, byrow = TRUE)
  exp(logits) / rowSums(exp(logits))
}

simulate <- function(case) {
  set.seed(case$seed)
  class <- sample.int(length(case$sizes), case$n, replace = TRUE, prob = case$sizes)
  group <- NULL
  if (!is.null(case$group_size)) {
    group <- rep(seq_len(case$n / case$group_size), each = case$group_size)
    group_class <- sample.int(length(case$group_classes), max(group), replace = TRUE,
                              prob = case$group_classes)
    class <- vapply(group_class[group], function(h) {
      sample.int(ncol(case$shares), 1L, prob = case$shares[h, ])
    }, integer(1))
  }
  columns <- lapply(case$indicators, function(indicator) {
    switch(indicator$type,
      continuous = stats::rnorm(case$n, indicator$mean[class], indicator$sd[class]),
      count = stats::rpois(case$n, indicator$mean[class]),
      negbin = stats::rnbinom(case$n, size = 1 / indicator$dispersion[class],
                              mu = indicator$mean[class]),
      ordinal = {
        probabilities <- ordinal_probabilities(indicator$intercepts, indicator$locations)
        vapply(class, function(h) sample.int(ncol(probabilities), 1L,
                                             prob = probabilities[h, ]), integer(1))
      })
  })
  frame <- data.frame(id = seq_len(case$n), columns)
  if (!is.null(group)) frame$g <- group
  frame
}

syntax <- function(name, case) {
  vars <- names(case$indicators)
  types <- vapply(case$indicators, `[[`, character(1), "type")
  lg_type <- c(continuous = "continuous", ordinal = "ordinal", count = "poisson",
               negbin = "poisson overdispersed")[types]
  continuous <- vars[types == "continuous"]
  # Overdispersed counts take a variance (dispersion) equation like continuous
  # indicators: by class, or shared when the case says so.
  overdispersed <- vars[types == "negbin"]
  shared <- identical(case$dispersion, "equal")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261001 sets=50 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;",
    "   output", "      parameters=first standarderrors profile;",
    sprintf("   outfile '%s_posteriors.txt' classification keep id;", name),
    "variables",
    if (!is.null(case$group_size)) "   groupid g;",
    sprintf("   dependent %s;", paste(paste(vars, lg_type), collapse = ", ")),
    if (is.null(case$group_size)) sprintf("   latent Cluster nominal %d;", length(case$sizes)) else
      c("   latent", sprintf("      GClass group nominal %d,", length(case$group_classes)),
        sprintf("      Cluster nominal %d;", length(case$sizes))),
    "equations",
    if (is.null(case$group_size)) "   Cluster <- 1;" else
      c("   GClass <- 1;", "   Cluster <- 1 + GClass;"),
    sprintf("   %s <- 1 + Cluster;", vars),
    if (length(continuous) > 0L) sprintf("   %s | Cluster;", continuous),
    if (length(overdispersed) > 0L)
      sprintf(if (shared) "   %s;" else "   %s | Cluster;", overdispersed),
    "end model")
}

invisible(lapply(names(cases), function(name) {
  utils::write.table(simulate(cases[[name]]), file.path(kit_dir, paste0(name, ".dat")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(syntax(name, cases[[name]]), file.path(kit_dir, paste0(name, ".lgs")))
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
file.copy(list.files(kit_dir, full.names = TRUE), wine_dir, overwrite = TRUE)
cat("kit written:", length(cases), "cases\n")
