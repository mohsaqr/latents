# Builds the Latent GOLD kit for the group-class families. Run from the
# project root: Rscript equivalence/latentgold-families/make-kit.R

kit_dir <- file.path("equivalence", "latentgold-families", "kit")
wine_dir <- file.path(Sys.getenv("HOME"), ".wine-latentgold", "drive_c", "lgfamilies")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wine_dir, showWarnings = FALSE, recursive = TRUE)

cases <- list(
  g01 = list(family = "additive", between = "varying", H = 2L, d = 1L, seed = 101L,
             means = rbind(-1, 1), T = rbind(0.3, 0.6), W = rbind(0.8, 0.8)),
  g02 = list(family = "additive", between = "equal", H = 2L, d = 1L, seed = 102L,
             means = rbind(-1, 1), T = rbind(0.4, 0.4), W = rbind(0.8, 0.8)),
  g03 = list(family = "dispersion", between = "equal", H = 2L, d = 1L, seed = 103L,
             means = rbind(0.5, 0.5), T = rbind(0.4, 0.4), W = rbind(0.4, 2.0)),
  g04 = list(family = "additive_dispersion", between = "varying", H = 2L, d = 1L,
             seed = 104L, means = rbind(-1, 1), T = rbind(0.3, 0.6),
             W = rbind(0.5, 1.5)),
  g05 = list(family = "additive_dispersion", between = "equal", H = 2L, d = 1L,
             seed = 105L, means = rbind(-1, 1), T = rbind(0.4, 0.4),
             W = rbind(0.5, 1.5)),
  g06 = list(family = "additive", between = "varying", H = 3L, d = 1L, seed = 106L,
             means = rbind(-2, 0, 2), T = rbind(0.3, 0.4, 0.5), W = rbind(1, 1, 1)),
  g07 = list(family = "additive", between = "varying", H = 2L, d = 2L, seed = 107L,
             means = rbind(c(-1, 0.5), c(1, -0.5)), T = rbind(c(0.3, 0.4), c(0.5, 0.2)),
             W = rbind(c(0.8, 1.1), c(0.8, 1.1))),
  g08 = list(family = "dispersion", between = "equal", H = 2L, d = 2L, seed = 108L,
             means = rbind(c(0.5, -0.5), c(0.5, -0.5)), T = rbind(c(0.4, 0.3), c(0.4, 0.3)),
             W = rbind(c(0.4, 0.5), c(1.8, 2.2))))

simulate <- function(case) {
  set.seed(case$seed)
  sizes <- rep_len(c(4L, 6L, 8L, 10L), 150L)
  classes <- sample.int(case$H, length(sizes), replace = TRUE)
  rows <- lapply(seq_along(sizes), function(j) {
    h <- classes[j]
    intercept <- stats::rnorm(case$d, case$means[h, ], sqrt(case$T[h, ]))
    ratings <- sweep(matrix(stats::rnorm(sizes[j] * case$d), sizes[j], case$d) %*%
                       diag(sqrt(case$W[h, ]), case$d), 2L, intercept, "+")
    colnames(ratings) <- paste0("y", seq_len(case$d))
    data.frame(g = j, group_class = h, ratings)
  })
  data <- do.call(rbind, rows)
  cbind(id = seq_len(nrow(data)), data)
}

syntax <- function(name, case) {
  vars <- paste0("y", seq_len(case$d))
  factors <- paste0("F", seq_len(case$d))
  class_means <- case$family != "dispersion"
  class_within <- case$family != "additive"
  class_between <- identical(case$between, "varying")
  nodes <- if (case$d == 1L) 50L else 15L
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-010 emtolerance=1e-008 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261001 sets=20 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    sprintf("   quadrature nodes=%d;", nodes), "   missing excludeall;",
    "   output", "      parameters=first standarderrors profile estimatedvalues=model;",
    sprintf("   outfile '%s_posteriors.txt' classification keep id g;", name),
    "variables", "   groupid g;",
    sprintf("   dependent %s;", paste(paste(vars, "continuous"), collapse = ", ")),
    "   latent",
    sprintf("      GClass group nominal %d,", case$H),
    paste0("      ", paste(paste(factors, "group continuous"), collapse = ",\n      "), ";"),
    "equations", "   GClass <- 1;",
    sprintf("   %s%s;", factors, if (class_between) " | GClass" else ""),
    sprintf("   %s <- 1%s + (1) %s;", vars, if (class_means) " + GClass" else "", factors),
    sprintf("   %s%s;", vars, if (class_within) " | GClass" else ""),
    "end model")
}

invisible(lapply(names(cases), function(name) {
  data <- simulate(cases[[name]])
  utils::write.table(data, file.path(kit_dir, paste0(name, ".dat")), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  writeLines(syntax(name, cases[[name]]), file.path(kit_dir, paste0(name, ".lgs")))
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
file.copy(list.files(kit_dir, full.names = TRUE), wine_dir, overwrite = TRUE)
cat("kit written:", length(cases), "cases\n")
