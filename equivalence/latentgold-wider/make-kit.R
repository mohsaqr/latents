# Builds the Latent GOLD kit for the outputs no Latent GOLD comparison covered
# before: mixture regression (estimates, model-based and robust standard
# errors; Gaussian, Poisson and logistic), and sampling weights in lpa(),
# lta() and two-level multilpa(). Run from the project root:
#   Rscript equivalence/latentgold-wider/make-kit.R
# then equivalence/latentgold-wider/run-kit.sh, then compare.R.
#
# Writes kit/<case>.dat and kit/<case>.lgs, the same files to the Wine drive,
# and kit/cases.rds with each case's simulated data and settings. The latents
# fits are made in compare.R, from those data, so the comparison always reads
# the package as it is now.

kit_dir <- file.path("equivalence", "latentgold-wider", "kit")
wine_dir <- file.path(Sys.getenv("HOME"), ".wine-latentgold", "drive_c", "lgwider")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wine_dir, showWarnings = FALSE, recursive = TRUE)

# Weighted cases ask for Latent GOLD's robust (sandwich) standard errors: with a
# sampling weight it otherwise reports model-based ones, which latents
# refuses for weighted fits.
#
# Sampling weights in (0.4, 2.5), unrelated to class, as a survey design would
# give: the weighted and unweighted fits differ, and both are identified.
.draw_weights <- function(n) round(stats::runif(n, 0.4, 2.5), 4)

simulate_regression <- function(seed, family, n = 800L) {
  set.seed(seed)
  z <- round(stats::rnorm(n), 4)
  x <- round(stats::rnorm(n), 4)
  share <- stats::plogis(-0.4 + 0.8 * z)
  class <- 1L + (stats::runif(n) < share)
  intercept <- c(1, -0.5)[class]
  slope <- c(0.5, 1.4)[class]
  eta <- intercept + slope * x
  y <- switch(family,
    gaussian = round(eta + stats::rnorm(n, sd = c(0.6, 1.1)[class]), 4),
    poisson = stats::rpois(n, exp(0.3 * eta + c(0.2, 0.9)[class])),
    binomial = as.integer(stats::runif(n) < stats::plogis(1.2 * eta)))
  data.frame(id = seq_len(n), y = y, x = x, z = z, w = .draw_weights(n))
}

# A logistic mixture needs several outcomes per class assignment (latents
# refuses one binary outcome per unit, which is not identified): five binary
# responses per person, the class and z at the person level.
simulate_repeated_binary <- function(seed, persons = 300L, repeats = 5L) {
  set.seed(seed)
  z <- round(stats::rnorm(persons), 4)
  class <- 1L + (stats::runif(persons) < stats::plogis(-0.3 + 0.9 * z))
  long <- data.frame(id = rep(seq_len(persons), each = repeats),
                     x = round(stats::rnorm(persons * repeats), 4),
                     z = rep(z, each = repeats), class = rep(class, each = repeats))
  eta <- c(1, -0.8)[long$class] + c(0.6, 1.5)[long$class] * long$x
  long$y <- as.integer(stats::runif(nrow(long)) < stats::plogis(eta))
  long$w <- 1
  long[c("id", "y", "x", "z", "w")]
}

simulate_profiles <- function(seed, n = 900L) {
  set.seed(seed)
  class <- sample.int(3L, n, replace = TRUE, prob = c(0.45, 0.35, 0.2))
  means <- rbind(c(-1, 0, 1), c(1, 1, -0.5), c(0, -1.5, 0.5))
  sds <- rbind(c(0.7, 0.8, 0.6), c(0.9, 0.6, 0.8), c(0.6, 0.7, 1.0))
  values <- means[class, ] + matrix(stats::rnorm(n * 3L), n, 3L) * sds[class, ]
  data.frame(id = seq_len(n), y1 = round(values[, 1L], 4), y2 = round(values[, 2L], 4),
             y3 = round(values[, 3L], 4), w = .draw_weights(n))
}

simulate_transitions <- function(seed, persons = 400L, occasions = 4L) {
  set.seed(seed)
  transition <- rbind(c(0.85, 0.15), c(0.25, 0.75))
  state <- matrix(0L, persons, occasions)
  state[, 1L] <- 1L + (stats::runif(persons) < 0.4)
  invisible(lapply(seq_len(occasions - 1L), function(t) {
    state[, t + 1L] <<- 1L + (stats::runif(persons) < transition[state[, t], 2L])
  }))
  means <- rbind(c(-0.8, 0.5), c(0.9, -0.6))
  long <- data.frame(id = rep(seq_len(persons), each = occasions),
                     time = rep(seq_len(occasions), times = persons),
                     state = as.vector(t(state)))
  long$y1 <- round(means[long$state, 1L] + stats::rnorm(nrow(long), sd = 0.7), 4)
  long$y2 <- round(means[long$state, 2L] + stats::rnorm(nrow(long), sd = 0.8), 4)
  long$w <- rep(.draw_weights(persons), each = occasions)
  long$state <- NULL
  long
}

simulate_two_level <- function(seed, groups = 120L) {
  set.seed(seed)
  sizes <- rep_len(c(6L, 8L, 10L, 12L), groups)
  group_class <- 1L + (stats::runif(groups) < 0.45)
  mix <- rbind(c(0.8, 0.2), c(0.3, 0.7))
  rows <- lapply(seq_len(groups), function(j) {
    profile <- 1L + (stats::runif(sizes[j]) < mix[group_class[j], 2L])
    data.frame(g = j,
               y1 = round(c(-0.7, 0.8)[profile] + stats::rnorm(sizes[j], sd = 0.7), 4),
               y2 = round(c(0.6, -0.5)[profile] + stats::rnorm(sizes[j], sd = 0.8), 4))
  })
  data <- do.call(rbind, rows)
  data$w <- .draw_weights(groups)[data$g]
  cbind(id = seq_len(nrow(data)), data)
}

cases <- list(
  r01 = list(kind = "regression", family = "gaussian", se = "standard",
             data = simulate_regression(20261101L, "gaussian")),
  r02 = list(kind = "regression", family = "gaussian", se = "robust",
             data = simulate_regression(20261101L, "gaussian")),
  r03 = list(kind = "regression", family = "poisson", se = "standard",
             data = simulate_regression(20261103L, "poisson")),
  r04 = list(kind = "regression", family = "binomial", se = "standard",
             data = simulate_repeated_binary(20261104L)),
  r05 = list(kind = "regression", family = "gaussian", se = "weighted",
             data = simulate_regression(20261105L, "gaussian")),
  w01 = list(kind = "profiles", data = simulate_profiles(20261111L)),
  w02 = list(kind = "transitions", data = simulate_transitions(20261112L)),
  w03 = list(kind = "two_level", data = simulate_two_level(20261113L)))

.header <- function(name, standard_errors = "standarderrors") {
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=2;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261101 sets=50 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;",
    "   output", sprintf("      parameters=first %s;", standard_errors),
    sprintf("   outfile '%s_posteriors.txt' classification keep id;", name))
}

syntax <- function(name, case) {
  switch(case$kind,
    regression = {
      scale <- switch(case$family, gaussian = "continuous", poisson = "poisson",
                      binomial = "nominal coding=first")
      errors <- switch(case$se, standard = "standarderrors",
                       robust = "standarderrors=robust", weighted = "standarderrors=robust")
      c(.header(name, errors), "variables",
        if (identical(case$family, "binomial")) "   caseid id;",
        if (identical(case$se, "weighted")) "   samplingweight w rescale;",
        sprintf("   dependent y %s;", scale), "   independent x, z;",
        "   latent Cluster nominal 2 coding=first;",
        "equations", "   Cluster <- 1 + z;", "   y <- 1 | Cluster + x | Cluster;",
        if (identical(case$family, "gaussian")) "   y | Cluster;",
        "end model")
    },
    profiles = c(.header(name, "standarderrors=robust"), "variables", "   samplingweight w rescale;",
                 "   dependent y1 continuous, y2 continuous, y3 continuous;",
                 "   latent Cluster nominal 3;", "equations", "   Cluster <- 1;",
                 sprintf("   %s <- 1 + Cluster;", c("y1", "y2", "y3")),
                 sprintf("   %s | Cluster;", c("y1", "y2", "y3")), "end model"),
    transitions = c(.header(name, "standarderrors=robust"), "variables", "   caseid id;",
                    "   samplingweight w rescale;",
                    "   dependent y1 continuous, y2 continuous;",
                    "   latent State nominal dynamic 2;", "equations",
                    "   State[=0] <- 1;", "   State <- 1 | State[-1];",
                    sprintf("   %s <- 1 + State;", c("y1", "y2")),
                    sprintf("   %s | State;", c("y1", "y2")), "end model"),
    two_level = c(sub("keep id;", "keep id g;", .header(name, "standarderrors=robust"),
                      fixed = TRUE),
                  "variables", "   groupid g;", "   samplingweight w rescale;",
                  "   dependent y1 continuous, y2 continuous;",
                  "   latent GClass group nominal 2, Cluster nominal 2;", "equations",
                  "   GClass <- 1;", "   Cluster <- 1 + GClass;",
                  sprintf("   %s <- 1 + Cluster;", c("y1", "y2")),
                  sprintf("   %s | Cluster;", c("y1", "y2")), "end model"))
}

invisible(lapply(names(cases), function(name) {
  lines <- syntax(name, cases[[name]])
  data <- cases[[name]]$data
  lapply(c(kit_dir, wine_dir), function(directory) {
    utils::write.table(data, file.path(directory, paste0(name, ".dat")), sep = "\t",
                       row.names = FALSE, quote = FALSE)
    writeLines(lines, file.path(directory, paste0(name, ".lgs")))
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("wrote %d cases to %s and %s\n", length(cases), kit_dir, wine_dir))
