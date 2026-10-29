# Builds the Latent GOLD kit for the families no Latent GOLD comparison covered:
# growth mixtures (LCGA and random-effect GMM), single-level lca() and lpa()
# without weights, lta() with Gaussian, binary and covariate-driven
# transitions, and class enumeration for lca() and lpa(). Run from the project
# root:
#   Rscript equivalence/latentgold-coverage/make-kit.R
# then equivalence/latentgold-coverage/run-kit.sh, then compare.R.
#
# Writes kit/<case>.dat and kit/<case>.lgs, the same files to the Wine drive,
# and kit/cases.rds with each case's data and settings. The latents fits are
# made in compare.R from those data, so the comparison reads the package as it
# is now.

kit_dir <- file.path("equivalence", "latentgold-coverage", "kit")
wine_dir <- file.path(Sys.getenv("HOME"), ".wine-latentgold", "drive_c", "lgcoverage")
dir.create(kit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wine_dir, showWarnings = FALSE, recursive = TRUE)

.binary_items <- function(class, probabilities) {
  values <- vapply(seq_len(ncol(probabilities)), function(item) {
    as.integer(stats::runif(length(class)) < probabilities[class, item])
  }, integer(length(class)))
  colnames(values) <- paste0("u", seq_len(ncol(probabilities)))
  as.data.frame(values)
}

simulate_lca_binary <- function(seed, n = 800L) {
  set.seed(seed)
  class <- sample.int(3L, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  probabilities <- rbind(c(0.9, 0.85, 0.8, 0.85, 0.9, 0.8),
                         c(0.2, 0.15, 0.25, 0.8, 0.85, 0.75),
                         c(0.1, 0.15, 0.2, 0.1, 0.2, 0.15))
  cbind(id = seq_len(n), .binary_items(class, probabilities))
}

simulate_lca_nominal <- function(seed, n = 600L) {
  set.seed(seed)
  class <- 1L + (stats::runif(n) < 0.4)
  probabilities <- list(rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6)),
                        rbind(c(0.6, 0.3, 0.1), c(0.2, 0.2, 0.6)),
                        rbind(c(0.2, 0.7, 0.1), c(0.6, 0.1, 0.3)),
                        rbind(c(0.5, 0.4, 0.1), c(0.1, 0.4, 0.5)))
  items <- lapply(probabilities, function(p) {
    vapply(class, function(k) sample.int(3L, 1L, prob = p[k, ]), integer(1))
  })
  names(items) <- paste0("u", seq_along(items))
  data.frame(id = seq_len(n), items)
}

simulate_lca_covariate <- function(seed, n = 700L) {
  set.seed(seed)
  x <- round(stats::rnorm(n), 4)
  class <- 1L + (stats::runif(n) < stats::plogis(-0.3 + 1.1 * x))
  probabilities <- rbind(c(0.85, 0.8, 0.9, 0.75, 0.8), c(0.2, 0.25, 0.15, 0.3, 0.2))
  cbind(id = seq_len(n), .binary_items(class, probabilities), x = x)
}

simulate_profiles <- function(seed, n = 600L, equal = FALSE) {
  set.seed(seed)
  class <- sample.int(3L, n, replace = TRUE, prob = c(0.45, 0.35, 0.2))
  means <- rbind(c(-1, 0, 1), c(1, 1.2, -0.6), c(0, -1.5, 0.4))
  sds <- if (equal) matrix(c(0.7, 0.8, 0.6), 3L, 3L, byrow = TRUE) else
    rbind(c(0.7, 0.8, 0.6), c(0.9, 0.6, 0.8), c(0.6, 0.7, 1.0))
  values <- means[class, ] + matrix(stats::rnorm(n * 3L), n, 3L) * sds[class, ]
  data.frame(id = seq_len(n), y1 = round(values[, 1L], 4), y2 = round(values[, 2L], 4),
             y3 = round(values[, 3L], 4))
}

simulate_correlated <- function(seed, n = 600L) {
  set.seed(seed)
  class <- 1L + (stats::runif(n) < 0.45)
  means <- rbind(c(-1, 0.5), c(1, -0.5))
  covariances <- list(matrix(c(0.6, 0.3, 0.3, 0.7), 2L), matrix(c(0.9, -0.4, -0.4, 0.6), 2L))
  values <- t(vapply(class, function(k) {
    means[k, ] + drop(t(chol(covariances[[k]])) %*% stats::rnorm(2L))
  }, numeric(2)))
  data.frame(id = seq_len(n), y1 = round(values[, 1L], 4), y2 = round(values[, 2L], 4))
}

.markov_states <- function(persons, occasions, initial, transition) {
  state <- matrix(0L, persons, occasions)
  state[, 1L] <- 1L + (stats::runif(persons) < initial)
  invisible(lapply(seq_len(occasions - 1L), function(t) {
    state[, t + 1L] <<- 1L + (stats::runif(persons) < transition(state[, t], t))
  }))
  state
}

simulate_lta <- function(seed, kind, persons = 400L, occasions = 4L) {
  set.seed(seed)
  z <- round(stats::rnorm(persons), 4)
  transition <- if (identical(kind, "covariate")) {
    function(previous, t) stats::plogis(ifelse(previous == 1L, -1.6 + 0.9 * z, 1.2 - 0.8 * z))
  } else function(previous, t) c(0.15, 0.75)[previous]
  state <- .markov_states(persons, occasions, 0.4, transition)
  long <- data.frame(id = rep(seq_len(persons), each = occasions),
                     time = rep(seq_len(occasions), times = persons),
                     state = as.vector(t(state)), z = rep(z, each = occasions))
  if (identical(kind, "binary")) {
    probabilities <- rbind(c(0.85, 0.8, 0.9), c(0.15, 0.25, 0.2))
    long <- cbind(long, .binary_items(long$state, probabilities))
  } else {
    means <- rbind(c(-0.8, 0.5), c(0.9, -0.6))
    long$y1 <- round(means[long$state, 1L] + stats::rnorm(nrow(long), sd = 0.7), 4)
    long$y2 <- round(means[long$state, 2L] + stats::rnorm(nrow(long), sd = 0.8), 4)
  }
  long$state <- NULL
  long
}

# Persons in two trajectory classes, five equally spaced waves (time 0 to 4).
# `random`: "none" (LCGA), "intercept" or "slope" (random intercept and slope).
simulate_growth <- function(seed, random, persons = 400L, waves = 5L) {
  set.seed(seed)
  class <- 1L + (stats::runif(persons) < 0.4)
  intercept <- c(1, 3)[class]
  slope <- c(0.6, -0.4)[class]
  covariance <- list(matrix(c(0.5, 0.05, 0.05, 0.08), 2L), matrix(c(0.8, -0.1, -0.1, 0.12), 2L))
  effects <- t(vapply(class, function(k) {
    switch(random, none = c(0, 0),
           intercept = c(stats::rnorm(1L, sd = 0.7), 0),
           slope = drop(t(chol(covariance[[k]])) %*% stats::rnorm(2L)))
  }, numeric(2)))
  time <- seq_len(waves) - 1L
  long <- data.frame(id = rep(seq_len(persons), each = waves), time = rep(time, persons))
  person <- long$id
  long$y <- round(intercept[person] + effects[person, 1L] +
                    (slope[person] + effects[person, 2L]) * long$time +
                    stats::rnorm(nrow(long), sd = c(0.5, 0.7)[class[person]]), 4)
  long
}

.to_wide <- function(long) {
  wide <- stats::reshape(long, idvar = "id", timevar = "time", direction = "wide", sep = "")
  names(wide) <- c("id", paste0("y", seq_len(ncol(wide) - 1L)))
  row.names(wide) <- NULL
  wide
}

# Mixture regressions beyond the wider kit's Gaussian, Poisson and logistic:
# negative binomial (NB2, dispersion by class) and cumulative-logit ordinal.
simulate_family_regression <- function(seed, family, n = 900L) {
  set.seed(seed)
  x <- round(stats::rnorm(n), 4)
  class <- 1L + (stats::runif(n) < 0.4)
  eta <- c(0.4, 1.6)[class] + c(0.5, -0.4)[class] * x
  y <- switch(family,
    negative_binomial = stats::rnbinom(n, size = c(2, 5)[class], mu = exp(eta)),
    # Opposite slopes keep the two classes apart; with similar slopes a class
    # collapses onto a deterministic response (a separation boundary).
    ordinal = {
      cuts <- rbind(c(-1.2, 0, 1.2), c(-0.8, 0.4, 1.6))
      latent <- c(1.5, -1.5)[class] * x + stats::rlogis(n)
      1L + vapply(seq_len(n), function(i) sum(latent[[i]] > cuts[class[[i]], ]), integer(1))
    })
  data.frame(id = seq_len(n), y = y, x = x)
}

lca_binary <- simulate_lca_binary(20261201L)
profiles <- simulate_profiles(20261211L)
cases <- c(
  list(
    g01 = list(kind = "lcga", data = simulate_growth(20261221L, "none")),
    g02 = list(kind = "gmm_slope", data = simulate_growth(20261222L, "slope")),
    g03 = list(kind = "gmm_intercept", data = simulate_growth(20261223L, "intercept")),
    c01 = list(kind = "lca", items = paste0("u", 1:6), classes = 3L, data = lca_binary),
    c02 = list(kind = "lca", items = paste0("u", 1:4), classes = 2L,
               data = simulate_lca_nominal(20261202L)),
    c03 = list(kind = "lca_covariate", items = paste0("u", 1:5), classes = 2L,
               data = simulate_lca_covariate(20261203L)),
    p01 = list(kind = "lpa", model = "VVI", classes = 3L, data = profiles),
    p02 = list(kind = "lpa", model = "EEI", classes = 3L,
               data = simulate_profiles(20261212L, equal = TRUE)),
    p03 = list(kind = "lpa", model = "VVV", classes = 2L, data = simulate_correlated(20261213L)),
    t01 = list(kind = "lta", items = c("y1", "y2"), categorical = FALSE,
               data = simulate_lta(20261231L, "gaussian")),
    t02 = list(kind = "lta", items = c("u1", "u2", "u3"), categorical = TRUE,
               data = simulate_lta(20261232L, "binary")),
    t03 = list(kind = "lta_covariate", items = c("y1", "y2"), categorical = FALSE,
               data = simulate_lta(20261233L, "covariate")),
    r06 = list(kind = "nb_regression", data = simulate_family_regression(20261241L, "negative_binomial")),
    r07 = list(kind = "ordinal_regression",
               data = simulate_family_regression(20261242L, "ordinal", n = 1500L))),
  stats::setNames(lapply(1:4, function(k) list(kind = "lca_enumeration", items = paste0("u", 1:6),
                                               classes = k, data = lca_binary)),
                  sprintf("e1k%d", 1:4)),
  stats::setNames(lapply(1:4, function(k) list(kind = "lpa_enumeration", model = "VVI",
                                               classes = k, data = profiles)),
                  sprintf("e2k%d", 1:4)))

.header <- function(name, keep = "id") {
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", name), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=1;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261201 sets=50 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;",
    "   output", "      parameters=first standarderrors;",
    sprintf("   outfile '%s_posteriors.txt' classification keep %s;", name, keep))
}

.lpa_lines <- function(vars, model, classes) {
  variances <- switch(model, EEI = sprintf("   %s;", vars),
                      sprintf("   %s | Cluster;", vars))
  c(sprintf("   latent Cluster nominal %d coding=first;", classes), "equations",
    "   Cluster <- 1;",
    if (classes > 1L) sprintf("   %s <- 1 + Cluster;", vars) else sprintf("   %s <- 1;", vars),
    variances,
    if (identical(model, "VVV")) "   y1 <-> y2 | Cluster;")
}

syntax <- function(name, case) {
  data_vars <- setdiff(names(case$data), c("id", "time"))
  switch(case$kind,
    lcga = c(.header(name), "variables", "   caseid id;", "   dependent y continuous;",
             "   independent time;", "   latent Class nominal 2 coding=first;", "equations",
             "   Class <- 1;", "   y <- 1 | Class + time | Class;", "   y | Class;",
             "end model"),
    gmm_slope = , gmm_intercept = {
      waves <- seq_len(ncol(.to_wide(case$data)) - 1L)
      random_slope <- identical(case$kind, "gmm_slope")
      c(.header(name), "variables",
        sprintf("   dependent %s;", paste(sprintf("y%d continuous", waves), collapse = ", ")),
        "   latent Class nominal 2 coding=first, F1 continuous, F2 continuous;", "equations",
        "   Class <- 1;",
        if (random_slope) c("   F1 | Class;", "   F2 | Class;", "   F1 <-> F2 | Class;") else
          c("   F1;", "   (0) F2;", "   (0) F1 <-> F2;"),
        "   F1 <- 1 | Class;", "   F2 <- 1 | Class;",
        sprintf("   y%d <- (1) F1 + (%d) F2;", waves, waves - 1L),
        sprintf("   (e) y%d | Class;", waves), "end model")
    },
    nb_regression = , ordinal_regression = c(
      .header(name), "variables",
      sprintf("   dependent y %s;", if (identical(case$kind, "nb_regression")) "poisson overdispersed" else
        "cumlogit"),
      "   independent x;", "   latent Cluster nominal 2 coding=first;", "equations",
      "   Cluster <- 1;", "   y <- 1 | Cluster + x | Cluster;",
      if (identical(case$kind, "nb_regression")) "   y | Cluster;", "end model"),
    lca = , lca_enumeration = c(
      .header(name), "variables",
      sprintf("   dependent %s;", paste(sprintf("%s nominal coding=first", case$items), collapse = ", ")),
      sprintf("   latent Cluster nominal %d coding=first;", case$classes), "equations",
      "   Cluster <- 1;",
      sprintf("   %s <- 1 %s;", case$items, if (case$classes > 1L) "| Cluster" else ""),
      "end model"),
    lca_covariate = c(
      .header(name), "variables",
      sprintf("   dependent %s;", paste(sprintf("%s nominal coding=first", case$items), collapse = ", ")),
      "   independent x;",
      sprintf("   latent Cluster nominal %d coding=first;", case$classes), "equations",
      "   Cluster <- 1 + x;", sprintf("   %s <- 1 | Cluster;", case$items), "end model"),
    lpa = , lpa_enumeration = c(
      .header(name), "variables",
      sprintf("   dependent %s;", paste(sprintf("%s continuous", data_vars), collapse = ", ")),
      .lpa_lines(data_vars, case$model, case$classes), "end model"),
    lta = , lta_covariate = {
      scale <- if (case$categorical) "nominal coding=first" else "continuous"
      covariate <- identical(case$kind, "lta_covariate")
      c(.header(name), "variables", "   caseid id;",
        sprintf("   dependent %s;", paste(sprintf("%s %s", case$items, scale), collapse = ", ")),
        if (covariate) "   independent z;",
        "   latent State nominal dynamic 2 coding=first;", "equations",
        "   State[=0] <- 1;",
        if (covariate) "   State <- 1 | State[-1] + z | State[-1];" else "   State <- 1 | State[-1];",
        if (case$categorical) sprintf("   %s <- 1 | State;", case$items) else
          c(sprintf("   %s <- 1 + State;", case$items), sprintf("   %s | State;", case$items)),
        "end model")
    })
}

.lg_data <- function(case) {
  if (case$kind %in% c("gmm_slope", "gmm_intercept")) .to_wide(case$data) else case$data
}

invisible(lapply(names(cases), function(name) {
  lines <- syntax(name, cases[[name]])
  data <- .lg_data(cases[[name]])
  lapply(c(kit_dir, wine_dir), function(directory) {
    utils::write.table(data, file.path(directory, paste0(name, ".dat")), sep = "\t",
                       row.names = FALSE, quote = FALSE)
    writeLines(lines, file.path(directory, paste0(name, ".lgs")))
  })
}))
saveRDS(cases, file.path(kit_dir, "cases.rds"))
cat(sprintf("wrote %d cases to %s and %s\n", length(cases), kit_dir, wine_dir))
