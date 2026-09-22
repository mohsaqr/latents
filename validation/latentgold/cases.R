# The Latent GOLD comparison cases: data, this package's fit, and the matching
# Latent GOLD syntax, defined once.
#
# Sourced by make-kit.R (which writes the kit) and compare.R (which reads the
# returned output). Each case is built from a fixed seed, so the data written
# into the kit and the data the targets were computed from are the same rows.
#
# Latent GOLD is the reference implementation of the multilevel latent class
# model (Vermunt 2003) and of the three-step corrections (Vermunt 2010; Bakk,
# Tekle & Vermunt 2013). See README.md for which feature each case anchors.

# ---------------------------------------------------------------------------
# Latent GOLD syntax
# ---------------------------------------------------------------------------

#' The options block shared by every model
#'
#' Two settings carry the comparison and must not be relaxed:
#' - `bayes ... = 0`. Latent GOLD's default Bayes constants put a prior on the
#'   variances and response probabilities, so its default estimates are
#'   posterior modes, not maximum likelihood. This package is ML.
#' - many random start sets with a tight tolerance, so both programs report the
#'   global optimum rather than whichever local one a few starts reach.
#'
#' @param outfile Name of the posterior file to write, or `NULL`.
#' @param keep Variables to carry into the posterior file.
#' @param missing `"excludeall"` (listwise) or `"includeall"` (FIML).
#' @param step3 A Step-3 specification such as `"modal bch"`, or `NULL`.
#' @return A character vector of syntax lines.
.lg_options <- function(outfile = NULL, keep = "id", missing = "excludeall",
                        step3 = NULL) {
  c("options",
    "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-010 emtolerance=1e-008 emiterations=5000 nriterations=500;",
    "   startvalues",
    "      seed=20260922 sets=100 tolerance=1e-005 iterations=250;",
    "   bayes",
    "      categorical=0 variances=0 latent=0 poisson=0;",
    sprintf("   missing %s;", missing),
    if (!is.null(step3)) sprintf("   step3 %s;", step3),
    "   output",
    "      parameters=first standarderrors profile probmeans=posterior",
    "      bivariateresiduals estimatedvalues=model;",
    if (!is.null(outfile))
      sprintf("   outfile '%s' classification keep %s;", outfile,
              paste(keep, collapse = " ")))
}

#' One complete Latent GOLD syntax file
#'
#' @param infile The tab-separated data file the model reads.
#' @param models A list of character vectors, one per `model ... end model`
#'   block, each holding that block's `title`, `options`, `variables` and
#'   `equations` lines.
#' @return A character vector: the file's lines.
.lg_file <- function(infile, models) {
  c("//LG6.1//",
    "version = 6.1",
    sprintf("infile '%s'", infile),
    "",
    unlist(lapply(models, function(block) c("model", block, "end model", ""))))
}

#' A dependent-variable declaration
#' @param continuous,nominal Names of continuous and nominal indicators.
#' @return One `dependent` line.
.lg_dependents <- function(continuous = character(), nominal = character()) {
  declared <- c(sprintf("%s continuous", continuous), sprintf("%s nominal", nominal))
  sprintf("   dependent %s;", paste(declared, collapse = ", "))
}

#' Measurement equations for a set of indicators
#'
#' @param continuous,nominal Indicator names.
#' @param latent The latent variable they load on.
#' @param variances `"varying"` for class-specific residual variances,
#'   `"equal"` for one variance per indicator.
#' @param covariances Pairs of continuous indicators whose residual covariance
#'   is estimated, class-specific, as a two-column character matrix.
#' @return Equation lines.
.lg_measurement <- function(continuous = character(), nominal = character(),
                            latent = "Cluster", variances = "varying",
                            covariances = NULL) {
  c(sprintf("   %s <- 1 + %s;", c(continuous, nominal), latent),
    if (identical(variances, "varying") && length(continuous) > 0L)
      sprintf("   %s | %s;", continuous, latent),
    if (!is.null(covariances))
      sprintf("   %s <-> %s | %s;", covariances[, 1L], covariances[, 2L], latent))
}

#' Predictors as terms of an equation, or nothing at all
#'
#' `paste0(" + ", character(0), collapse = "")` is `" + "`, not `""`: paste
#' treats a zero-length argument as `""` unless `recycle0 = TRUE`. That wrote
#' `GClass <- 1 + ;` into every case without covariates, which Latent GOLD
#' rejects with "expecting variable, encountered \";\"".
#'
#' @param predictors Character vector of predictor names, possibly empty.
#' @return `""`, or ` + a + b` for the predictors given.
.lg_terms <- function(predictors) {
  if (length(predictors) == 0L) return("")
  paste0(" + ", predictors, collapse = "")
}

#' A two-level cluster model: profiles within discrete group classes
#'
#' @param name Case name, used for the title and the posterior file.
#' @param continuous,nominal Indicator names.
#' @param variances,covariances See [.lg_measurement()].
#' @param missing See [.lg_options()].
#' @param profile_covariates Individual-level predictors of the profile, with a
#'   slope common to every group class, as this package fits them.
#' @param group_covariates Group-level predictors of the group class.
#' @param n_profiles,n_group_classes Class counts.
#' @return One model block.
.lg_twolevel <- function(name, continuous = character(), nominal = character(),
                         variances = "varying", covariances = NULL,
                         missing = "excludeall", profile_covariates = character(),
                         group_covariates = character(), n_profiles = 2L,
                         n_group_classes = 2L) {
  predictors <- c(profile_covariates, group_covariates)
  c(sprintf("title '%s';", name),
    .lg_options(outfile = sprintf("%s_posteriors.txt", name), keep = c("id", "g"),
                missing = missing),
    "variables",
    "   groupid g;",
    if (length(predictors) > 0L)
      sprintf("   independent %s;", paste(predictors, collapse = ", ")),
    .lg_dependents(continuous, nominal),
    "   latent",
    sprintf("      GClass group nominal %d,", n_group_classes),
    sprintf("      Cluster nominal %d;", n_profiles),
    "equations",
    sprintf("   GClass <- 1%s;", .lg_terms(group_covariates)),
    sprintf("   Cluster <- 1 + GClass%s;", .lg_terms(profile_covariates)),
    .lg_measurement(continuous, nominal, "Cluster", variances, covariances))
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#' Evaluate `code` under a fixed seed, leaving the caller's RNG state as it was
#'
#' @param seed RNG seed.
#' @param code An expression, evaluated lazily after the seed is set.
#' @return The value of `code`.
.lg_with_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed
  on.exit(if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE, after = FALSE)
  set.seed(seed)
  code
}

#' Two-level data: profiles within discrete group classes
#'
#' Group class 1 is dominated by profile 1 and class 2 by profile 2, so both
#' levels are identified and the posteriors at each are well away from 0.5.
#'
#' @param seed RNG seed; the data are a pure function of it.
#' @param n_groups,group_size Layout.
#' @param means A profile-by-indicator matrix of means.
#' @param correlation Within-profile residual correlation of the first two
#'   indicators.
#' @return A data frame with `id`, `g`, `profile`, `group_class` and the
#'   indicators `y1`, `y2`, ...
.lg_twolevel_data <- function(seed, n_groups = 60L, group_size = 10L,
                              means = rbind(c(0, 0, 0), c(2, 1.5, 1)),
                              correlation = 0) {
  .lg_with_seed(seed, {
    n <- n_groups * group_size
    group_class <- rep(sample(1:2, n_groups, replace = TRUE), each = group_size)
    profile <- ifelse(stats::runif(n) < ifelse(group_class == 1L, 0.8, 0.2), 1L, 2L)
    noise <- matrix(stats::rnorm(n * ncol(means)), n)
    # A shared term gives the first two indicators the requested correlation.
    noise[, 2L] <- correlation * noise[, 1L] + sqrt(1 - correlation^2) * noise[, 2L]
    indicators <- means[profile, , drop = FALSE] + noise
    colnames(indicators) <- sprintf("y%d", seq_len(ncol(means)))
    data.frame(id = seq_len(n), g = rep(seq_len(n_groups), each = group_size),
               profile = profile, group_class = group_class, indicators)
  })
}

#' Categorical items drawn from profile-specific response probabilities
#' @param profile Integer profile of each row.
#' @param probabilities A list of profile-by-category matrices, one per item.
#' @return A data frame of integer-coded items `c1`, `c2`, ...
.lg_items <- function(profile, probabilities) {
  items <- lapply(probabilities, function(item) {
    vapply(profile, function(k) sample.int(ncol(item), 1L, prob = item[k, ]),
           integer(1))
  })
  stats::setNames(as.data.frame(items), sprintf("c%d", seq_along(probabilities)))
}

# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

#' Fit this package the way every case needs it fitted
#'
#' Fifty starts, so the target is the global optimum, and a tolerance far
#' tighter than the package's default. At `tol = 1e-10` the posteriors still
#' sat up to 1.7e-4 from Latent GOLD's, and tightening Latent GOLD did not move
#' them at all; at `tol = 1e-14` the same comparisons fall to about 2e-6. The
#' gap was this package's stopping rule, not a difference between the two
#' programs, so the target is converged rather than the tolerance widened.
#' Non-convergence stops the kit rather than shipping an unconverged target.
#'
#' @param ... Passed to [multilpa()].
#' @return The fit.
.lg_fit <- function(...) {
  fit <- multilpa(..., n_starts = 50L, max_iter = 200000L, tol = 1e-14, seed = 1L)
  stopifnot("a comparison target did not converge" = isTRUE(fit$converged))
  fit
}

#' Every comparison case
#'
#' @return A named list of cases. Each has `description`, `anchors` (what it
#'   validates), `data`, `fit` (a function returning this package's fit),
#'   `files` (named list: Latent GOLD syntax file name -> lines, in the order
#'   they must be estimated) and `data_file`.
lg_cases <- function() {
  continuous <- .lg_twolevel_data(seed = 101L)
  correlated <- .lg_twolevel_data(seed = 103L, means = rbind(c(0, 0), c(2, 1.5)),
                                  correlation = 0.5)
  dependent <- .lg_twolevel_data(seed = 108L, n_groups = 40L)
  # y1 and y2 share a term the profiles do not explain: local dependence.
  dependent$y2 <- dependent$y2 + 0.8 * (dependent$y1 - ave(dependent$y1, dependent$profile))

  categorical <- .lg_twolevel_data(seed = 104L, n_groups = 80L, group_size = 8L)
  item <- rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6))
  categorical <- cbind(categorical[c("id", "g", "profile", "group_class")],
                       .lg_with_seed(1104L, .lg_items(categorical$profile,
                                                      list(item, item[, 3:1], item, item))))
  mixed <- .lg_twolevel_data(seed = 105L, means = rbind(c(0, 0), c(2, 1.5)))
  mixed <- cbind(mixed, .lg_with_seed(1105L, .lg_items(mixed$profile,
                                                       list(item, item[, 3:1]))))

  missing_data <- continuous
  # Ten per cent of indicator values missing completely at random, never a
  # whole row, so every row still enters the FIML likelihood.
  holes <- .lg_with_seed(1106L, matrix(stats::runif(nrow(missing_data) * 3L) < 0.10,
                                       ncol = 3L))
  holes[rowSums(holes) == 3L, 3L] <- FALSE
  missing_data[c("y1", "y2", "y3")][holes] <- NA

  covariates <- continuous
  covariates$x <- .lg_with_seed(1107L, stats::rnorm(nrow(covariates)))
  covariates$w <- rep(.lg_with_seed(1108L, stats::rnorm(60L)), each = 10L)

  steps <- .lg_step_data(seed = 109L)
  sequences <- .lg_sequence_data(seed = 110L)

  cont3 <- c("y1", "y2", "y3")
  list(
    c01_continuous_varying = list(
      description = "Two-level, 2 profiles x 2 group classes, three continuous indicators, class-specific variances.",
      anchors = c("likelihood", "posteriors", "group posteriors", "means and variances", "standard errors"),
      data = continuous,
      fit = function() .lg_fit(continuous, cont3, "g", n_profiles = 2, n_group_classes = 2),
      files = list(c01_continuous_varying.lgs = .lg_file("c01.dat", list(
        .lg_twolevel("c01_continuous_varying", continuous = cont3)))),
      data_file = "c01.dat"),
    c02_continuous_equal = list(
      description = "As c01 with one residual variance per indicator, shared across profiles.",
      anchors = c("likelihood", "posteriors", "group posteriors", "means and variances", "standard errors"),
      data = continuous,
      fit = function() .lg_fit(continuous, cont3, "g", n_profiles = 2, n_group_classes = 2,
                               variance_model = "equal"),
      files = list(c02_continuous_equal.lgs = .lg_file("c01.dat", list(
        .lg_twolevel("c02_continuous_equal", continuous = cont3, variances = "equal")))),
      data_file = "c01.dat"),
    c03_full_covariance = list(
      description = "Two-level, two correlated continuous indicators, class-specific full covariance.",
      anchors = c("likelihood", "posteriors", "group posteriors"),
      data = correlated,
      fit = function() .lg_fit(correlated, c("y1", "y2"), "g", n_profiles = 2,
                               n_group_classes = 2, covariance_model = "full"),
      files = list(c03_full_covariance.lgs = .lg_file("c03.dat", list(
        .lg_twolevel("c03_full_covariance", continuous = c("y1", "y2"),
                     covariances = matrix(c("y1", "y2"), 1L))))),
      data_file = "c03.dat"),
    c04_categorical = list(
      description = "Two-level, four nominal items with three categories each.",
      anchors = c("likelihood", "posteriors", "group posteriors"),
      data = categorical,
      fit = function() .lg_fit(categorical, sprintf("c%d", 1:4), "g", n_profiles = 2,
                               n_group_classes = 2, categorical = sprintf("c%d", 1:4)),
      files = list(c04_categorical.lgs = .lg_file("c04.dat", list(
        .lg_twolevel("c04_categorical", nominal = sprintf("c%d", 1:4))))),
      data_file = "c04.dat"),
    c05_mixed = list(
      description = "Two-level, two continuous indicators and two nominal items.",
      anchors = c("likelihood", "posteriors", "group posteriors"),
      data = mixed,
      fit = function() .lg_fit(mixed, c("y1", "y2", "c1", "c2"), "g", n_profiles = 2,
                               n_group_classes = 2, categorical = c("c1", "c2")),
      files = list(c05_mixed.lgs = .lg_file("c05.dat", list(
        .lg_twolevel("c05_mixed", continuous = c("y1", "y2"), nominal = c("c1", "c2"))))),
      data_file = "c05.dat"),
    c06_missing_fiml = list(
      description = "c01 with 10% of indicator values missing completely at random, full-information ML.",
      anchors = c("likelihood", "posteriors", "group posteriors"),
      data = missing_data,
      fit = function() .lg_fit(missing_data, cont3, "g", n_profiles = 2,
                               n_group_classes = 2, missing = "fiml"),
      files = list(c06_missing_fiml.lgs = .lg_file("c06.dat", list(
        .lg_twolevel("c06_missing_fiml", continuous = cont3, missing = "includeall")))),
      data_file = "c06.dat"),
    c07_covariates = list(
      description = "c01 with an individual predictor of the profile and a group predictor of the group class.",
      anchors = c("likelihood", "posteriors", "group posteriors", "standard errors"),
      data = covariates,
      fit = function() .lg_fit(covariates, cont3, "g", n_profiles = 2, n_group_classes = 2,
                               profile_covariates = "x", group_covariates = "w"),
      files = list(c07_covariates.lgs = .lg_file("c07.dat", list(
        .lg_twolevel("c07_covariates", continuous = cont3,
                     profile_covariates = "x", group_covariates = "w")))),
      data_file = "c07.dat"),
    c08_bivariate_residuals = list(
      description = "Two-level, y1 and y2 locally dependent within profile; y3 independent.",
      anchors = c("likelihood", "posteriors", "bivariate residual ranking"),
      data = dependent,
      fit = function() .lg_fit(dependent, cont3, "g", n_profiles = 2, n_group_classes = 2),
      files = list(c08_bivariate_residuals.lgs = .lg_file("c08.dat", list(
        .lg_twolevel("c08_bivariate_residuals", continuous = cont3)))),
      data_file = "c08.dat"),
    c09_three_step = .lg_step_case(steps),
    c10_lta = .lg_sequence_case(sequences, n_group_classes = 1L),
    c11_lta_mixture = .lg_sequence_case(sequences, n_group_classes = 2L))
}

# ---------------------------------------------------------------------------
# Three-step
# ---------------------------------------------------------------------------

#' Single-level data for the three-step corrections
#'
#' Membership follows `plogis(-0.3 + 1.2 * x)`, and the distal outcome is 10
#' in profile 1 and 0 in profile 2, so both packages can also be scored
#' against the generating values.
#'
#' @param seed RNG seed.
#' @return A data frame with `id`, `g`, `profile`, `x`, `y1`, `y2`, `distal`.
.lg_step_data <- function(seed) {
  .lg_with_seed(seed, {
    n <- 600L
    x <- stats::rnorm(n)
    profile <- ifelse(stats::runif(n) < stats::plogis(-0.3 + 1.2 * x), 1L, 2L)
    data.frame(id = seq_len(n), g = seq_len(n), profile = profile, x = x,
               y1 = stats::rnorm(n, ifelse(profile == 1L, 2, 0)),
               y2 = stats::rnorm(n, ifelse(profile == 1L, 1.5, 0)),
               distal = stats::rnorm(n, ifelse(profile == 1L, 10, 0), 2))
  })
}

#' The three-step case: one measurement model, then four Step-3 analyses
#'
#' Latent GOLD's Step-3 reads the posterior file the step-one model writes, so
#' the step-one file must be estimated first, and `files` lists it first. The mapping
#' from this package's methods to Latent GOLD's Step-3 options:
#'
#' | this package | Latent GOLD |
#' |---|---|
#' | `three_step(method = "bch")` | `step3 modal bch` |
#' | `three_step(method = "modal")` | `step3 modal none` |
#' | `three_step(method = "proportional")` | `step3 proportional none` |
#' | `r3step()` | `step3 modal ml` |
#'
#' @param data From [.lg_step_data()].
#' @return One case.
.lg_step_case <- function(data) {
  posterior_columns <- "( Cluster#1 Cluster#2 )"
  step3_model <- function(name, step3, dependent = NULL, predictor = NULL) {
    c(sprintf("title '%s';", name),
      .lg_options(step3 = step3),
      "variables",
      if (!is.null(predictor)) sprintf("   independent %s;", predictor),
      if (!is.null(dependent)) sprintf("   dependent %s continuous;", dependent),
      sprintf("   latent Cluster nominal posterior = %s;", posterior_columns),
      "equations",
      sprintf("   Cluster <- 1%s;", if (is.null(predictor)) "" else paste0(" + ", predictor)),
      if (!is.null(dependent)) sprintf("   %s <- 1 + Cluster;", dependent))
  }
  step1 <- c("title 'c09_step1';",
             .lg_options(outfile = "c09_step1_posteriors.txt",
                         keep = c("id", "x", "distal")),
             "variables",
             .lg_dependents(c("y1", "y2")),
             "   latent Cluster nominal 2;",
             "equations",
             "   Cluster <- 1;",
             .lg_measurement(c("y1", "y2")))
  list(
    description = "Single-level 2-profile measurement model, then Step-3 for a distal outcome (BCH, modal, proportional) and for a covariate (ML).",
    anchors = c("likelihood", "posteriors", "classification errors", "BCH means", "modal means",
                "proportional means", "R3STEP slope"),
    data = data,
    fit = function() {
      # A single-level fit is deliberate here: Step-3 is defined on it, and
      # this is the one warning the call is expected to raise.
      withCallingHandlers(
        .lg_fit(data, c("y1", "y2"), id = NULL, n_profiles = 2, n_group_classes = 1),
        multilpa_single_level = function(w) invokeRestart("muffleWarning"))
    },
    # One model per file, so each listing holds one result and none has to be
    # split apart by title. Step one must run first: the rest read its file.
    files = list(
      c09_step1.lgs = .lg_file("c09.dat", list(step1)),
      c09_distal_bch.lgs = .lg_file("c09_step1_posteriors.txt", list(
        step3_model("c09_distal_bch", "modal bch", dependent = "distal"))),
      c09_distal_modal.lgs = .lg_file("c09_step1_posteriors.txt", list(
        step3_model("c09_distal_modal", "modal none", dependent = "distal"))),
      c09_distal_proportional.lgs = .lg_file("c09_step1_posteriors.txt", list(
        step3_model("c09_distal_proportional", "proportional none", dependent = "distal"))),
      c09_covariate_ml.lgs = .lg_file("c09_step1_posteriors.txt", list(
        step3_model("c09_covariate_ml", "modal ml", predictor = "x")))),
    data_file = "c09.dat")
}

# ---------------------------------------------------------------------------
# Latent transitions
# ---------------------------------------------------------------------------

#' Sequences for the latent transition cases
#'
#' 120 sequences of four occasions. Sequence class 1 stays in its state with
#' probability 0.85, class 2 with 0.55, so a two-class mixture is identified.
#'
#' @param seed RNG seed.
#' @return One row per occasion: `id`, `g` (the sequence), `t`, `state`,
#'   `sequence_class`, `y1`, `y2`.
.lg_sequence_data <- function(seed) {
  .lg_with_seed(seed, {
    n_sequences <- 120L
    n_occasions <- 4L
    sequence_class <- sample(1:2, n_sequences, replace = TRUE)
    stay <- c(0.85, 0.55)[sequence_class]
    first <- sample(1:2, n_sequences, replace = TRUE)
    # Each column is the previous one, switched with probability 1 - stay.
    states <- Reduce(function(previous, occasion) {
      ifelse(stats::runif(n_sequences) < stay, previous, 3L - previous)
    }, seq_len(n_occasions - 1L), accumulate = TRUE, init = first)
    state <- as.vector(t(do.call(cbind, states)))
    n <- n_sequences * n_occasions
    data.frame(id = seq_len(n), g = rep(seq_len(n_sequences), each = n_occasions),
               t = rep(seq_len(n_occasions), n_sequences), state = state,
               sequence_class = rep(sequence_class, each = n_occasions),
               y1 = stats::rnorm(n, ifelse(state == 1L, 0, 2)),
               y2 = stats::rnorm(n, ifelse(state == 1L, 0, 1.5)))
  })
}

#' A latent transition case, with or without sequence-level classes
#'
#' Latent GOLD's Markov module: `State` is dynamic, its initial distribution
#' and transition matrix each depend on the sequence class `GClass` when there
#' is one, and the measurement is invariant over occasions, as [lta()] fits it.
#'
#' @param data From [.lg_sequence_data()].
#' @param n_group_classes 1 for a plain latent transition model, 2 for the
#'   mixture.
#' @return One case.
.lg_sequence_case <- function(data, n_group_classes) {
  name <- if (n_group_classes == 1L) "c10_lta" else "c11_lta_mixture"
  mixture <- n_group_classes > 1L
  block <- c(sprintf("title '%s';", name),
             .lg_options(outfile = sprintf("%s_posteriors.txt", name), keep = c("id", "t")),
             "variables",
             "   caseid g;",
             .lg_dependents(c("y1", "y2")),
             "   latent",
             if (mixture) "      GClass nominal 2,",
             "      State nominal dynamic 2;",
             "equations",
             if (mixture) "   GClass <- 1;",
             sprintf("   State[=0] <- 1%s;", if (mixture) " | GClass" else ""),
             sprintf("   State <- 1 | State[-1]%s;", if (mixture) " GClass" else ""),
             .lg_measurement(c("y1", "y2"), latent = "State"))
  list(
    description = if (mixture)
      "Latent transitions with two sequence classes, each with its own initial distribution and transition matrix." else
      "Latent transitions: two states, four occasions, first-order homogeneous, invariant measurement.",
    anchors = c("likelihood", "state posteriors", if (mixture) "sequence-class posteriors"),
    data = data,
    fit = function() {
      fit <- lta(data, c("y1", "y2"), "g", n_profiles = 2, time = "t",
                 n_group_classes = n_group_classes, n_starts = 50L,
                 max_iter = 200000L, tol = 1e-14, seed = 1L)
      stopifnot("a comparison target did not converge" = isTRUE(fit$converged))
      fit
    },
    files = stats::setNames(list(.lg_file("c10.dat", list(block))), sprintf("%s.lgs", name)),
    data_file = "c10.dat")
}
