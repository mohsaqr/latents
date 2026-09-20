# Mplus User's Guide example 10.7: two-level latent class analysis with
# categorical latent class indicators and a between-level categorical latent
# variable. This is the nonparametric two-level mixture of Vermunt (2003) in
# Mplus's own notation -- cw(4) profiles measured by ten three-category items,
# cb(5) group classes that shift profile prevalence through `cw#1-cw#3 ON cb`,
# and a measurement model shared across group classes.
#
# Artifacts in validation/mplus/twolevel-lca/ are the genuine Mplus 8.8 files
# from https://www.statmodel.com/usersguide/chap10/chapter10.zip
# (sha256 ae29930be639669311d8652af790c531a861d8aa5414783766cd175f7c32b1d0).
# Reference values are parsed from ex10.7.out rather than transcribed.

#' Read the Mplus output shipped with example 10.7
#' @param path Path to the gzipped Mplus output file.
#' @return Character vector of output lines.
.mplus_output_lines <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection)
}

#' Pull a labelled scalar from an Mplus fit-statistics block
#' @param lines Mplus output lines.
#' @param label Text preceding the value.
#' @return The single numeric value reported on that line.
.mplus_scalar <- function(lines, label) {
  matched <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (length(matched) < 1L) stop(sprintf("Mplus output has no line matching '%s'.", label))
  value <- utils::tail(strsplit(trimws(matched[[1L]]), " +")[[1L]], 1L)
  as.numeric(value)
}

#' Parse the RESULTS IN PROBABILITY SCALE block
#'
#' @param lines Mplus output lines.
#' @param items Character vector of item names, upper case as Mplus prints them.
#' @param n_categories Number of categories per item.
#' @return A named list of `pattern` (a two-column matrix of cb and cw indices)
#'   and `probabilities` (a patterns-by-item-by-category array).
.mplus_probability_scale <- function(lines, items, n_categories) {
  start <- grep("^RESULTS IN PROBABILITY SCALE", lines)
  if (length(start) != 1L) stop("Expected exactly one probability-scale block.")
  block <- lines[seq(start, length(lines))]
  headers <- grep("^Latent Class Pattern ", block)
  pattern <- do.call(rbind, lapply(headers, function(index) {
    as.integer(strsplit(trimws(sub("^Latent Class Pattern ", "", block[[index]])), " +")[[1L]])
  }))
  probabilities <- vapply(headers, function(index) {
    window <- block[seq(index, min(index + 6L * length(items), length(block)))]
    vapply(items, function(item) {
      at <- grep(sprintf("^ %s$", item), window)[[1L]]
      rows <- window[at + seq_len(n_categories)]
      as.numeric(sub("^ +Category [0-9]+ +", "", rows))
    }, numeric(n_categories))
  }, matrix(0, n_categories, length(items)))
  list(pattern = pattern,
       probabilities = aperm(probabilities, c(3L, 2L, 1L)))
}

#' Parse the categorical latent variable logits into prevalence probabilities
#' @param lines Mplus output lines.
#' @param n_profiles Number of within-level classes.
#' @param n_group_classes Number of between-level classes.
#' @return A list with `group_probabilities` and `profile_probabilities`.
.mplus_prevalence <- function(lines, n_profiles, n_group_classes) {
  start <- grep("^Categorical Latent Variables", lines)[[1L]]
  block <- lines[seq(start, length(lines))]
  value_of <- function(name, within) {
    hits <- grep(sprintf("^ +%s +-?[0-9]", name), block)
    if (length(hits) < 1L) stop(sprintf("no estimate printed for %s", name))
    chosen <- if (within) hits[[1L]] else hits[[length(hits)]]
    as.numeric(strsplit(trimws(block[[chosen]]), " +")[[1L]][[2L]])
  }
  intercepts <- vapply(seq_len(n_profiles - 1L),
                       \(s) value_of(sprintf("CW#%d", s), within = TRUE), numeric(1))
  slope_start <- grep("^ CW#1 +ON", block)[[1L]]
  slopes <- vapply(seq_len(n_profiles - 1L), function(s) {
    at <- grep(sprintf("^ CW#%d +ON", s), block)[[1L]]
    window <- block[at + seq_len(n_group_classes)]
    vapply(seq_len(n_group_classes - 1L), function(j) {
      as.numeric(strsplit(trimws(window[[j]]), " +")[[1L]][[2L]])
    }, numeric(1))
  }, numeric(n_group_classes - 1L))
  means_at <- grep("^ Means", block)
  means_at <- means_at[means_at > slope_start][[1L]]
  group_logits <- vapply(seq_len(n_group_classes - 1L), function(j) {
    as.numeric(strsplit(trimws(block[[means_at + j]]), " +")[[1L]][[2L]])
  }, numeric(1))
  softmax <- function(logits) exp(logits) / sum(exp(logits))
  profile_probabilities <- t(vapply(seq_len(n_group_classes), function(j) {
    dummy <- as.numeric(seq_len(n_group_classes - 1L) == j)
    softmax(c(intercepts + as.numeric(crossprod(slopes, dummy)), 0))
  }, numeric(n_profiles)))
  list(group_probabilities = softmax(c(group_logits, 0)),
       profile_probabilities = profile_probabilities)
}

#' Renormalize probabilities that a program printed to fixed precision
#'
#' Published output rounds each cell, so rows sum to one only up to the printed
#' precision. This checks the discrepancy is no larger than that rounding can
#' explain, then renormalizes, so that rounded published values can be fed back
#' in as an exact parameter set.
#'
#' @param block A matrix of probabilities with rows summing to one before rounding.
#' @param precision The printed resolution, for example 1e-3 for three decimals.
#' @return The matrix with rows renormalized to sum to one.
.renormalize_printed <- function(block, precision) {
  drift <- abs(rowSums(block) - 1)
  if (any(drift > ncol(block) * precision / 2)) {
    stop("Row sums depart from one by more than the printed rounding can explain.")
  }
  block / rowSums(block)
}

suite_mplus_twolevel_lca <- function() {
  artifacts <- file.path("validation", "mplus", "twolevel-lca")
  require_suite_files(
    file.path(artifacts, c("ex10.7.out.gz", "ex10.7.dat.gz")),
    reason = "the retained Mplus 8.8 run of User's Guide example 10.7 and its data")
  items <- paste0("u", seq_len(10))
  n_profiles <- 4L
  n_group_classes <- 5L
  n_categories <- 3L

  lines <- .mplus_output_lines(file.path(artifacts, "ex10.7.out.gz"))
  stopifnot(
    "Mplus did not terminate normally" =
      any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
    "Mplus did not replicate the best loglikelihood" =
      any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)),
    "Mplus reported an error" = !any(grepl("*** ERROR", lines, fixed = TRUE))
  )

  observed <- utils::read.table(gzfile(file.path(artifacts, "ex10.7.dat.gz")),
                                col.names = c(items, "dumb", "dumw", "clus"))
  fit <- multilpa(observed, vars = items, id = "clus",
                  n_profiles = n_profiles, n_group_classes = n_group_classes,
                  categorical = items, n_starts = 20L, seed = 20260919L,
                  tol = 1e-12, max_iter = 20000L)

  scale <- .mplus_probability_scale(lines, toupper(items), n_categories)
  # The measurement model is shared across group classes, so any one group
  # class carries it. Verify that before relying on it.
  first_group <- which(scale$pattern[, 1L] == 1L)
  profile_order_mplus <- scale$pattern[first_group, 2L]
  spread <- max(vapply(seq_len(n_profiles), function(s) {
    rows <- which(scale$pattern[, 2L] == s)
    max(apply(scale$probabilities[rows, , , drop = FALSE], c(2L, 3L),
              function(v) diff(range(v))))
  }, numeric(1)))
  stopifnot("Mplus did not hold the measurement model fixed across group classes" =
              spread < 1e-9)
  mplus_responses <- lapply(seq_along(items), function(item) {
    .renormalize_printed(
      matrix(scale$probabilities[first_group[order(profile_order_mplus)], item, ],
             n_profiles, n_categories), 1e-3)
  })
  prevalence <- .mplus_prevalence(lines, n_profiles, n_group_classes)

  # Score Mplus's own published estimates under the multilpa likelihood. This
  # separates agreement of the likelihood function from agreement of the
  # optimizer: the published parameters are rounded to three decimals, so
  # evaluating there must fall slightly below Mplus's reported maximum.
  at_mplus <- multilpa(observed, vars = items, id = "clus",
                       n_profiles = n_profiles, n_group_classes = n_group_classes,
                       categorical = items, n_starts = 1L, max_iter = 0L,
                       start = starting_values(list(
                         profile_probabilities = prevalence$profile_probabilities,
                         group_probabilities = prevalence$group_probabilities,
                         response_probabilities = mplus_responses)))

  profile_order <- align_classes(flatten_blocks(mplus_responses),
                                 flatten_blocks(fit$response_probabilities))
  responses <- lapply(fit$response_probabilities,
                      function(block) unname(as.matrix(block))[profile_order, , drop = FALSE])
  group_order <- align_classes(prevalence$profile_probabilities,
                               unname(fit$profile_probabilities)[, profile_order, drop = FALSE])
  profile_probabilities <- unname(fit$profile_probabilities)[group_order, profile_order, drop = FALSE]
  group_probabilities <- unname(fit$group_probabilities)[group_order]

  cell_names <- sprintf("item-response P(%s = %d | profile %d)",
                        rep(items, each = n_profiles * n_categories),
                        rep(rep(seq_len(n_categories), each = n_profiles), length(items)),
                        rep(seq_len(n_profiles), n_categories * length(items)))

  # These four are compared on a relative scale. Mplus prints them to three
  # decimals, but multilpa does not reproduce those digits: it converges a
  # little past Mplus's reported maximum, to -98736.121 against -98736.127. On a
  # log-likelihood of magnitude 98,736 an absolute gap of 0.006 is the wrong
  # thing to report, and "reproduces the printed digits" is the wrong question;
  # what matters is that the two agree to within 1e-7 of their own size.
  rbind(
    compare_values("maximised log-likelihood",
                   .mplus_scalar(lines, "H0 Value"), as.numeric(logLik(fit)),
                   tolerance = 1e-6, scale = "relative"),
    compare_values("free parameters",
                   .mplus_scalar(lines, "Number of Free Parameters"),
                   fit$n_parameters, tolerance = 0.5),
    compare_values("Akaike (AIC)", .mplus_scalar(lines, "Akaike (AIC)"),
                   fit$aic, tolerance = 1e-6, scale = "relative"),
    compare_values("Bayesian (BIC)", .mplus_scalar(lines, "Bayesian (BIC)"),
                   fit$bic_individual, tolerance = 1e-6, scale = "relative"),
    compare_values("log-likelihood at Mplus's published estimates",
                   .mplus_scalar(lines, "H0 Value"), as.numeric(logLik(at_mplus)),
                   tolerance = 1e-6, scale = "relative"),
    # Prevalence is recovered from logits printed to three decimals whose
    # reported standard errors are 0.22 to 0.41, so this tolerance sits about
    # two orders of magnitude below the sampling uncertainty Mplus reports.
    compare_values(sprintf("group-class probability %d", seq_len(n_group_classes)),
                   prevalence$group_probabilities, group_probabilities,
                   tolerance = 5e-3),
    compare_values(sprintf("profile prevalence P(profile %d | group class %d)",
                           rep(seq_len(n_profiles), each = n_group_classes),
                           rep(seq_len(n_group_classes), n_profiles)),
                   as.numeric(prevalence$profile_probabilities),
                   as.numeric(profile_probabilities), tolerance = 5e-3),
    compare_values(cell_names,
                   as.numeric(vapply(mplus_responses, identity,
                                     matrix(0, n_profiles, n_categories))),
                   as.numeric(vapply(responses, identity,
                                     matrix(0, n_profiles, n_categories))),
                   tolerance = 1.5e-3)
  ) |> transform(source = "Mplus 8.8 (UG example 10.7)", dataset = "ex10.7",
                 precision = "printed")
}
