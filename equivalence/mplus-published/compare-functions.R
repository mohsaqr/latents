# Functions for compare.R: lca(), lpa(), lta() and the Vuong-Lo-Mendell-Rubin
# statistics against genuine Mplus output printed on the UCLA OARC Mplus
# pages (Mplus textbook examples of Applied Latent Class Analysis, Hagenaars
# and McCutcheon (eds.) 2002, chapters 2, 3 and 11; of Dayton (1998) Latent
# Class Scaling Analysis; and the Mplus data analysis example "Latent Class
# Analysis"). The printed blocks below are copied from those pages verbatim;
# README.md lists the pages and PROVENANCE.txt the data files.
#
# Tolerances, fixed before the comparison: half a unit of the last printed
# digit plus 1e-6 of the value. Parameter counts and modal class counts
# exactly.
#
# Standard errors: Mplus estimates mixtures by MLR by default, a sandwich
# without a small-sample factor (CR0); one page (Dayton's cheating example)
# printed ML's, from the observed information (see compare_dayton()). latents' robust covariance carries the CR1
# factor n / (n - 1) of sandwich::vcovCL() (.inference_cluster_factor() in
# R/inference-service.R), so its standard errors are multiplied by
# sqrt((n - 1) / n) before the comparison. Mplus holds a threshold that
# reaches +-15 fixed; parameter_inference(boundary = "fix") holds the same
# probability at its bound, and such a parameter's standard error is 0.
#
# Not compared: Mplus's TECH11 p-values, mean and standard deviation (its H0
# is re-estimated from the K-class model with its first class deleted and its
# weights are built differently; see test-features-vlmr.R), the chi-square
# fit tests (latents has none for a whole table), and Mplus's entropy of a
# model with several latent class variables (computed over their joint
# pattern; latents' entropy of a transition fit is over occasions).

# permutations() from the wider Latent GOLD kit.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

.ucla_pages <- c(
  ch2 = "UCLA OARC, Mplus textbook examples: McCutcheon (2002), Applied Latent Class Analysis, ch. 2",
  ch3 = "UCLA OARC, Mplus textbook examples: Vermunt and Magidson (2002), Applied Latent Class Analysis, ch. 3, Table 1",
  ch11 = "UCLA OARC, Mplus textbook examples: Langeheine and van de Pol (2002), Applied Latent Class Analysis, ch. 11, Table 3",
  dayton = "UCLA OARC, Mplus textbook examples: Dayton (1998), Latent Class Scaling Analysis",
  dae = "UCLA OARC, Mplus data analysis examples: Latent Class Analysis")

#' One comparison row of a printed value
#' @param printed The value as Mplus printed it (text), or a number with
#'   `digits` given.
#' @noRd
.mplus_row <- function(case, fn, page, quantity, printed, obtained, digits = NULL) {
  digits <- digits %||% nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- if (is.na(digits)) 0 else 0.5 * 10^(-digits) + 1e-6 * abs(value)
  difference <- abs(obtained - value)
  data.frame(case = case, reference = "Mplus", printed_in = page, quantity = quantity,
             reference_value = value, latents = obtained, difference = difference,
             tolerance = tolerance, n_values = 1L, kind = "commercial program",
             agrees = difference <= tolerance, fn = fn)
}

#' An exact comparison (a parameter count or a modal count)
#' @noRd
.mplus_exact <- function(case, fn, page, quantity, printed, obtained) {
  .mplus_row(case, fn, page, quantity, as.character(printed), as.numeric(obtained), digits = NA)
}

#' The fit statistics Mplus prints: log likelihood, parameters, BIC, entropy
#' @param printed Named character vector: `ll`, `npar`, `bic`, optionally
#'   `entropy`.
#' @noRd
.mplus_fit_rows <- function(case, fn, page, fit, printed) {
  entropy <- if (!is.na(printed["entropy"])) {
    table <- get_results(fit, "entropy")
    .mplus_row(case, fn, page, "entropy", printed[["entropy"]],
               table$relative_entropy[table$level == "individuals"])
  }
  rbind(.mplus_row(case, fn, page, "log likelihood", printed[["ll"]], as.numeric(stats::logLik(fit))),
        .mplus_exact(case, fn, page, "number of parameters", printed[["npar"]], fit$n_parameters),
        .mplus_row(case, fn, page, "BIC", printed[["bic"]], stats::BIC(fit)),
        entropy)
}

#' Parse Mplus's RESULTS IN PROBABILITY SCALE block
#' @param text The block as printed ("Latent Class k", an indicator name,
#'   "Category j  estimate  s.e.  est./s.e." lines).
#' @return One row per class, indicator and category: printed estimate and
#'   standard error (text).
#' @noRd
.mplus_probability_block <- function(text) {
  lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  class <- cumsum(grepl("^Latent Class [0-9]+$", lines))
  indicator_line <- !grepl("^(Latent Class|Category)", lines)
  indicator <- cumsum(indicator_line)
  names_seen <- lines[indicator_line]
  rows <- grepl("^Category", lines)
  fields <- strsplit(lines[rows], "[[:space:]]+")
  data.frame(class = class[rows], indicator = tolower(names_seen[indicator[rows]]),
             category = vapply(fields, `[[`, "", 2L),
             estimate = vapply(fields, `[[`, "", 3L), standard_error = vapply(fields, `[[`, "", 4L),
             stringsAsFactors = FALSE)
}

#' latents' response probabilities with standard errors of Mplus's estimator
#' @param estimator "MLR" (Mplus's default for mixtures: the robust
#'   covariance without the CR1 factor) or "ML" (observed information).
#' @return One row per profile, indicator and category.
#' @noRd
.latents_probabilities <- function(fit, n, estimator = c("MLR", "ML")) {
  estimator <- match.arg(estimator)
  inference <- parameter_inference(fit, vcov_type = if (estimator == "MLR") "robust" else "observed",
                                   boundary = "fix")
  responses <- get_results(fit, "responses")
  term <- paste0(responses$indicator, ":", responses$category)
  at <- match(paste(paste0("profile_", responses$profile), term),
              paste(inference$outcome, inference$term)[inference$level == "measurement"])
  measurement <- inference[inference$level == "measurement", ]
  # A probability held at its bound is fixed: no row, standard error 0.
  se <- ifelse(is.na(at), 0, measurement$standard_error[at]) *
    if (estimator == "MLR") sqrt((n - 1) / n) else 1
  data.frame(profile = responses$profile, indicator = responses$indicator,
             category = as.character(responses$category), probability = responses$probability,
             standard_error = se, stringsAsFactors = FALSE)
}

#' The class order matching Mplus's printed probabilities
#' @return `perm`, where Mplus class c is latents profile perm[c].
#' @noRd
.mplus_class_order <- function(printed, latents, categories) {
  k <- max(printed$class)
  orders <- permutations(k)
  gap <- vapply(orders, function(perm) {
    obtained <- latents$probability[match(
      paste(perm[printed$class], printed$indicator, categories[printed$category_position]),
      paste(latents$profile, latents$indicator, latents$category))]
    max(abs(obtained - as.numeric(printed$estimate)))
  }, numeric(1))
  orders[[which.min(gap)]]
}

#' Rows for every printed probability and its standard error
#' @param categories latents' category labels in Mplus's order (Mplus's
#'   Category 1 is the lowest value).
#' @noRd
.mplus_probability_rows <- function(case, fn, page, fit, n, text, categories,
                                    estimator = "MLR") {
  printed <- .mplus_probability_block(text)
  printed$category_position <- as.integer(printed$category)
  latents <- .latents_probabilities(fit, n, estimator)
  perm <- .mplus_class_order(printed, latents, categories)
  at <- match(paste(perm[printed$class], printed$indicator, categories[printed$category_position]),
              paste(latents$profile, latents$indicator, latents$category))
  stopifnot("every printed probability has a latents counterpart" = !anyNA(at))
  label <- sprintf("class %d %s category %d", printed$class, printed$indicator,
                   printed$category_position)
  estimates <- do.call(rbind, lapply(seq_len(nrow(printed)), function(i)
    .mplus_row(case, fn, page, paste("probability", label[[i]]), printed$estimate[[i]],
               latents$probability[[at[[i]]]])))
  errors <- do.call(rbind, lapply(seq_len(nrow(printed)), function(i)
    .mplus_row(case, fn, page, paste("s.e. probability", label[[i]]), printed$standard_error[[i]],
               latents$standard_error[[at[[i]]]])))
  structure(rbind(estimates, errors), perm = perm)
}

#' Rows for Mplus's "FINAL CLASS COUNTS ... BASED ON THE ESTIMATED MODEL"
#' @noRd
.mplus_count_rows <- function(case, fn, page, fit, n, perm, printed) {
  shares <- get_results(fit, "profile_probabilities")
  do.call(rbind, lapply(seq_along(printed), function(c)
    .mplus_row(case, fn, page, sprintf("class %d count", c), printed[[c]],
               n * shares$probability[shares$profile == perm[[c]]])))
}

#' A frequency table expanded to one row per person, indicators as factors
#' @noRd
.expand_frequencies <- function(file, names, items, weight, keep = NULL) {
  table <- utils::read.csv(file, header = FALSE, col.names = names)
  if (!is.null(keep)) table <- table[keep(table), , drop = FALSE]
  persons <- table[rep(seq_len(nrow(table)), table[[weight]]), items, drop = FALSE]
  persons[] <- lapply(persons, factor)
  rownames(persons) <- NULL
  persons
}

.fit_lca <- function(data, k) {
  lca(data, names(data), n_classes = k, n_starts = 50L, seed = 1L, tol = 1e-14, max_iter = 100000L)
}

# ---------------------------------------------------------------------------
# McCutcheon (2002), ALCA ch. 2: Ego's Dilemma (Stouffer and Toby) and
# abortion approval

#' Ego's Dilemma (page59_a.dat): one and two classes, Tables 2 to 4
#'
#' The file holds two groups of 216 (the "group" column, used by Table 10's
#' two-group models); Tables 2 to 4 are of group 1 alone: Mplus prints 216
#' persons (155.68287 + 60.31713), and the log likelihood of group 1 is the
#' printed one (Stouffer and Toby's data, as in Stata's Example 50g).
#' @noRd
compare_mccutcheon_ego <- function(directory) {
  page <- .ucla_pages[["ch2"]]
  data <- .expand_frequencies(file.path(directory, "data", "page59_a.dat"),
                              c("a", "b", "c", "d", "group", "freq"), c("a", "b", "c", "d"),
                              "freq", keep = function(x) x$group == 1L)
  n <- nrow(data)
  one <- .fit_lca(data, 1L)
  two <- .fit_lca(data, 2L)
  probabilities <- .mplus_probability_rows("mccutcheon_ego_k2", "lca", page, two, n, "
Latent Class 1
 A
    Category 1         0.286    0.042      6.841
    Category 2         0.714    0.042     17.045
 B
    Category 1         0.646    0.049     13.175
    Category 2         0.354    0.049      7.220
 C
    Category 1         0.670    0.051     13.140
    Category 2         0.330    0.051      6.461
 D
    Category 1         0.868    0.039     22.325
    Category 2         0.132    0.039      3.406
Latent Class 2
 A
    Category 1         0.007    0.025      0.269
    Category 2         0.993    0.025     39.267
 B
    Category 1         0.073    0.068      1.088
    Category 2         0.927    0.068     13.716
 C
    Category 1         0.060    0.067      0.896
    Category 2         0.940    0.067     13.985
 D
    Category 1         0.231    0.098      2.351
    Category 2         0.769    0.098      7.833", c("1", "2"))
  perm <- attr(probabilities, "perm")
  # Table 3 (MODEL RESULTS): thresholds logit P(category 1) and the class
  # logit log(p1 / p2), with delta-method errors from the probabilities.
  latents <- .latents_probabilities(two, n)
  thresholds <- list(
    `class 1 threshold a` = c("-0.913", "0.205"), `class 1 threshold b` = c("0.601", "0.214"),
    `class 1 threshold c` = c("0.710", "0.231"), `class 1 threshold d` = c("1.880", "0.338"),
    `class 2 threshold a` = c("-4.983", "3.741"), `class 2 threshold b` = c("-2.535", "0.992"),
    `class 2 threshold c` = c("-2.747", "1.187"), `class 2 threshold d` = c("-1.203", "0.553"))
  threshold_rows <- do.call(rbind, lapply(names(thresholds), function(name) {
    parts <- strsplit(name, " ", fixed = TRUE)[[1L]]
    row <- latents[latents$profile == perm[[as.integer(parts[[2L]])]] &
                     latents$indicator == parts[[4L]] & latents$category == "1", ]
    p <- row$probability
    rbind(.mplus_row("mccutcheon_ego_k2", "lca", page, paste("estimate", name),
                     thresholds[[name]][[1L]], stats::qlogis(p)),
          .mplus_row("mccutcheon_ego_k2", "lca", page, paste("s.e.", name),
                     thresholds[[name]][[2L]], row$standard_error / (p * (1 - p))))
  }))
  inference <- parameter_inference(two, vcov_type = "robust", boundary = "fix")
  share <- inference[inference$level == "profile", ]
  p1 <- share$estimate[share$outcome == sprintf("profile_%d", perm[[1L]])]
  se1 <- share$standard_error[share$outcome == sprintf("profile_%d", perm[[1L]])] *
    sqrt((n - 1) / n)
  rbind(.mplus_fit_rows("mccutcheon_ego_k1", "lca", page, one,
                        c(ll = "-543.650", npar = "4", bic = "1108.801", entropy = NA)),
        .mplus_fit_rows("mccutcheon_ego_k2", "lca", page, two,
                        c(ll = "-504.468", npar = "9", bic = "1057.313", entropy = "0.719")),
        .mplus_count_rows("mccutcheon_ego_k2", "lca", page, two, n, perm,
                          c("155.68287", "60.31713")),
        probabilities, threshold_rows,
        .mplus_row("mccutcheon_ego_k2", "lca", page, "estimate class logit 1", "0.948",
                   stats::qlogis(p1)),
        .mplus_row("mccutcheon_ego_k2", "lca", page, "s.e. class logit 1", "0.300",
                   se1 / (p1 * (1 - p1))))
}

#' Abortion approval (page75.dat): two classes, Table 8 model H1
#' @noRd
compare_mccutcheon_abortion <- function(directory) {
  data <- .expand_frequencies(file.path(directory, "data", "page75.dat"),
                              c("a", "b", "c", "d", "freq"), c("a", "b", "c", "d"), "freq")
  .mplus_fit_rows("mccutcheon_abortion_k2", "lca", .ucla_pages[["ch2"]], .fit_lca(data, 2L),
                  c(ll = "-2773.793", npar = "9", bic = "5614.323", entropy = "0.925"))
}

# ---------------------------------------------------------------------------
# Dayton (1998): pleural thickening, academic cheating, cheating by gender

#' One unconstrained two-class model of Dayton's: fit, counts, probabilities
#' @noRd
.compare_dayton <- function(case, file, items, fit_printed, counts, text, estimator = "MLR") {
  page <- .ucla_pages[["dayton"]]
  data <- .expand_frequencies(file, c(items, "freq"), items, "freq")
  n <- nrow(data)
  fit <- .fit_lca(data, 2L)
  # Mplus's Category 1 is the lowest value (0 in these files).
  categories <- levels(data[[items[[1L]]]])
  stopifnot("every indicator has the same two categories" =
              all(vapply(data, \(x) identical(levels(x), categories), logical(1))))
  probabilities <- .mplus_probability_rows(case, "lca", page, fit, n, text, categories, estimator)
  rbind(.mplus_fit_rows(case, "lca", page, fit, fit_printed),
        if (length(counts)) .mplus_count_rows(case, "lca", page, fit, n,
                                              attr(probabilities, "perm"), counts),
        probabilities)
}

compare_dayton <- function(directory) {
  path <- function(file) file.path(directory, "data", file)
  rbind(
    .compare_dayton("dayton_pleural_k2", path("table3_1.dat"), c("a", "b", "c"),
                    c(ll = "-891.143", npar = "7", bic = "1834.321", entropy = "0.949"),
                    c("92.18101", "1599.81899"), "
Latent Class 1
 A
    Category 1         0.251    0.067      3.727
    Category 2         0.749    0.067     11.127
 B
    Category 1         0.356    0.067      5.357
    Category 2         0.644    0.067      9.670
 C
    Category 1         0.235    0.067      3.531
    Category 2         0.765    0.067     11.495
Latent Class 2
 A
    Category 1         0.990    0.003    285.156
    Category 2         0.010    0.003      2.871
 B
    Category 1         0.965    0.005    194.842
    Category 2         0.035    0.005      7.157
 C
    Category 1         0.989    0.004    271.551
    Category 2         0.011    0.004      3.000"),
    # The printed input names no estimator, but this page's standard errors
    # are ML's (observed information): all sixteen agree with it to the
    # printed digit, four with MLR's. The other pages' are MLR's.
    .compare_dayton("dayton_cheating_k2", path("table3_4.dat"), c("a", "b", "c", "d"),
                    c(ll = "-440.027", npar = "9", bic = "931.941", entropy = "0.737"),
                    character(), estimator = "ML", text = "
Latent Class 1
 A
    Category 1         0.423    0.180      2.346
    Category 2         0.577    0.180      3.199
 B
    Category 1         0.411    0.175      2.350
    Category 2         0.589    0.175      3.369
 C
    Category 1         0.784    0.085      9.272
    Category 2         0.216    0.085      2.555
 D
    Category 1         0.624    0.100      6.218
    Category 2         0.376    0.100      3.752
Latent Class 2
 A
    Category 1         0.983    0.029     34.212
    Category 2         0.017    0.029      0.578
 B
    Category 1         0.971    0.030     31.849
    Category 2         0.029    0.030      0.959
 C
    Category 1         0.963    0.015     63.319
    Category 2         0.037    0.015      2.439
 D
    Category 1         0.818    0.026     31.148
    Category 2         0.182    0.026      6.928"),
    .compare_dayton("dayton_female_k2", path("table6_1_f-1.dat"), c("a", "b", "c", "d"),
                    c(ll = "-273.207", npar = "9", bic = "593.151", entropy = "0.761"),
                    c("26.27241", "153.72759"), "
Latent Class 1
 A
    Category 1         0.356    0.301      1.183
    Category 2         0.644    0.301      2.144
 B
    Category 1         0.307    0.259      1.185
    Category 2         0.693    0.259      2.671
 C
    Category 1         0.812    0.120      6.766
    Category 2         0.188    0.120      1.564
 D
    Category 1         0.625    0.147      4.241
    Category 2         0.375    0.147      2.544
Latent Class 2
 A
    Category 1         0.980    0.043     22.900
    Category 2         0.020    0.043      0.467
 B
    Category 1         0.936    0.053     17.775
    Category 2         0.064    0.053      1.210
 C
    Category 1         0.941    0.021     44.072
    Category 2         0.059    0.021      2.762
 D
    Category 1         0.791    0.035     22.358
    Category 2         0.209    0.035      5.912"),
    .compare_dayton("dayton_male_k2", path("table6_1_m-1.dat"), c("a", "b", "c", "d"),
                    c(ll = "-156.177", npar = "9", bic = "356.634", entropy = "0.797"),
                    c("117.08534", "19.91466"), "
Latent Class 1
 A
    Category 1         0.977    0.027     36.663
    Category 2         0.023    0.027      0.846
 B
    Category 1         1.000    0.000      0.000
    Category 2         0.000    0.000      0.000
 C
    Category 1         0.994    0.022     45.732
    Category 2         0.006    0.022      0.259
 D
    Category 1         0.855    0.039     21.915
    Category 2         0.145    0.039      3.710
Latent Class 2
 A
    Category 1         0.430    0.206      2.087
    Category 2         0.570    0.206      2.771
 B
    Category 1         0.548    0.197      2.776
    Category 2         0.452    0.197      2.289
 C
    Category 1         0.682    0.129      5.293
    Category 2         0.318    0.129      2.472
 D
    Category 1         0.546    0.148      3.686
    Category 2         0.454    0.148      3.069"))
}

# ---------------------------------------------------------------------------
# UCLA's data analysis example: nine drinking items, three classes, TECH11

#' The drinking example (lca1.dat): three classes and TECH11 for 2 against 3
#'
#' Printed: each class's probability of "yes" (Mplus's category 2, coded 1)
#' to three decimals, the modal and the posterior class counts, and TECH11's
#' H0 log likelihood, twice the log likelihood difference, the difference in
#' parameters and the adjusted LRT statistic.
#' @noRd
compare_drinking <- function(directory) {
  page <- .ucla_pages[["dae"]]
  items <- paste0("item", 1:9)
  data <- utils::read.csv(file.path(directory, "data", "lca1.dat"), header = FALSE,
                          col.names = c("id", items))
  data[items] <- lapply(data[items], factor)
  grid <- enumerate_lca(data, items, n_classes = 2:3, n_starts = 50L, seed = 1L, tol = 1e-14,
                        max_iter = 100000L)
  three <- get_candidate(grid, "classes3")
  two <- get_candidate(grid, "classes2")
  yes <- rbind(
    item1 = c("0.908", "0.312", "0.923"), item2 = c("0.337", "0.164", "0.546"),
    item3 = c("0.067", "0.036", "0.426"), item4 = c("0.065", "0.056", "0.418"),
    item5 = c("0.219", "0.044", "0.765"), item6 = c("0.320", "0.183", "0.471"),
    item7 = c("0.113", "0.098", "0.512"), item8 = c("0.140", "0.110", "0.619"),
    item9 = c("0.325", "0.188", "0.349"))
  responses <- get_results(three, "responses")
  responses <- responses[responses$category == "1", ]
  gap <- vapply(permutations(3L), function(perm) max(abs(vapply(seq_len(3L), function(c)
    responses$probability[match(paste(perm[[c]], items), paste(responses$profile,
                                                               responses$indicator))],
    numeric(9)) - matrix(as.numeric(yes), 9L))), numeric(1))
  perm <- permutations(3L)[[which.min(gap)]]
  probability_rows <- do.call(rbind, lapply(seq_len(3L), function(c) do.call(rbind, lapply(items,
    function(item) .mplus_row("ucla_drinking_k3", "lca", page,
                              sprintf("probability class %d %s yes", c, item), yes[item, c],
                              responses$probability[responses$profile == perm[[c]] &
                                                      responses$indicator == item])))))
  classification <- get_results(three, "classification")
  classification <- classification[classification$level == "individuals", ]
  modal <- c(646L, 288L, 66L)
  posterior <- c("557.56836", "363.13989", "79.29175")
  count_rows <- do.call(rbind, lapply(seq_len(3L), function(c) {
    row <- classification[classification$class == perm[[c]], ]
    rbind(.mplus_exact("ucla_drinking_k3", "lca", page, sprintf("class %d modal count", c),
                       modal[[c]], row$n_modal),
          .mplus_row("ucla_drinking_k3", "lca", page, sprintf("class %d posterior count", c),
                     posterior[[c]], row$estimated_n))
  }))
  lmr <- get_results(grid, "lmr")
  lmr <- lmr[lmr$null_candidate == "classes2" & lmr$alternative_candidate == "classes3", ]
  stopifnot("one VLMR row for 2 against 3 classes" = nrow(lmr) == 1L)
  case <- "ucla_drinking_vlmr_k2_k3"
  rbind(probability_rows, count_rows,
        .mplus_row(case, "bootstrap_lrt", page, "TECH11 H0 log likelihood", "-4251.208",
                   as.numeric(stats::logLik(two))),
        .mplus_row(case, "bootstrap_lrt", page, "TECH11 2 times the log likelihood difference",
                   "39.025", lmr$lr),
        .mplus_exact(case, "bootstrap_lrt", page, "TECH11 difference in the number of parameters",
                     10L, lmr$parameters_added),
        .mplus_row(case, "bootstrap_lrt", page, "TECH11 adjusted LRT value", "38.468",
                   lmr$adjusted_lr))
}

# ---------------------------------------------------------------------------
# Vermunt and Magidson (2002), ALCA ch. 3, Table 1: the diabetes data

#' Ten latent profile models of EXAM1.dat (Reaven and Miller's diabetes data,
#' the version mclust shipped before 6.0), Mplus's class-dependent and
#' class-independent, unrestricted and diagonal covariances at one to three
#' clusters (Models 1, 2, 3 and 4; the one-cluster models 2:1 and 4:1 are
#' 1:1 and 3:1)
#' @noRd
compare_diabetes <- function(directory) {
  data <- utils::read.table(file.path(directory, "data", "EXAM1.dat"), header = TRUE)
  table <- data.frame(
    model = c("VVV", "VVV", "VVV", "EEE", "EEE", "VVI", "VVI", "VVI", "EEI", "EEI"),
    k = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 3L, 2L, 3L),
    ll = c("-2545.828", "-2355.905", "-2303.492", "-2473.144", "-2417.341", "-2750.135",
           "-2445.312", "-2364.137", "-2559.816", "-2464.383"),
    npar = c("9", "19", "29", "13", "17", "6", "13", "20", "10", "14"),
    bic = c("5136.446", "4806.369", "4751.309", "5010.985", "4919.287", "5530.130",
            "4955.322", "4827.809", "5169.399", "4998.440"),
    entropy = c(NA, "0.935", "0.855", "0.980", "0.993", NA, "0.980", "0.902", "0.994", "0.994"),
    stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(table)), function(i) {
    row <- table[i, ]
    fit <- lpa(data, c("y1", "y2", "y3"), n_profiles = row$k, model = row$model, n_starts = 50L,
               seed = 1L, tol = 1e-14, max_iter = 200000L)
    .mplus_fit_rows(sprintf("vermunt_diabetes_%s_k%d", row$model, row$k), "lpa",
                    .ucla_pages[["ch3"]], fit,
                    c(ll = row$ll, npar = row$npar, bic = row$bic, entropy = row$entropy))
  }))
}

# ---------------------------------------------------------------------------
# Langeheine and van de Pol (2002), ALCA ch. 11: latent Markov chains

#' Models 7a and 7b of Table 3: two latent states measured by one binary
#' response at each of five occasions (chap11.dat; Mplus's input reads
#' a b c d e count and leaves out the file's group column, so the two groups'
#' counts are pooled: 5147 persons), measurement invariant, transitions
#' homogeneous (7a) or per occasion (7b)
#' @noRd
compare_latent_markov <- function(directory) {
  table <- utils::read.csv(file.path(directory, "data", "chap11.dat"), header = FALSE,
                           col.names = c("a", "b", "c", "d", "e", "count", "group"))
  occasions <- c("a", "b", "c", "d", "e")
  wide <- table[rep(seq_len(nrow(table)), table$count), occasions]
  wide$id <- seq_len(nrow(wide))
  long <- stats::reshape(wide, direction = "long", varying = occasions, v.names = "y",
                         timevar = "time", times = seq_along(occasions), idvar = "id")
  long$y <- factor(long$y)
  long <- long[order(long$id, long$time), ]
  page <- .ucla_pages[["ch11"]]
  fit <- function(transitions) {
    lta(long, "y", id = "id", n_profiles = 2L, time = "time", categorical = "y",
        transitions = transitions, n_starts = 20L, seed = 1L, tol = 1e-12)
  }
  rbind(.mplus_fit_rows("langeheine_lm_homogeneous_k2", "lta", page, fit("homogeneous"),
                        c(ll = "-15378.476", npar = "5", bic = "30799.683", entropy = NA)),
        # Model 7b's fit has a transition at zero: lta() warns that it is a
        # boundary fit (its errors are withheld; none is compared here).
        .mplus_fit_rows("langeheine_lm_occasion_k2", "lta", page, fit("occasion"),
                        c(ll = "-15325.635", npar = "11", bic = "30745.279", entropy = NA)))
}

compare_all_mplus <- function(directory) {
  rbind(compare_mccutcheon_ego(directory), compare_mccutcheon_abortion(directory),
        compare_dayton(directory), compare_drinking(directory), compare_diabetes(directory),
        compare_latent_markov(directory))
}
