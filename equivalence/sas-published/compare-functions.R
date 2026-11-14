# Functions for compare.R: lca(), lpa(), enumerate_lpa() and
# mixture_regression() against SAS output printed in tutorials and in SAS's own
# documentation:
#
# - PROC LCA (Lanza, Dziak, Huang, Wagner and Collins, Methodology Center,
#   Penn State) as printed by UCLA OARC (the data analysis example "Latent
#   Class Analysis", lca1.dat; the textbook examples of McCutcheon (1987),
#   Latent Class Analysis, pages 16, 33 and 47) and by Kleinman and Horton,
#   SAS and R, Example 8.21 (the HELP study, help.csv).
# - PROC FMM of SAS/STAT 14.2 (User's Guide, chapter 40, Getting Started):
#   the galaxies velocities (Roeder 1990) and "Student"'s yeast cell counts.
#
# README.md lists the pages and PROVENANCE.txt the data files. The printed
# blocks below are copied from those pages verbatim.
#
# Tolerances, fixed before the comparison: half a unit of the last printed
# digit plus 1e-6 of the value. Parameter counts and degrees of freedom
# exactly.
#
# PROC LCA prints its criteria on the deviance scale: G2 = 2 (L_saturated -
# L), AIC = G2 + 2p, BIC = G2 + p log n, CAIC = G2 + p (log n + 1), adjusted
# BIC = G2 + p log((n + 2) / 24), degrees of freedom = cells - p - 1. latents
# has no whole-table G2; it is computed here from latents' log likelihood and
# the saturated log likelihood of the observed response patterns, and the
# criteria from that G2 and latents' parameter count. Its "Entropy" is the
# relative entropy 1 - S / (n log K). Its standard errors invert the Hessian
# of the log likelihood (users' guide 1.3.2, section 4.3), so latents'
# observed-information errors are compared unscaled.
#
# PROC FMM's mixing logits are GLogit(Prob) = log(p_k / p_K), the last
# printed component the reference; their errors are the delta method on
# latents' observed-information covariance. Its AICC is -2L + 2pn / (n - p - 1).
#
# Not compared: PROC LCA's raw entropy and the McCutcheon models with
# restrictions (pages 33 to 70: fixed or equated probabilities, which latents
# does not fit); PROC FMM's Pearson statistic and maximum gradients, and the
# unequal-variance grid above three components (that likelihood is unbounded:
# with four to seven components latents puts a component on one or two
# velocities, while PROC FMM printed its three-component value for each).

# permutations() and delta_terms() from the wider Latent GOLD kit.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

# Reference names carry no comma: known-differences.csv matches them exactly.
.sas_sources <- c(
  dae = "SAS PROC LCA (UCLA OARC data analysis example)",
  mccutcheon = "SAS PROC LCA (UCLA OARC textbook examples of McCutcheon 1987)",
  help = "SAS PROC LCA (SAS and R Example 8.21)",
  fmm = "SAS/STAT 14.2 PROC FMM (User's Guide Getting Started)")

#' One comparison row of a printed value
#' @param printed The value as SAS printed it (text); `digits = NA` compares
#'   exactly.
#' @noRd
.sas_row <- function(case, fn, source, quantity, printed, obtained, digits = NULL) {
  digits <- digits %||% nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- if (is.na(digits)) 0 else 0.5 * 10^(-digits) + 1e-6 * abs(value)
  difference <- abs(obtained - value)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = difference, tolerance = tolerance, n_values = 1L,
             kind = "commercial program", agrees = difference <= tolerance, fn = fn)
}

.sas_exact <- function(case, fn, source, quantity, printed, obtained) {
  .sas_row(case, fn, source, quantity, as.character(printed), as.numeric(obtained), digits = NA)
}

# ---------------------------------------------------------------------------
# PROC LCA

#' Parse PROC LCA's "Gamma estimates" and "Rho estimates" blocks
#'
#' @param gamma The gamma block: a "Class:" line, the estimates, optionally
#'   their errors in parentheses.
#' @param rho The rho block: "Response category j:" headings, "Class:" lines,
#'   "item : estimates" lines, each optionally followed by its errors.
#' @return One row per printed parameter: `parameter` ("gamma" or "rho"),
#'   `class`, `item`, `category`, `estimate` and `standard_error` (text, NA
#'   when not printed).
#' @noRd
.proc_lca_estimates <- function(gamma, rho) {
  lines_of <- function(text) {
    lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1L]])
    lines[nzchar(lines)]
  }
  values <- function(line) strsplit(trimws(gsub("[()]", "", line)), "[[:space:]]+")[[1L]]
  gamma_lines <- lines_of(gamma)
  gamma_lines <- gamma_lines[!grepl("^Class:", gamma_lines)]
  estimates <- values(gamma_lines[[1L]])
  errors <- if (length(gamma_lines) > 1L) values(gamma_lines[[2L]]) else
    rep(NA_character_, length(estimates))
  gamma_rows <- data.frame(parameter = "gamma", class = seq_along(estimates), item = NA_character_,
                           category = NA_integer_, estimate = estimates, standard_error = errors)
  rho_lines <- lines_of(rho)
  rho_lines <- rho_lines[!grepl("^Class:", rho_lines)]
  category <- cumsum(grepl("^Response category", rho_lines))
  item_line <- grepl(":", rho_lines, fixed = TRUE) & !grepl("^Response category", rho_lines)
  rho_rows <- do.call(rbind, lapply(which(item_line), function(i) {
    parts <- strsplit(rho_lines[[i]], ":", fixed = TRUE)[[1L]]
    estimates <- values(parts[[2L]])
    next_line <- if (i < length(rho_lines)) rho_lines[[i + 1L]] else ""
    errors <- if (startsWith(next_line, "(")) values(next_line) else
      rep(NA_character_, length(estimates))
    data.frame(parameter = "rho", class = seq_along(estimates), item = trimws(parts[[1L]]),
               category = category[[i]], estimate = estimates, standard_error = errors)
  }))
  rbind(gamma_rows, rho_rows)
}

#' Persons as rows, items as factors with PROC LCA's response codes as levels
#' @param table A data frame of response patterns (codes 1, 2, ...).
#' @param weight The frequency column, or NULL for one row per person.
#' @noRd
.proc_lca_persons <- function(table, items, categories, weight = NULL) {
  persons <- if (is.null(weight)) table[items] else
    table[rep(seq_len(nrow(table)), table[[weight]]), items, drop = FALSE]
  persons[] <- Map(function(x, k) factor(x, levels = seq_len(k)), persons, categories)
  rownames(persons) <- NULL
  persons
}

.proc_lca_fit <- function(persons, k) {
  lca(persons, names(persons), n_classes = k, n_starts = 50L, seed = 1L, tol = 1e-14,
      max_iter = 100000L)
}

#' PROC LCA's fit statistics from latents' log likelihood
#' @param printed Named character vector: `ll`, `g2`, `aic`, `bic`, `df`,
#'   optionally `caic`, `abic`, `entropy`.
#' @noRd
.proc_lca_fit_rows <- function(case, source, fit, persons, categories, printed) {
  n <- nrow(persons)
  log_likelihood <- as.numeric(stats::logLik(fit))
  p <- fit$n_parameters
  patterns <- table(do.call(paste, persons))
  saturated <- sum(patterns * log(patterns / n))
  g2 <- 2 * (saturated - log_likelihood)
  obtained <- c(ll = log_likelihood, g2 = g2, aic = g2 + 2 * p, bic = g2 + p * log(n),
                caic = g2 + p * (log(n) + 1), abic = g2 + p * log((n + 2) / 24),
                entropy = if (nrow(get_results(fit, "profile_probabilities")) > 1L) {
                  entropy <- as.data.frame(get_results(fit, "entropy"))
                  entropy$relative_entropy[entropy$level == "individuals"]
                } else NA_real_)
  labels <- c(ll = "log likelihood", g2 = "G-squared", aic = "AIC", bic = "BIC", caic = "CAIC",
              abic = "adjusted BIC", entropy = "entropy")
  shown <- setdiff(names(printed), "df")
  rbind(do.call(rbind, lapply(shown, function(name)
          .sas_row(case, "lca", source, labels[[name]], printed[[name]], obtained[[name]]))),
        .sas_exact(case, "lca", source, "degrees of freedom", printed[["df"]],
                   prod(categories) - p - 1))
}

#' latents' class shares and response probabilities with observed-information
#' errors
#' @noRd
.proc_lca_latents <- function(fit) {
  inference <- as.data.frame(parameter_inference(fit, vcov_type = "observed", boundary = "fix"))
  responses <- as.data.frame(get_results(fit, "responses"))
  measurement <- inference[inference$level == "measurement", ]
  at <- match(paste0("profile_", responses$profile, " ", responses$indicator, ":",
                     responses$category),
              paste(measurement$outcome, measurement$term))
  shares <- as.data.frame(get_results(fit, "profile_probabilities"))
  profile <- inference[inference$level == "profile", ]
  share_at <- match(paste0("profile_", shares$profile), profile$outcome)
  # A probability held at its bound is fixed: no row, standard error 0.
  rbind(data.frame(parameter = "gamma", profile = shares$profile, item = NA_character_,
                   category = NA_integer_, estimate = shares$probability,
                   standard_error = ifelse(is.na(share_at), 0, profile$standard_error[share_at])),
        data.frame(parameter = "rho", profile = responses$profile, item = responses$indicator,
                   category = as.integer(as.character(responses$category)),
                   estimate = responses$probability,
                   standard_error = ifelse(is.na(at), 0, measurement$standard_error[at])))
}

#' Rows for every printed gamma and rho (and error, when printed)
#' @noRd
.proc_lca_estimate_rows <- function(case, source, fit, gamma, rho) {
  printed <- .proc_lca_estimates(gamma, rho)
  latents <- .proc_lca_latents(fit)
  key <- function(profile, rows) paste(rows$parameter, profile, rows$item, rows$category)
  k <- max(printed$class)
  # PROC LCA's class c is latents' profile perm[c]: the order closest in rho.
  orders <- permutations(k)
  gaps <- vapply(orders, function(perm) {
    at <- match(key(perm[printed$class], printed), key(latents$profile, latents))
    max(abs(latents$estimate[at] - as.numeric(printed$estimate)))
  }, numeric(1))
  perm <- orders[[which.min(gaps)]]
  at <- match(key(perm[printed$class], printed), key(latents$profile, latents))
  stopifnot("every printed parameter has a latents counterpart" = !anyNA(at))
  label <- ifelse(printed$parameter == "gamma", sprintf("gamma class %d", printed$class),
                  sprintf("rho class %d %s category %d", printed$class, printed$item,
                          printed$category))
  estimates <- do.call(rbind, lapply(seq_len(nrow(printed)), function(i)
    .sas_row(case, "lca", source, label[[i]], printed$estimate[[i]], latents$estimate[[at[[i]]]])))
  with_error <- which(!is.na(printed$standard_error))
  errors <- do.call(rbind, lapply(with_error, function(i)
    .sas_row(case, "lca", source, paste("s.e.", label[[i]]), printed$standard_error[[i]],
             latents$standard_error[[at[[i]]]])))
  rbind(estimates, errors)
}

# PROC LCA's printed output, case by case: fit statistics (text), the gamma
# block and the rho block, copied verbatim from the pages (README.md).
.proc_lca_printed <- list(
  proclca_drinking_k3 = list(
    fit = c(ll = "-4231.70", g2 = "319.96", aic = "377.96", bic = "520.28", caic = "549.28",
            abic = "428.17", entropy = "0.55", df = "482"),
    gamma = "
Class:                     1          2          3
                      0.5573     0.0793     0.3634
                     (0.1752)   (0.0243)   (0.1839)",
    rho = "
  Response category  1:
Class:                     1          2          3
  item1       :       0.9084     0.9233     0.3125
                     (0.0878)   (0.0411)   (0.2117)
  item2       :       0.3375     0.5461     0.1641
                     (0.0569)   (0.0782)   (0.0415)
  item3       :       0.0668     0.4264     0.0357
                     (0.0175)   (0.0938)   (0.0159)
  item4       :       0.0655     0.4180     0.0560
                     (0.0190)   (0.0817)   (0.0168)
  item5       :       0.2192     0.7655     0.0444
                     (0.0532)   (0.1000)   (0.0371)
  item6       :       0.3197     0.4710     0.1830
                     (0.0432)   (0.0763)   (0.0453)
  item7       :       0.1128     0.5124     0.0977
                     (0.0213)   (0.0929)   (0.0220)
  item8       :       0.1399     0.6193     0.1099
                     (0.0272)   (0.1081)   (0.0224)
  item9       :       0.3249     0.3488     0.1878
                     (0.0400)   (0.0698)   (0.0402)
  Response category  2:
Class:                     1          2          3
  item1       :       0.0916     0.0767     0.6875
                     (0.0878)   (0.0411)   (0.2117)
  item2       :       0.6625     0.4539     0.8359
                     (0.0569)   (0.0782)   (0.0415)
  item3       :       0.9332     0.5736     0.9643
                     (0.0175)   (0.0938)   (0.0159)
  item4       :       0.9345     0.5820     0.9440
                     (0.0190)   (0.0817)   (0.0168)
  item5       :       0.7808     0.2345     0.9556
                     (0.0532)   (0.1000)   (0.0371)
  item6       :       0.6803     0.5290     0.8170
                     (0.0432)   (0.0763)   (0.0453)
  item7       :       0.8872     0.4876     0.9023
                     (0.0213)   (0.0929)   (0.0220)
  item8       :       0.8601     0.3807     0.8901
                     (0.0272)   (0.1081)   (0.0224)
  item9       :       0.6751     0.6512     0.8122
                     (0.0400)   (0.0698)   (0.0402)"),
  mccutcheon_p16_k2 = list(
    fit = c(ll = "-565.59", g2 = "0.00", aic = "14.00", bic = "39.93", df = "0"),
    gamma = "
Class:                     1          2
                      0.5002     0.4998",
    rho = "
  Response category  1:
Class:                     1          2
  a           :       0.3330     0.6671
  b           :       0.1995     0.7006
  c           :       0.0013     0.9990
  Response category  2:
Class:                     1          2
  a           :       0.6670     0.3329
  b           :       0.8005     0.2994
  c           :       0.9987     0.0010"),
  mccutcheon_p33_k1 = list(
    fit = c(ll = "-2872.23", g2 = "257.26", aic = "269.26", bic = "299.81", df = "29"),
    gamma = "
Class:                     1
                      1.0000",
    rho = "
  Response category  1:
Class:                     1
  p           :       0.7646
  a           :       0.5200
  u           :       0.8153
  c           :       0.8386
  Response category  2:
Class:                     1
  p           :       0.0865
  a           :       0.4800
  u           :       0.1847
  c           :       0.1323
  Response category  3:
Class:                     1
  p           :       0.1489
  c           :       0.0291"),
  mccutcheon_p33_k2 = list(
    fit = c(ll = "-2783.27", g2 = "79.34", aic = "105.34", bic = "171.53", df = "22"),
    gamma = "
Class:                     1          2
                      0.8077     0.1923",
    rho = "
  Response category  1:
Class:                     1          2
  p           :       0.8953     0.2155
  a           :       0.6367     0.0298
  u           :       0.8327     0.7422
  c           :       0.8840     0.6478
  Response category  2:
Class:                     1          2
  p           :       0.0579     0.2066
  a           :       0.3633     0.9702
  u           :       0.1673     0.2578
  c           :       0.1043     0.2498
  Response category  3:
Class:                     1          2
  p           :       0.0468     0.5780
  c           :       0.0117     0.1024"),
  mccutcheon_p47_k1 = list(
    fit = c(ll = "-2346.67", g2 = "262.27", aic = "270.27", bic = "291.25", df = "11"),
    gamma = "
Class:                     1
                      1.0000",
    rho = "
  Response category  1:
Class:                     1
  w           :       0.9643
  a           :       0.9244
  i           :       0.6377
  v           :       0.2846
  Response category  2:
Class:                     1
  w           :       0.0357
  a           :       0.0756
  i           :       0.3623
  v           :       0.7154"),
  proclca_help_k3 = list(
    fit = c(ll = "-1032.48", g2 = "1.22", aic = "29.22", bic = "86.15", caic = "100.15",
            abic = "41.72", entropy = "0.94", df = "1"),
    gamma = "
Class:                     1          2          3
                      0.2163     0.0785     0.7052",
    rho = "
  Response category  1:
Class:                     1          2          3
  homeless    :       0.2703     1.0000     0.5625
  cesdcut     :       0.1154     0.4214     0.1678
  satreat     :       0.0004     0.0000     1.0000
  linkstatus  :       0.6029     1.0000     0.5855
  Response category  2:
Class:                     1          2          3
  homeless    :       0.7297     0.0000     0.4375
  cesdcut     :       0.8846     0.5786     0.8322
  satreat     :       0.9996     1.0000     0.0000
  linkstatus  :       0.3971     0.0000     0.4145"))

#' The persons of each PROC LCA case, items as factors with PROC LCA's codes
#'
#' UCLA's data analysis example: lca1.dat recoded by UCLA's DATA step
#' (it(i) = 2 - it(i), so 1 = yes). McCutcheon (1987) on UCLA's pages: page 16
#' (three items), page 33 (Table 3.2: four items with 3, 2, 2, 3 categories)
#' and page 47 (Table 4.2: four items), UCLA's SAS data sets written to CSV
#' unchanged (PROVENANCE.txt). SAS and R, Example 8.21: the book's DATA step
#' adds 1 to homeless, satreat and linkstatus and sets cesdcut = (cesd > 20) +
#' 1; PROC LCA analysed 431 subjects, the 453 of help.csv with linkstatus
#' observed (poLCA's listwise deletion in the same example gives the same log
#' likelihood).
#' @return A list of data frames named by case data set: `drinking`,
#'   `page16`, `page33`, `page47`, `help`.
#' @noRd
proc_lca_data <- function(directory) {
  path <- function(file) file.path(directory, "data", file)
  items <- paste0("item", 1:9)
  drinking <- utils::read.csv(path("lca1.dat"), header = FALSE, col.names = c("id", items))
  drinking[items] <- lapply(drinking[items], function(x) 2L - x)
  page <- function(name) utils::read.csv(path(sprintf("mccutcheon-%s.csv", name)))
  help <- utils::read.csv(path("help.csv"))
  help <- data.frame(homeless = help$homeless + 1L, cesdcut = as.integer(help$cesd > 20) + 1L,
                     satreat = help$satreat + 1L, linkstatus = help$linkstatus + 1L)
  help <- help[stats::complete.cases(help), ]
  data <- list(
    drinking = .proc_lca_persons(drinking, items, rep(2L, 9L)),
    page16 = .proc_lca_persons(page("page16"), c("a", "b", "c"), c(2L, 2L, 2L), "w"),
    page33 = .proc_lca_persons(page("page33"), c("p", "a", "u", "c"), c(3L, 2L, 2L, 3L), "wt"),
    page47 = .proc_lca_persons(page("page47"), c("w", "a", "i", "v"), rep(2L, 4L), "wt"),
    help = .proc_lca_persons(help, names(help), rep(2L, 4L)))
  stopifnot("1000 drinking subjects, 300 on page 16, 431 HELP subjects" =
              nrow(data$drinking) == 1000L && nrow(data$page16) == 300L &&
              nrow(data$help) == 431L)
  data
}

# Each PROC LCA case: its data set, number of classes and source.
.proc_lca_cases <- data.frame(
  case = names(.proc_lca_printed),
  data = c("drinking", "page16", "page33", "page33", "page47", "help"),
  k = c(3L, 2L, 1L, 2L, 1L, 3L),
  source = .sas_sources[c("dae", "mccutcheon", "mccutcheon", "mccutcheon", "mccutcheon", "help")],
  row.names = NULL)

#' Every PROC LCA case: fit statistics, gammas, rhos and printed errors
#' @noRd
compare_proc_lca <- function(directory) {
  data <- proc_lca_data(directory)
  do.call(rbind, lapply(seq_len(nrow(.proc_lca_cases)), function(i) {
    case <- .proc_lca_cases$case[[i]]
    persons <- data[[.proc_lca_cases$data[[i]]]]
    printed <- .proc_lca_printed[[case]]
    fit <- .proc_lca_fit(persons, .proc_lca_cases$k[[i]])
    categories <- vapply(persons, nlevels, integer(1))
    source <- .proc_lca_cases$source[[i]]
    rbind(.proc_lca_fit_rows(case, source, fit, persons, categories, printed$fit),
          .proc_lca_estimate_rows(case, source, fit, printed$gamma, printed$rho))
  }))
}

# ---------------------------------------------------------------------------
# PROC FMM

#' The galaxies velocities as the manual's DATA step reads them (v = velocity
#' / 1000)
#' @noRd
read_galaxies <- function(directory) {
  velocity <- scan(file.path(directory, "data", "galaxies.txt"), quiet = TRUE)
  stopifnot("82 galaxies" = length(velocity) == 82L)
  data.frame(v = velocity / 1000)
}

.fmm_fit <- function(data, k, model) {
  lpa(data, "v", n_profiles = k, model = model, n_starts = 50L, seed = 1L, tol = 1e-14,
      max_iter = 100000L)
}

#' A printed PROC FMM normal mixture: criteria, means, variances, mixing
#' probabilities and logits
#' @param printed List: `criteria` (named text: `m2ll`, `aic`, `aicc`, `bic`,
#'   `npar`), `means` and `variances` (a list of `c(estimate, s.e.)` per
#'   printed component), `mixing` (a list of `c(probability, glogit, s.e.)`,
#'   the last `c(probability)`).
#' @noRd
.fmm_normal_rows <- function(case, fit, printed) {
  source <- .sas_sources[["fmm"]]
  n <- 82
  inference <- parameter_inference(fit, vcov_type = "observed")
  rows <- as.data.frame(inference)
  covariance <- attr(inference, "covariance")
  theta <- stats::setNames(rows$estimate, rownames(covariance)[seq_len(nrow(rows))])
  stopifnot("the covariance names the estimates" =
              identical(names(theta), rownames(covariance)))
  measurement <- rows[rows$level == "measurement", ]
  means <- measurement[measurement$parameter == "mean", ]
  variances <- measurement[measurement$parameter == "variance", ]
  # PROC FMM's component g is the latents profile with the closest mean.
  profile_of <- vapply(printed$means, function(mean) {
    means$outcome[[which.min(abs(means$estimate - as.numeric(mean[[1L]])))]]
  }, character(1))
  stopifnot("each printed component is a different profile" = !anyDuplicated(profile_of))
  k <- length(profile_of)
  log_likelihood <- as.numeric(stats::logLik(fit))
  p <- fit$n_parameters
  criteria <- printed$criteria
  obtained <- c(m2ll = -2 * log_likelihood, aic = -2 * log_likelihood + 2 * p,
                aicc = -2 * log_likelihood + 2 * p * n / (n - p - 1),
                bic = -2 * log_likelihood + p * log(n))
  labels <- c(m2ll = "-2 log likelihood", aic = "AIC", aicc = "AICC", bic = "BIC")
  criterion_rows <- rbind(
    do.call(rbind, lapply(names(labels), function(name)
      .sas_row(case, "lpa", source, labels[[name]], criteria[[name]], obtained[[name]]))),
    .sas_exact(case, "lpa", source, "number of parameters", criteria[["npar"]], p))
  parameter_rows <- do.call(rbind, lapply(seq_len(k), function(g) {
    mean <- means[means$outcome == profile_of[[g]], ]
    # An equal-variance model has one variance, outcome "shared".
    variance <- variances[variances$outcome %in% c(profile_of[[g]], "shared"), ]
    rbind(.sas_row(case, "lpa", source, sprintf("component %d mean", g), printed$means[[g]][[1L]],
                   mean$estimate),
          .sas_row(case, "lpa", source, sprintf("s.e. component %d mean", g),
                   printed$means[[g]][[2L]], mean$standard_error),
          .sas_row(case, "lpa", source, sprintf("component %d variance", g),
                   printed$variances[[g]][[1L]], variance$estimate),
          .sas_row(case, "lpa", source, sprintf("s.e. component %d variance", g),
                   printed$variances[[g]][[2L]], variance$standard_error))
  }))
  share <- function(profile) sprintf("profile.probability.%s.group_class_1", profile)
  logits <- delta_terms(function(values) {
    stats::setNames(vapply(seq_len(k - 1L), function(g)
      log(values[[share(profile_of[[g]])]] / values[[share(profile_of[[k]])]]), numeric(1)),
      sprintf("glogit %d", seq_len(k - 1L)))
  }, theta, covariance)
  mixing_rows <- do.call(rbind, lapply(seq_len(k), function(g) {
    probability <- .sas_row(case, "lpa", source, sprintf("component %d mixing probability", g),
                            printed$mixing[[g]][[1L]], theta[[share(profile_of[[g]])]])
    if (g == k) return(probability)
    rbind(probability,
          .sas_row(case, "lpa", source, sprintf("component %d GLogit", g),
                   printed$mixing[[g]][[2L]], logits$estimate[[g]]),
          .sas_row(case, "lpa", source, sprintf("s.e. component %d GLogit", g),
                   printed$mixing[[g]][[3L]], logits$standard_error[[g]]))
  }))
  rbind(criterion_rows, parameter_rows, mixing_rows)
}

#' The galaxies: Figures 40.14 and 40.17, the best model of each grid in full
#' and, for equal variances, every grid model's -2 log likelihood and
#' parameter count
#' @noRd
compare_fmm_galaxies <- function(directory) {
  source <- .sas_sources[["fmm"]]
  data <- read_galaxies(directory)
  varying <- .fmm_fit(data, 3L, "VII")
  equal <- .fmm_fit(data, 4L, "EII")
  grid_equal <- enumerate_lpa(data, "v", n_profiles = 3:7, model = "EII", n_starts = 50L,
                              seed = 1L, tol = 1e-14, max_iter = 100000L)
  equal_m2ll <- c(`3` = "478.74", `4` = "416.49", `5` = "416.49", `6` = "416.49", `7` = "416.49")
  grid_rows <- do.call(rbind, lapply(3:7, function(k) {
    fit_equal <- .grid_candidate(grid_equal, k)
    rbind(.sas_row(sprintf("fmm_galaxies_equal_k%d", k), "enumeration", source,
                   "-2 log likelihood", equal_m2ll[[as.character(k)]],
                   -2 * as.numeric(stats::logLik(fit_equal))),
          .sas_exact(sprintf("fmm_galaxies_equal_k%d", k), "enumeration", source,
                     "number of parameters", 2L * k, fit_equal$n_parameters))
  }))
  rbind(
    .fmm_normal_rows("fmm_galaxies_varying_k3", varying, list(
      criteria = c(m2ll = "406.96", aic = "422.96", aicc = "424.94", bic = "442.22", npar = "8"),
      means = list(c("9.7101", "0.1597"), c("33.0444", "0.5322"), c("21.4039", "0.2597")),
      variances = list(c("0.1785", "0.09542"), c("0.8496", "0.6937"), c("4.8567", "0.8098")),
      mixing = list(c("0.0854", "-2.3308", "0.3959"), c("0.0366", "-3.1781", "0.5893"),
                    c("0.8781")))),
    .fmm_normal_rows("fmm_galaxies_equal_k4", equal, list(
      criteria = c(m2ll = "416.49", aic = "432.49", aicc = "434.47", bic = "451.75", npar = "8"),
      means = list(c("23.5058", "0.3460"), c("33.0440", "0.7610"), c("20.0086", "0.3029"),
                   c("9.7103", "0.4981")),
      variances = rep(list(c("1.7354", "0.3905")), 4L),
      mixing = list(c("0.3503", "1.4118", "0.4497"), c("0.0366", "-0.8473", "0.6901"),
                    c("0.5277", "1.8216", "0.4205"), c("0.0854")))),
    grid_rows)
}

#' The fit of an lpa enumeration grid with k profiles
#' @noRd
.grid_candidate <- function(grid, k) {
  candidates <- as.data.frame(get_results(grid, "candidates"))
  name <- candidates$candidate[candidates$n_profiles == k]
  stopifnot("one candidate with k profiles" = length(name) == 1L)
  get_candidate(grid, name)
}

#' "Student"'s yeast cell counts: a two-component binomial(5) mixture
#' (Figures 40.2 and 40.3)
#'
#' The six cells of Table 40.1, 400 squares; PROC FMM's FREQ statement is
#' the same likelihood as one row per square.
#' @noRd
compare_fmm_yeast <- function(directory) {
  source <- .sas_sources[["fmm"]]
  cells <- utils::read.csv(file.path(directory, "data", "yeast.csv"))
  squares <- data.frame(count = rep(cells$count, cells$f), n = 5L)
  stopifnot("400 squares, 273 cells" = nrow(squares) == 400L && sum(squares$count) == 273L)
  fit <- mixture_regression(cbind(count, n - count) ~ 1, squares, n_classes = 2L,
                            family = "binomial", n_starts = 50L, seed = 1L, tol = 1e-14,
                            max_iter = 100000L)
  theta <- stats::coef(fit)
  errors <- sqrt(diag(stats::vcov(fit)))
  intercept <- function(class) sprintf("coefficient.class_%d.(Intercept)", class)
  # PROC FMM's component 1 is the one with the smaller success probability.
  order <- order(theta[c(intercept(1L), intercept(2L))])
  membership <- "membership.class_2.(Intercept)"
  # latents' membership logit is log(p2 / p1); PROC FMM's is log(p_1 / p_2)
  # of its own components.
  logit <- if (identical(order, c(1L, 2L))) -theta[[membership]] else theta[[membership]]
  share_1 <- stats::plogis(logit)
  case <- "fmm_yeast_k2"
  row <- function(quantity, printed, obtained)
    .sas_row(case, "mixture_regression", source, quantity, printed, obtained)
  assignments <- as.data.frame(get_results(fit, "assignments"))
  posterior_columns <- sprintf("probability_class_%d", order)
  stopifnot("the assignments hold each class's posterior" =
              all(posterior_columns %in% names(assignments)))
  first <- match(cells$count, squares$count)
  posterior_1 <- c("0.98606", "0.91089", "0.59638", "0.17598", "0.02994", "0.00444")
  posterior_2 <- c("0.01394", "0.08911", "0.40362", "0.82402", "0.97006", "0.99556")
  rbind(
    row("component 1 intercept", "-2.2316", theta[[intercept(order[[1L]])]]),
    row("s.e. component 1 intercept", "0.1522", errors[[intercept(order[[1L]])]]),
    row("component 1 success probability", "0.09695", stats::plogis(theta[[intercept(order[[1L]])]])),
    row("component 2 intercept", "-0.2974", theta[[intercept(order[[2L]])]]),
    row("s.e. component 2 intercept", "0.3655", errors[[intercept(order[[2L]])]]),
    row("component 2 success probability", "0.4262", stats::plogis(theta[[intercept(order[[2L]])]])),
    row("component 1 mixing probability", "0.8799", share_1),
    row("component 1 logit", "1.9913", logit),
    row("s.e. component 1 logit", "0.5725", errors[[membership]]),
    row("component 2 mixing probability", "0.1201", 1 - share_1),
    row("component 1 predicted probability", "0.096951", stats::plogis(theta[[intercept(order[[1L]])]])),
    row("component 2 predicted probability", "0.42620", stats::plogis(theta[[intercept(order[[2L]])]])),
    do.call(rbind, lapply(seq_len(nrow(cells)), function(i) rbind(
      row(sprintf("posterior component 1 count %d", cells$count[[i]]), posterior_1[[i]],
          assignments[[posterior_columns[[1L]]]][[first[[i]]]]),
      row(sprintf("posterior component 2 count %d", cells$count[[i]]), posterior_2[[i]],
          assignments[[posterior_columns[[2L]]]][[first[[i]]]])))))
}

compare_all_sas <- function(directory) {
  rbind(compare_proc_lca(directory), compare_fmm_galaxies(directory),
        compare_fmm_yeast(directory))
}
