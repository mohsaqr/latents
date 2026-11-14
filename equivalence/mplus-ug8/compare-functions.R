# Functions for compare.R: the growth and latent transition examples of the
# Mplus 8.8 User's Guide, chapter 8 (statmodel.com, chapter8.zip): input,
# data and Mplus's own output are kept in listings/ (LISTINGS.md5), and every
# compared value is read from the output file, not copied by hand.
#
# - ex8.9, a latent class growth analysis of a binary outcome at four
#   occasions: mixture_regression(family = "binomial", class_level = "group").
#   Mplus's logit P(u = 1) = -threshold + I + S * time, with one threshold for
#   every occasion and class and I fixed at 0 in its class 2, is latents'
#   class intercept and slope: threshold = -intercept of Mplus's class 2,
#   I = the intercepts' difference, S = each slope.
# - ex8.12, a hidden Markov model (one binary indicator at four occasions,
#   homogeneous transitions): lta(transitions = "homogeneous").
# - ex8.14, two occasions, five binary indicators, three states, a covariate
#   on the first state and on every origin's transitions:
#   lta(initial_covariates = "x", transition_covariates = "x").
# - ex8.15, a mover-stayer model at three occasions (the initial state
#   distribution differs between movers and stayers):
#   lta(mover_stayer = TRUE, transitions = "occasion").
#
# Compared: log likelihood, parameters, AIC, BIC; ex8.9's entropy, class
# proportions, modal counts and every estimate with its standard error; each
# state's item probabilities (with their errors in ex8.12); the transition
# tables of ex8.12 and ex8.14 (ex8.14's averaged over the persons'
# covariate, as both programs print it).
#
# Tolerances, fixed before the comparison: half a unit of the last printed
# digit plus 1e-6 of the value; counts exactly. Class sizes are compared as
# the proportions Mplus prints (five decimals), not its counts to 1e-5
# persons. Standard errors: Mplus's MLR, a sandwich without a small-sample
# factor; latents' robust errors carry the CR1 factor n / (n - 1)
# (.inference_cluster_factor()), so they are scaled by sqrt((n - 1) / n),
# n the persons.
#
# Not compared: Mplus's entropy of several latent class variables (over
# their joint pattern), ex8.15's transition tables and class counts (Mplus's
# mix movers and stayers; latents reports each), the chi-square tests.

# permutations() from the wider Latent GOLD kit.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

.ug8_reference <- "Mplus"

#' One comparison row; `printed` is the text Mplus printed
#' @noRd
.ug8_row <- function(case, fn, quantity, printed, obtained, exact = FALSE) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- if (exact) 0 else 0.5 * 10^(-digits) + 1e-6 * abs(value)
  difference <- abs(obtained - value)
  data.frame(case = case, reference = .ug8_reference,
             printed_in = sprintf("Mplus 8.8 User's Guide, %s.out", sub("_.*$", "", case)),
             quantity = quantity, reference_value = value, latents = obtained,
             difference = difference, tolerance = tolerance, n_values = 1L,
             kind = "commercial program", agrees = difference <= tolerance, fn = fn)
}

#' The lines of an Mplus output file
#' @noRd
.ug8_lines <- function(directory, example) {
  readLines(file.path(directory, "listings", paste0(example, ".out")), warn = FALSE)
}

#' The value printed after a label (the first such line)
#' @noRd
.ug8_value <- function(lines, label) {
  hit <- grep(label, lines, fixed = TRUE)[[1L]]
  fields <- strsplit(trimws(lines[[hit]]), "[[:space:]]+")[[1L]]
  fields[[length(fields)]]
}

#' Log likelihood, parameters, AIC, BIC (and entropy) rows
#' @noRd
.ug8_fit_rows <- function(case, fn, lines, fit, entropy = NULL) {
  log_likelihood <- as.numeric(stats::logLik(fit))
  p <- attr(stats::logLik(fit), "df") %||% fit$n_parameters
  # BIC on the persons: lta()'s bic_groups; a group-level mixture regression's BIC().
  bic <- if (identical(fn, "lta")) get_results(fit, "information_criteria")$bic_groups else
    stats::BIC(fit)
  rbind(.ug8_row(case, fn, "log likelihood", .ug8_value(lines, "H0 Value"), log_likelihood),
        .ug8_row(case, fn, "number of parameters", .ug8_value(lines, "Number of Free Parameters"),
                 p, exact = TRUE),
        .ug8_row(case, fn, "AIC", .ug8_value(lines, "Akaike (AIC)"), -2 * log_likelihood + 2 * p),
        .ug8_row(case, fn, "BIC", .ug8_value(lines, "Bayesian (BIC)"), bic),
        if (!is.null(entropy)) .ug8_row(case, fn, "entropy", .ug8_value(lines, "Entropy"), entropy))
}

#' Mplus's "MEAN/PROBABILITY PROFILES FOR <variable>" block
#' @return One row per indicator, category and class: printed probability.
#' @noRd
.ug8_profiles <- function(lines, variable) {
  start <- grep(sprintf("MEAN/PROBABILITY PROFILES FOR %s$", variable), lines)[[1L]]
  following <- lines[(start + 1L):length(lines)]
  end <- grep("MEAN/PROBABILITY PROFILES FOR|^[A-Z]", following)[[1L]] - 1L
  block <- trimws(following[seq_len(end)])
  block <- block[nzchar(block) & !grepl("^(Latent class|[0-9 ]+$)", block)]
  indicator <- cumsum(!grepl("^Category", block))
  names_seen <- block[!grepl("^Category", block)]
  rows <- grepl("^Category", block)
  fields <- strsplit(block[rows], "[[:space:]]+")
  do.call(rbind, lapply(seq_along(fields), function(i) {
    values <- fields[[i]][-(1:2)]
    data.frame(indicator = tolower(names_seen[indicator[rows][[i]]]),
               category = as.integer(fields[[i]][[2L]]), class = seq_along(values),
               printed = values, stringsAsFactors = FALSE)
  }))
}

#' Mplus's transition table "<from> Classes (Rows) by <to> Classes (Columns)"
#' @noRd
.ug8_transitions <- function(lines, from, to) {
  start <- grep(sprintf("%s Classes (Rows) by %s Classes (Columns)", from, to), lines,
                fixed = TRUE)[[1L]]
  following <- trimws(lines[(start + 1L):(start + 12L)])
  rows <- following[grepl("^[0-9]+ +[0-9]+[.]", following)]
  do.call(rbind, lapply(rows, function(row) {
    fields <- strsplit(row, "[[:space:]]+")[[1L]]
    data.frame(from = as.integer(fields[[1L]]), to = seq_along(fields[-1L]), printed = fields[-1L],
               stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# Latent transition models

#' Read a wide UG data file into lta()'s long form
#' @param times,items Occasion and indicator numbers: column u<time><item>
#'   (u<time> when there is one indicator).
#' @noRd
.ug8_long <- function(directory, example, times, items, extra = character()) {
  wide <- utils::read.table(file.path(directory, "listings", paste0(example, ".dat")))
  single <- length(items) == 1L
  columns <- as.vector(outer(items, times, function(i, t) if (single) paste0("u", t) else
    paste0("u", t, i)))
  names(wide)[seq_along(columns)] <- columns
  if (length(extra)) names(wide)[length(columns) + seq_along(extra)] <- extra
  wide$id <- seq_len(nrow(wide))
  indicators <- paste0("v", items)
  long <- do.call(rbind, lapply(times, function(t) {
    sources <- if (single) paste0("u", t) else paste0("u", t, items)
    occasion <- wide[c("id", extra, sources)]
    names(occasion)[length(extra) + 1L + seq_along(items)] <- indicators
    occasion$time <- t
    occasion
  }))
  long <- long[order(long$id, long$time), ]
  long[indicators] <- lapply(long[indicators], factor)
  rownames(long) <- NULL
  long
}

#' lta()'s responses with profiles numbered (a general transition fit labels
#' them "profile_1", ...)
#' @noRd
.ug8_responses <- function(fit) {
  responses <- as.data.frame(get_results(fit, "responses"))
  responses$profile <- as.integer(sub("^profile_", "", responses$profile))
  responses
}

#' latents' state order matching Mplus's printed profiles
#' @return `perm`: Mplus class c is latents profile perm[c].
#' @noRd
.ug8_state_order <- function(profiles, responses) {
  k <- max(profiles$class)
  category_one <- profiles[profiles$category == 1L, ]
  gaps <- vapply(permutations(k), function(perm) {
    obtained <- responses$probability[match(
      paste(perm[category_one$class], category_one$indicator, "0"),
      paste(responses$profile, responses$indicator, responses$category))]
    max(abs(obtained - as.numeric(category_one$printed)))
  }, numeric(1))
  permutations(k)[[which.min(gaps)]]
}

#' Item probability rows (Mplus's category 1 is the response 0)
#' @noRd
.ug8_profile_rows <- function(case, profiles, responses, perm, rename) {
  do.call(rbind, lapply(seq_len(nrow(profiles)), function(i) {
    row <- profiles[i, ]
    obtained <- responses$probability[responses$profile == perm[[row$class]] &
                                        responses$indicator == rename(row$indicator) &
                                        responses$category == as.character(row$category - 1L)]
    .ug8_row(case, "lta", sprintf("probability class %d %s category %d", row$class,
                                  row$indicator, row$category), row$printed, obtained)
  }))
}

#' Transition rows: Mplus's from/to classes against latents' profiles
#' @noRd
.ug8_transition_rows <- function(case, printed, transitions, perm, label) {
  do.call(rbind, lapply(seq_len(nrow(printed)), function(i) {
    row <- printed[i, ]
    obtained <- transitions$probability[
      transitions$from %in% c(perm[[row$from]], sprintf("profile_%d", perm[[row$from]])) &
        transitions$to %in% c(perm[[row$to]], sprintf("profile_%d", perm[[row$to]]))]
    stopifnot("one latents transition per printed cell" = length(obtained) == 1L)
    .ug8_row(case, "lta", sprintf("transition %s %d to %d", label, row$from, row$to),
             row$printed, obtained)
  }))
}

compare_ex8_12 <- function(directory) {
  lines <- .ug8_lines(directory, "ex8.12")
  data <- .ug8_long(directory, "ex8.12", 1:4, 1L)
  fit <- lta(data, "v1", id = "id", n_profiles = 2L, time = "time", categorical = "v1",
             transitions = "homogeneous", n_starts = 20L, seed = 1L, tol = 1e-12)
  n <- length(unique(data$id))
  responses <- .ug8_responses(fit)
  profiles <- .ug8_profiles(lines, "C1")
  profiles$indicator <- "v1"
  perm <- .ug8_state_order(profiles, responses)
  # RESULTS IN PROBABILITY SCALE, C1's part (the measurement is invariant).
  start <- grep("^Results for Class-specific Model Parts of C1$", lines)[[1L]]
  block <- trimws(lines[start + seq_len(10L)])
  class <- cumsum(grepl("^Latent Class C1#", block))
  category <- grepl("^Category", block)
  fields <- strsplit(block[category], "[[:space:]]+")
  inference <- parameter_inference(fit, vcov_type = "robust", boundary = "fix")
  measurement <- inference[inference$level == "measurement", ]
  error_rows <- do.call(rbind, lapply(seq_along(fields), function(i) {
    c <- class[category][[i]]
    j <- as.integer(fields[[i]][[2L]])
    se <- measurement$standard_error[measurement$outcome == sprintf("profile_%d", perm[[c]]) &
                                       measurement$term == sprintf("v1:%d", j - 1L)]
    .ug8_row("ex8.12_hmm", "lta", sprintf("s.e. probability class %d u category %d", c, j),
             fields[[i]][[4L]], se * sqrt((n - 1) / n))
  }))
  rbind(.ug8_fit_rows("ex8.12_hmm", "lta", lines, fit),
        .ug8_profile_rows("ex8.12_hmm", profiles, responses, perm, identity),
        error_rows,
        .ug8_transition_rows("ex8.12_hmm", .ug8_transitions(lines, "C1", "C2"),
                             get_results(fit, "transitions"), perm, "C1 C2"))
}

compare_ex8_14 <- function(directory) {
  lines <- .ug8_lines(directory, "ex8.14")
  data <- .ug8_long(directory, "ex8.14", 1:2, 1:5, "x")
  fit <- lta(data, paste0("v", 1:5), id = "id", n_profiles = 3L, time = "time",
             categorical = paste0("v", 1:5), initial_covariates = "x",
             transition_covariates = "x", n_starts = 20L, seed = 1L, tol = 1e-12)
  responses <- .ug8_responses(fit)
  profiles <- .ug8_profiles(lines, "C1")
  profiles$indicator <- sub("^u1", "v", profiles$indicator)
  perm <- .ug8_state_order(profiles, responses)
  rbind(.ug8_fit_rows("ex8.14_covariate", "lta", lines, fit),
        .ug8_profile_rows("ex8.14_covariate", profiles, responses, perm, identity),
        .ug8_transition_rows("ex8.14_covariate", .ug8_transitions(lines, "C1", "C2"),
                             get_results(fit, "transitions"), perm, "C1 C2"))
}

compare_ex8_15 <- function(directory) {
  lines <- .ug8_lines(directory, "ex8.15")
  data <- .ug8_long(directory, "ex8.15", 1:3, 1:5)
  # Movers' first transitions include one at zero: lta() warns that the fit
  # is on a boundary (its standard errors are withheld; none is compared).
  fit <- lta(data, paste0("v", 1:5), id = "id", n_profiles = 3L, time = "time",
             categorical = paste0("v", 1:5), mover_stayer = TRUE, transitions = "occasion",
             n_starts = 20L, seed = 1L, tol = 1e-12)
  responses <- .ug8_responses(fit)
  profiles <- .ug8_profiles(lines, "C1")
  profiles$indicator <- sub("^u1", "v", profiles$indicator)
  perm <- .ug8_state_order(profiles, responses)
  rbind(.ug8_fit_rows("ex8.15_mover_stayer", "lta", lines, fit),
        .ug8_profile_rows("ex8.15_mover_stayer", profiles, responses, perm, identity))
}

# ---------------------------------------------------------------------------
# ex8.9: a latent class growth analysis of a binary outcome

compare_ex8_9 <- function(directory) {
  lines <- .ug8_lines(directory, "ex8.9")
  wide <- utils::read.table(file.path(directory, "listings", "ex8.9.dat"))[1:4]
  names(wide) <- paste0("u", 1:4)
  wide$id <- seq_len(nrow(wide))
  long <- stats::reshape(wide, direction = "long", varying = paste0("u", 1:4), v.names = "u",
                         timevar = "time", times = 0:3, idvar = "id")
  long <- long[order(long$id, long$time), ]
  n <- nrow(wide)
  fit <- mixture_regression(u ~ time, long, n_classes = 2L, family = "binomial", id = "id",
                            class_level = "group", n_starts = 20L, seed = 1L, tol = 1e-12,
                            vcov_type = "robust")
  theta <- stats::coef(fit)
  v <- stats::vcov(fit) * (n - 1) / n
  # MODEL RESULTS: estimate and standard error of each free parameter.
  start <- grep("^MODEL RESULTS$", lines)[[1L]]
  results <- trimws(lines[start:(start + 60L)])
  printed_at <- function(pattern, which = 1L) {
    fields <- strsplit(results[grep(pattern, results)][[which]], "[[:space:]]+")[[1L]]
    fields[2:3]
  }
  printed <- list(`class 1 mean S` = printed_at("^S +-?[0-9]", 1L),
                  `class 2 mean S` = printed_at("^S +-?[0-9]", 2L),
                  `class 1 mean I` = printed_at("^I +-?[0-9]", 1L),
                  `threshold U$1` = printed_at("^U1[$]1 ", 1L),
                  `C#1` = printed_at("^C#1 "))
  threshold <- as.numeric(printed[["threshold U$1"]][[1L]])
  # Mplus's class 2 has I fixed at 0; it is the latents class whose
  # intercept is the negative of the printed threshold.
  intercepts <- theta[grep("^coefficient[.]class_[0-9]+[.][(]Intercept[)]$", names(theta))]
  latents_class <- function(name) as.integer(sub("^coefficient[.]class_([0-9]+)[.].*$", "\\1", name))
  second <- latents_class(names(intercepts)[[which.min(abs(intercepts + threshold))]])
  first <- setdiff(1:2, second)
  name <- function(class, term) sprintf("coefficient.class_%d.%s", class, term)
  share_logit <- grep("^membership[.]class_2[.][(]Intercept[)]$", names(theta), value = TRUE)
  # Each printed parameter as a linear function g'theta of latents' parameters.
  linear <- function(weights) {
    g <- stats::setNames(numeric(length(theta)), names(theta))
    g[names(weights)] <- weights
    c(estimate = sum(g * theta), se = sqrt(drop(t(g) %*% v %*% g)))
  }
  sign <- if (first == 1L) -1 else 1  # C#1 = log(p_first / p_second)
  terms <- list(
    `class 1 mean I` = linear(stats::setNames(c(1, -1), c(name(first, "(Intercept)"),
                                                         name(second, "(Intercept)")))),
    `class 1 mean S` = linear(stats::setNames(1, name(first, "time"))),
    `class 2 mean S` = linear(stats::setNames(1, name(second, "time"))),
    `threshold U$1` = linear(stats::setNames(-1, name(second, "(Intercept)"))),
    `C#1` = linear(stats::setNames(sign, share_logit)))
  estimate_rows <- do.call(rbind, lapply(names(printed), function(term) rbind(
    .ug8_row("ex8.9_lcga", "growth", paste("estimate", term), printed[[term]][[1L]],
             terms[[term]][["estimate"]]),
    .ug8_row("ex8.9_lcga", "growth", paste("s.e.", term), printed[[term]][[2L]],
             terms[[term]][["se"]]))))
  classes <- as.data.frame(get_results(fit, "classes"))
  share <- classes$share[order(c(first, second))]
  modal <- classes$n_assigned[order(c(first, second))]
  proportions <- grep("^BASED ON THE ESTIMATED MODEL$", lines)[[1L]]
  shares_printed <- vapply(1:2, function(c) {
    strsplit(trimws(lines[grep(sprintf("^ +%d +[0-9.]+ +[0-9.]+$", c), lines)[
      grep(sprintf("^ +%d +[0-9.]+ +[0-9.]+$", c), lines) > proportions][[1L]]]),
      "[[:space:]]+")[[1L]][[3L]]
  }, character(1))
  modal_start <- grep("^BASED ON THEIR MOST LIKELY LATENT CLASS MEMBERSHIP$", lines)[[1L]]
  modal_printed <- vapply(1:2, function(c) {
    candidates <- grep(sprintf("^ +%d +[0-9]+ +[0-9.]+$", c), lines)
    strsplit(trimws(lines[candidates[candidates > modal_start][[1L]]]), "[[:space:]]+")[[1L]][[2L]]
  }, character(1))
  rbind(.ug8_fit_rows("ex8.9_lcga", "growth", lines, fit,
                      entropy = as.data.frame(get_results(fit, "fit"))$entropy),
        do.call(rbind, lapply(1:2, function(c) rbind(
          .ug8_row("ex8.9_lcga", "growth", sprintf("class %d proportion", c), shares_printed[[c]],
                   share[[c]]),
          .ug8_row("ex8.9_lcga", "growth", sprintf("class %d modal count", c), modal_printed[[c]],
                   modal[[c]], exact = TRUE)))),
        estimate_rows)
}

compare_all_ug8 <- function(directory) {
  rbind(compare_ex8_9(directory), compare_ex8_12(directory), compare_ex8_14(directory),
        compare_ex8_15(directory))
}
