# Shared by make-kit.R, compare.R and equivalence/test-latentgold-step-one.R:
# the simulated cases, an independent likelihood for each step-one model, the
# fixed step-three parameters as functions of the step-one parameters, the
# Latent GOLD syntax and the comparison. Needs the package loaded (for the
# fits), numDeriv, and equivalence/latentgold/compare-functions.R sourced (for
# the listing parsers).
#
# The check (Bakk, Oberski and Vermunt 2014; Latent GOLD 6.1 Syntax Guide
# 4.3.12): Latent GOLD's step-three model reads the classification-error
# logits e and the class-size logits s as fixed parameters, together with
# their Jacobian towards the step-one parameters and the step-one covariance,
# and adds that uncertainty to its standard errors. The Jacobian and the
# covariance are computed here from a likelihood written independently of the
# package, so the comparison tests latents' whole step-one path against
# Latent GOLD's step-three propagation. Latent GOLD's own step-one run writes
# the same quantities for its own parameterization (writeclassificationerrors),
# which checks the covariance of e independently of this file.

step_cases <- c("s01", "s02", "s03")
step_outcomes <- c("z", "cat")

#' Simulate one case and fit its step-one model with latents
#'   s01: lpa(), 2 profiles, 2 Gaussian indicators (VVI), 600 rows.
#'   s02: lpa(), 3 profiles, 3 Gaussian indicators (VVI), 900 rows.
#'   s03: lca(), 3 classes, 6 binary items, 1000 rows.
#' Each has a numeric outcome z and a categorical outcome cat.
step_case <- function(name) {
  set.seed(switch(name, s01 = 20261111, s02 = 20261112, s03 = 20261113))
  n <- switch(name, s01 = 600L, s02 = 900L, s03 = 1000L)
  k <- if (identical(name, "s01")) 2L else 3L
  truth <- sample.int(k, n, replace = TRUE, prob = switch(name, s01 = c(0.55, 0.45),
                                                            c(0.4, 0.35, 0.25)))
  if (identical(name, "s03")) {
    item <- rbind(c(0.85, 0.8, 0.8, 0.75, 0.8, 0.7), c(0.2, 0.3, 0.8, 0.75, 0.25, 0.3),
                  c(0.15, 0.2, 0.25, 0.2, 0.2, 0.25))
    data <- as.data.frame(lapply(stats::setNames(seq_len(6L), paste0("u", seq_len(6L))),
                                 function(j) factor(stats::rbinom(n, 1L, item[truth, j]),
                                                    levels = 0:1)))
    indicators <- paste0("u", seq_len(6L))
  } else {
    centre <- switch(name, s01 = rbind(c(0, 0), c(1.8, 1.4)),
                     rbind(c(0, 0, 0), c(1.6, 1.2, 0.4), c(0.2, 1.8, 1.6)))
    spread <- switch(name, s01 = rbind(c(1, 1), c(1.2, 0.8)),
                     rbind(c(1, 1, 1), c(0.8, 1.2, 1), c(1, 0.7, 1.3)))
    data <- as.data.frame(lapply(stats::setNames(seq_len(ncol(centre)),
                                                 paste0("y", seq_len(ncol(centre)))),
                                 function(j) stats::rnorm(n, centre[truth, j],
                                                          sqrt(spread[truth, j]))))
    indicators <- names(data)
  }
  data$z <- stats::rnorm(n, c(0, 1, 0.5)[truth], c(1, 1.3, 0.8)[truth])
  categories <- rbind(c(0.6, 0.3, 0.1), c(0.2, 0.3, 0.5), c(0.3, 0.5, 0.2))
  data$cat <- factor(vapply(truth, function(class) {
    sample(c("low", "middle", "high"), 1L, prob = categories[class, ])
  }, character(1)), levels = c("low", "middle", "high"))
  fit <- if (identical(name, "s03")) {
    lca(data, indicators, n_classes = k, n_starts = 10, seed = 1, tol = 1e-12,
        max_iter = 20000)
  } else {
    lpa(data, indicators, n_profiles = k, model = "VVI", n_starts = 10, seed = 1,
        tol = 1e-12, max_iter = 20000)
  }
  list(name = name, data = data, fit = fit, indicators = indicators, k = k,
       binary = identical(name, "s03"))
}

#' The independent step-one model: per-unit log likelihood and posteriors
#'
#' Coordinates: class means (class by class) and log variances, or the logit
#' of P(u = 1) per class and item; then the size logits log(pi_k / pi_1),
#' k = 2..K.
step_one_model <- function(made) {
  k <- made$k
  y <- if (made$binary) {
    vapply(made$indicators, function(name) as.numeric(made$data[[name]] == "1"),
           numeric(nrow(made$data)))
  } else as.matrix(made$data[made$indicators])
  d <- ncol(y)
  blocks <- if (made$binary) k * d else 2L * k * d
  joint <- function(theta) {
    sizes <- c(0, theta[blocks + seq_len(k - 1L)])
    log_sizes <- sizes - log(sum(exp(sizes)))
    vapply(seq_len(k), function(class) {
      rows <- (class - 1L) * d + seq_len(d)
      density <- if (made$binary) {
        p <- stats::plogis(theta[rows])
        rowSums(y %*% diag(log(p), d) + (1 - y) %*% diag(log1p(-p), d))
      } else {
        means <- theta[rows]
        sds <- exp(0.5 * theta[k * d + rows])
        rowSums(vapply(seq_len(d), function(j) {
          stats::dnorm(y[, j], means[j], sds[j], log = TRUE)
        }, numeric(nrow(y))))
      }
      log_sizes[class] + density
    }, numeric(nrow(y)))
  }
  unit_log_likelihood <- function(theta) {
    values <- joint(theta)
    largest <- apply(values, 1L, max)
    largest + log(rowSums(exp(values - largest)))
  }
  posteriors <- function(theta) exp(joint(theta) - unit_log_likelihood(theta))
  list(unit_log_likelihood = unit_log_likelihood, posteriors = posteriors,
       log_likelihood = function(theta) sum(unit_log_likelihood(theta)))
}

#' latents' estimate in the independent coordinates, polished to their maximum
step_one_estimate <- function(made, model) {
  fit <- made$fit
  sizes <- as.vector(fit$profile_probabilities)
  start <- if (made$binary) {
    c(as.vector(t(vapply(made$indicators, function(name) {
      stats::qlogis(fit$response_probabilities[[name]][, "1"])
    }, numeric(made$k)))), log(sizes[-1L] / sizes[1L]))
  } else {
    c(as.vector(t(fit$means)), as.vector(t(log(fit$variances))), log(sizes[-1L] / sizes[1L]))
  }
  polished <- stats::optim(start, function(theta) -model$log_likelihood(theta), method = "BFGS",
                           control = list(reltol = 1e-15, maxit = 5000L))
  list(theta = polished$par, log_likelihood = -polished$value,
       covariance = solve(-numDeriv::hessian(model$log_likelihood, polished$par)))
}

#' The fixed step-three parameters from the units' posteriors
#'
#' e: Latent GOLD's `~err` logits, log D[t, s] / D[t, t] for s != t, true
#' class by true class; s: log(pi_k / pi_1) with pi the average posteriors.
#' The modal assignments are the data of step three and stay fixed.
step_values <- function(posteriors, modal, k) {
  errors <- vapply(seq_len(k), function(assigned) {
    colSums(posteriors * (modal == assigned))
  }, numeric(k)) / colSums(posteriors)
  logits <- unlist(lapply(seq_len(k), function(true_class) {
    log(errors[true_class, -true_class] / errors[true_class, true_class])
  }))
  sizes <- colMeans(posteriors)
  c(logits, log(sizes[-1L] / sizes[1L]))
}

#' The same, as functions of the independent step-one coordinates
step_fixed_values <- function(theta, model, modal, k) {
  step_values(model$posteriors(theta), modal, k)
}

#' The step-three data: id, latents' modal assignment W and the outcomes
step_data_file <- function(made) {
  data.frame(id = seq_len(nrow(made$data)), W = made$fit$subject_profiles,
             z = made$data$z, cat = as.integer(made$data$cat))
}

#' The file Latent GOLD reads the fixed values, Jacobian and covariance from
#'
#' Latent GOLD 6.1 adds the step-one uncertainty only when the file has a
#' Step1Variances block (here the same covariance as Step1Hessian): with
#' parameters, Jacobian and Hessian alone its standard errors are those of
#' plain fixed values, whatever the manual says. With Step1Gradients (one row
#' per case, in case-id order) it applies the second-order correction.
step_fixed_file <- function(values, jacobian, covariance, gradients = NULL) {
  rows <- function(matrix_form) {
    apply(matrix_form, 1L, function(row) paste(sprintf("%.17g", row), collapse = " "))
  }
  c("Step1Parameters", paste(sprintf("%.17g", values), collapse = " "),
    "Step1Jacobian", as.character(ncol(jacobian)), rows(jacobian),
    "Step1Hessian", rows(covariance), "Step1Variances", rows(covariance),
    if (!is.null(gradients)) c("Step1Gradients", rows(gradients)))
}

step_header <- function(name, infile, robust, extra = "") {
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s'", infile), "", "model",
    sprintf("title '%s';", name), "options", "   maxthreads=all;", "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261101 sets=20 tolerance=1e-005 iterations=100;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;", "   output",
    sprintf("      parameters=first standarderrors%s profile estimatedvalues=model%s;",
            if (robust) "=robust" else "", extra))
}

#' Step three: e and s fixed, read with their covariance (`estimated`, the
#' first-order correction), with it and the case gradients (`second`), or as
#' plain values (`fixed`); the outcome's class parameters free
step_three_syntax <- function(case, outcome, variant, robust, k, values) {
  name <- sprintf("%s_%s_%s%s", case, outcome, variant, if (robust) "_robust" else "")
  dependent <- if (identical(outcome, "z")) "z continuous" else "cat nominal coding=first"
  source <- switch(variant,
    estimated = sprintf("'%s_fixed.txt'", case),
    second = sprintf("'%s_fixed_gradients.txt'", case),
    sprintf("{%s}", paste(sprintf("%.17g", values), collapse = " ")))
  c(step_header(name, sprintf("%s.dat", case), robust),
    "variables", "   caseid id;", sprintf("   dependent W nominal, %s;", dependent),
    sprintf("   latent Cluster nominal %d coding=first;", k), "equations",
    "   Cluster <- (s) 1;", "   W <- (e~err) 1 | Cluster;",
    sprintf("   %s <- 1 + Cluster;", outcome),
    if (identical(outcome, "z")) "   z | Cluster;",
    sprintf("   e s = %s;", source), "end model")
}

#' Latent GOLD's own step one, writing its classification-error file
step_one_syntax <- function(made) {
  case <- made$name
  dependent <- paste(sprintf("%s %s", made$indicators,
                             if (made$binary) "nominal" else "continuous"), collapse = ", ")
  c(step_header(sprintf("%s_step1", case), sprintf("%s_step1.dat", case), FALSE,
                sprintf(" writeclassificationerrors='%s_lg_errors.txt'", case)),
    sprintf("   outfile '%s_lg_classes.txt' classification keep id;", case),
    "variables", "   caseid id;", sprintf("   dependent %s;", dependent),
    sprintf("   latent Cluster nominal %d;", made$k), "equations", "   Cluster <- 1;",
    sprintf("   %s <- 1 + Cluster;", made$indicators),
    if (!made$binary) sprintf("   %s | Cluster;", made$indicators), "end model")
}

# The comparison (see compare.R for what each row compares).

row <- function(case, outcome, step_one, quantity, ours, theirs, tolerance, note = "") {
  difference <- abs(ours - theirs)
  status <- if (is.na(theirs) && nzchar(note)) "not comparable" else
    if (is.na(theirs)) "missing" else
    if (isTRUE(difference <= tolerance)) "agree" else "differ"
  data.frame(case = case, outcome = outcome, step_one = step_one, quantity = quantity,
             latents = ours, latent_gold = theirs, difference = difference,
             tolerance = tolerance, status = status, note = note)
}

# Latent GOLD's file: the logits, the Jacobian and the step-one covariance.
read_errors_file <- function(path) {
  lines <- readLines(path)
  numbers <- function(text) as.numeric(strsplit(trimws(text), "[[:space:]]+")[[1L]])
  at <- function(keyword) match(keyword, trimws(lines))
  values <- numbers(lines[at("Step1Parameters") + 1L])
  width <- as.integer(lines[at("Step1Jacobian") + 1L])
  jacobian <- do.call(rbind, lapply(lines[at("Step1Jacobian") + 1L + seq_along(values)], numbers))
  covariance <- do.call(rbind, lapply(lines[at("Step1Hessian") + seq_len(width)], numbers))
  list(values = values, covariance = jacobian %*% covariance %*% t(jacobian))
}

# The ~err logits are ordered by true class, then by assigned class without
# the true one; `order` maps Latent GOLD's classes onto latents'.
error_positions <- function(k) {
  do.call(rbind, lapply(seq_len(k), function(true_class) {
    data.frame(true_class = true_class, assigned = setdiff(seq_len(k), true_class))
  }))
}

error_rows <- function(case, made, returned, independent) {
  path <- file.path(returned, sprintf("%s_lg_errors.txt", case))
  classes <- .lg_read_table(file.path(returned, sprintf("%s_lg_classes.txt", case)))
  lg_modal <- classes[["Cluster#"]]
  if (!file.exists(path) || is.null(lg_modal)) {
    return(row(case, "", "", "step-one error covariance", NA_real_, NA_real_, 0,
               "Latent GOLD's step-one files are missing."))
  }
  agreement <- table(factor(lg_modal, seq_len(made$k)),
                     factor(made$fit$subject_profiles, seq_len(made$k)))
  order <- apply(agreement, 1L, which.max)
  theirs <- read_errors_file(path)
  ours <- independent[[case]]
  positions <- error_positions(made$k)
  # Our logit for each of Latent GOLD's: the same pair of classes, relabelled.
  key <- paste(positions$true_class, positions$assigned)
  mapped <- match(paste(order[positions$true_class], order[positions$assigned]), key)
  count <- nrow(positions)
  our_covariance <- (ours$jacobian %*% ours$covariance %*% t(ours$jacobian))[
    mapped, mapped, drop = FALSE]
  # Latent GOLD's own search can end at another maximum (a heteroscedastic
  # Gaussian mixture has degenerate ones, a class of a few cases with a
  # near-zero variance); its file then describes another model.
  lg_log_likelihood <- .lg_listing_value(
    .lg_listing_lines(file.path(returned, sprintf("%s_step1.lst", case))), "Log-likelihood (LL)")$value
  regular <- isTRUE(all(sort(order) == seq_len(made$k))) &&
    isTRUE(abs(lg_log_likelihood - made$fit$log_likelihood) < 1e-3)
  note <- if (regular) "" else sprintf(paste(
    "Latent GOLD's own step one reached another maximum (LL %.4f against %.4f),",
    "so its file is not comparable."), lg_log_likelihood, made$fit$log_likelihood)
  scale <- max(abs(our_covariance))
  rbind(
    row(case, "", "", "step-one error logits (max abs difference)",
        0, if (regular) max(abs(theirs$values - ours$values[mapped])) else NA_real_,
        1e-3 * max(abs(ours$values[seq_len(count)])), note),
    row(case, "", "", "step-one error covariance (max abs difference)",
        0, if (regular) max(abs(theirs$covariance - our_covariance)) else NA_real_,
        1e-3 * scale, note))
}

wald_value <- function(lines, outcome) {
  fields <- strsplit(lines, "\t", fixed = TRUE)
  hit <- which(vapply(fields, function(cells) {
    length(cells) >= 3L && trimws(cells[[1L]]) %in% c(outcome, paste0(outcome, "(1)")) &&
      identical(trimws(cells[[2L]]), "<-") && identical(trimws(cells[[3L]]), "Cluster(1)")
  }, logical(1)))
  if (length(hit) == 0L) return(c(statistic = NA_real_, df = NA_real_))
  cells <- trimws(fields[[hit[[1L]]]][-(1:3)])
  numbers <- as.numeric(cells[vapply(cells, .lg_is_number, logical(1))])
  c(statistic = numbers[[2L]], df = numbers[[3L]])
}

step_rows <- function(case, made, outcome, step_one, returned) {
  path <- file.path(returned, sprintf("%s_%s_%s_robust.lst", case, outcome, step_one))
  if (!file.exists(path)) {
    return(row(case, outcome, step_one, "listing", NA_real_, NA_real_, 0, "No listing."))
  }
  lines <- .lg_listing_lines(path)
  values <- subset(.lg_estimated_values(lines), response == outcome & condition == "Cluster")
  ours <- three_step(made$fit, made$data, outcome, method = "ml", vcov_type = "independent",
                     step_one = step_one)
  overall <- three_step(made$fit, made$data, outcome, method = "ml",
                        vcov_type = "independent", contrast = "overall", step_one = step_one)
  label <- if (identical(outcome, "z")) rep("Mean", nrow(ours)) else
    as.character(match(ours$category, levels(made$data$cat)))
  cells <- Map(function(class, code, estimate, error) {
    hit <- subset(values, condition_value == as.character(class) & label == code)
    name <- sprintf("class %d%s", class, if (identical(code, "Mean")) "" else
      sprintf(", category %s", code))
    if (nrow(hit) != 1L) {
      return(row(case, outcome, step_one, paste("estimate:", name), estimate, NA_real_, 0,
                 "Not found under EstimatedValues-Model."))
    }
    rbind(row(case, outcome, step_one, paste("estimate:", name), estimate, hit$estimate,
              .lg_half_unit(hit$estimate_text) + 1e-4 * abs(hit$estimate)),
          row(case, outcome, step_one, paste("robust standard error:", name), error,
              hit$standard_error,
              .lg_half_unit(hit$standard_error_text) + 1e-4 * abs(hit$standard_error)))
  }, ours$class, label, ours$estimate, ours$standard_error)
  wald <- wald_value(lines, outcome)
  rbind(do.call(rbind, cells),
        row(case, outcome, step_one, "Wald: outcome <- Cluster", overall$statistic,
            wald[["statistic"]], 5e-5 + 1e-4 * abs(wald[["statistic"]])),
        row(case, outcome, step_one, "Wald df", overall$df, wald[["df"]], 0))
}

#' Every comparison row: Latent GOLD's step-one file, then each outcome's
#' step-three listings with the step-one estimates fixed and estimated
step_comparison <- function(returned, independent) {
  do.call(rbind, lapply(step_cases, function(case) {
    made <- step_case(case)
    rbind(error_rows(case, made, returned, independent),
          do.call(rbind, lapply(step_outcomes, function(outcome) {
            rbind(step_rows(case, made, outcome, "fixed", returned),
                  step_rows(case, made, outcome, "estimated", returned))
          })))
  }))
}
