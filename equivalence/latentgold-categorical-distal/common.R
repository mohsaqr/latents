# Shared by make-kit.R, compare.R and equivalence/test-latentgold-categorical-distal.R:
# the simulated cases, the Latent GOLD syntax, latents' tables and the
# comparison. Needs the package loaded and equivalence/latentgold/compare-functions.R
# sourced (for the listing parsers).

distal_cases <- c("d01", "d02")
distal_methods <- c(bch = "modal bch", modal = "modal none", proportional = "proportional none",
                    ml = "modal ml")

#' Simulate one case and fit its measurement model with latents
#'   d01: lpa(), two profiles, two Gaussian indicators, 600 rows, a distal
#'        outcome with three categories.
#'   d02: lca(), three classes, five binary items, 900 rows, a binary distal.
distal_case <- function(name) {
  set.seed(switch(name, d01 = 20261101, d02 = 20261102))
  if (identical(name, "d01")) {
    truth <- sample.int(2L, 600L, replace = TRUE, prob = c(0.55, 0.45))
    data <- data.frame(y1 = stats::rnorm(600L, c(0, 1.8)[truth]),
                       y2 = stats::rnorm(600L, c(0, 1.4)[truth]))
    probabilities <- rbind(c(0.6, 0.3, 0.1), c(0.2, 0.3, 0.5))
    data$distal <- vapply(truth, function(class) {
      sample(c("low", "middle", "high"), 1L, prob = probabilities[class, ])
    }, character(1))
    data$distal <- factor(data$distal, levels = c("low", "middle", "high"))
    fit <- lpa(data, c("y1", "y2"), n_profiles = 2, model = "VVI", n_starts = 5, seed = 1)
    return(list(data = data, fit = fit))
  }
  truth <- sample.int(3L, 900L, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  item_probability <- rbind(c(0.85, 0.8, 0.75, 0.8, 0.7), c(0.2, 0.25, 0.8, 0.75, 0.3),
                            c(0.15, 0.2, 0.2, 0.25, 0.2))
  items <- as.data.frame(lapply(stats::setNames(seq_len(5L), paste0("u", seq_len(5L))),
                                function(item) {
    factor(stats::rbinom(900L, 1L, item_probability[truth, item]), levels = 0:1)
  }))
  items$distal <- factor(stats::rbinom(900L, 1L, c(0.7, 0.4, 0.15)[truth]), levels = 0:1,
                         labels = c("no", "yes"))
  fit <- lca(items, paste0("u", seq_len(5L)), n_classes = 3, n_starts = 5, seed = 1)
  list(data = items, fit = fit)
}

#' The data file Latent GOLD reads: the step-one posteriors and the outcome
#' as its category code (1, 2, ... in the factor's level order)
distal_data_file <- function(made) {
  posteriors <- made$fit$subject_posteriors
  colnames(posteriors) <- sprintf("post%d", seq_len(ncol(posteriors)))
  data.frame(id = seq_len(nrow(posteriors)), posteriors,
             distal = as.integer(made$data$distal))
}

distal_syntax <- function(case, method, robust, n_classes) {
  posterior <- paste(sprintf("post%d", seq_len(n_classes)), collapse = " ")
  name <- sprintf("%s_%s%s", case, method, if (robust) "_robust" else "")
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s.dat'", case), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=all;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=5000 nriterations=500;",
    "   startvalues", "      seed=20261101 sets=20 tolerance=1e-005 iterations=100;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;", sprintf("   step3 %s;", distal_methods[[method]]),
    "   output",
    sprintf("      parameters=first standarderrors%s profile estimatedvalues=model;",
            if (robust) "=robust" else ""),
    "variables", "   dependent distal nominal coding=first;",
    sprintf("   latent Cluster nominal coding=first posterior = ( %s );", posterior),
    "equations",
    # The ML model fixes the class sizes at the average posterior
    # probabilities, as three_step(method = "ml") does for every outcome and
    # Mplus's DU3STEP does; Latent GOLD does so when the sizes have no equation.
    if (!identical(method, "ml")) "   Cluster <- 1;",
    "   distal <- 1 + Cluster;", "end model")
}

#' latents' tables for one case, every method, with the variance Latent
#' GOLD's robust standard errors use (independent units)
distal_reference <- function(made) {
  methods <- stats::setNames(names(distal_methods), names(distal_methods))
  run <- function(contrast) lapply(methods, function(method) {
    three_step(made$fit, made$data, "distal", method = method, vcov_type = "independent",
               contrast = contrast)
  })
  list(tables = run("none"), overall = run("overall"), n = nrow(made$data),
       categories = levels(made$data$distal))
}

# Latent GOLD's Wald test of `distal <- Cluster`, on the row of the first
# category's first class effect.
.lg_distal_wald <- function(lines) {
  fields <- strsplit(lines, "\t", fixed = TRUE)
  hit <- which(vapply(fields, function(row) {
    length(row) >= 9L && identical(trimws(row[1:3]), c("distal(1)", "<-", "Cluster(1)"))
  }, logical(1)))
  if (length(hit) != 1L) return(c(statistic = NA_real_, df = NA_real_))
  row <- fields[[hit]]
  c(statistic = .lg_as_number(row[[8L]]), df = .lg_as_number(row[[9L]]))
}

distal_compare_listing <- function(case, method, robust, reference, returned) {
  name <- sprintf("%s_%s%s", case, method, if (robust) "_robust" else "")
  path <- file.path(returned, paste0(name, ".lst"))
  standard <- if (robust) "robust" else "model-based"
  row <- function(quantity, ours, theirs, tolerance, note = "") {
    difference <- abs(ours - theirs)
    status <- if (is.na(tolerance)) "not comparable" else
      if (isTRUE(difference <= tolerance)) "agree" else "differ"
    data.frame(case = case, method = method, standard_errors = standard, quantity = quantity,
               latents = ours, latent_gold = theirs, difference = difference,
               tolerance = tolerance, status = status, note = note)
  }
  if (!file.exists(path)) {
    return(row("listing", NA_real_, NA_real_, 0, sprintf("No %s.lst.", name)))
  }
  lines <- .lg_listing_lines(path)
  values <- subset(.lg_estimated_values(lines), response == "distal" & condition == "Cluster")
  table <- reference$tables[[method]]
  # Latent GOLD reads the outcome as its category code.
  table$code <- as.character(match(table$category, reference$categories))
  cells <- Map(function(class, code, category, estimate, error) {
    hit <- subset(values, condition_value == as.character(class) & label == code)
    label <- sprintf("class %d, %s", class, category)
    if (nrow(hit) != 1L) return(row(paste("probability:", label), estimate, NA_real_, 0,
                                    "Not found under EstimatedValues-Model."))
    rbind(
      row(paste("probability:", label), estimate, hit$estimate,
          .lg_half_unit(hit$estimate_text) + 1e-8),
      if (robust && identical(method, "proportional")) {
        row(paste("robust standard error:", label), error, hit$standard_error, NA_real_,
            paste("Latent GOLD treats the posterior weights as frequency weights:",
                  "its value is sqrt(p (1 - p) / T_k * n / (n - 1)), reproduced to every",
                  "printed digit; latents' is the sandwich of the weighted mean."))
      } else if (robust) {
        # Both sandwiches carry the small-sample factor n / (n - 1).
        row(paste("robust standard error:", label), error, hit$standard_error,
            .lg_half_unit(hit$standard_error_text) + 1e-8,
            "latents' vcov_type = \"independent\".")
      } else {
        row(paste("standard error:", label), error, hit$standard_error, NA_real_,
            "This package's is a sandwich; Latent GOLD's is model-based.")
      })
  }, table$class, table$code, table$category, table$estimate, table$standard_error)
  wald <- .lg_distal_wald(lines)
  overall <- reference$overall[[method]]
  test <- if (robust && !identical(method, "proportional")) {
    rbind(row("Wald: distal <- Cluster", overall$statistic, wald[["statistic"]],
              5e-5 + 1e-6 * abs(wald[["statistic"]]),
              "latents' logit-scale Wald with vcov_type = \"independent\"."),
          row("Wald df", overall$df, wald[["df"]], 0))
  } else {
    row("Wald: distal <- Cluster", overall$statistic, wald[["statistic"]], NA_real_,
        if (robust) "Latent GOLD's treats the posterior weights as frequency weights." else
          "Latent GOLD's uses the model-based information; latents' the sandwich.")
  }
  rbind(do.call(rbind, cells), test)
}

#' Every listing of every case against latents' tables
distal_comparison <- function(references, returned) {
  do.call(rbind, lapply(names(references), function(case) {
    do.call(rbind, lapply(names(distal_methods), function(method) {
      rbind(distal_compare_listing(case, method, FALSE, references[[case]], returned),
            distal_compare_listing(case, method, TRUE, references[[case]], returned))
    }))
  }))
}
