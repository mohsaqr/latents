# Genuine Mplus mixture runs kept in the sibling repository JStats
# (validation/r-reference/), refitted with this package. The outputs and their
# data files are byte-for-byte copies in jstats-mplus/; PROVENANCE.md there
# records their source paths, md5 sums and the JStats commit.
#
# Every model is read from the run's own input echo (VARIABLE, ANALYSIS and
# MODEL; the free-text TITLE is ignored), so the data columns, indicators,
# measurement types, class predictors and auxiliary covariates come from the
# input Mplus ran, not from a copy typed into this file. Only the covariance
# structure of a profile model is stated below, and it is checked against the
# input's class-specific statements.
#
# Mplus CATEGORICAL items are fitted as unordered categorical items. With
# thresholds free in every class (no equality labels, as in every run here) a
# K-category item has K - 1 free thresholds per class, and ordered free
# thresholds are a one-to-one reparameterization of an unrestricted K-category
# probability vector, so the two likelihoods are the same function. For binary
# items this is the familiar identity; for polytomous items it holds for the
# same reason. NOMINAL is the multinomial logit, unrestricted per class.
#
# Tolerances. Mplus prints the log likelihood to three decimals (5e-4) and
# stops on its own convergence criteria, which leave up to about 1e-4 of the
# maximum; BIC is a linear function of it (5e-3). Entropy and AvePP are
# printed to three decimals (5e-4) and computed from Mplus's own posteriors,
# which sit up to 2e-3 from ours when Mplus stops short of the exact maximum;
# when Mplus saved the posteriors the measured difference replaces the 2e-3.
# Logit coefficients are printed to three decimals; 1e-3 more allows for where
# Mplus's convergence criterion stopped the parameters.

jm_settings <- list(n_starts = 50L, seed = 1L, tol = 1e-10, max_iter = 5000L)
jm_ll_tol <- mplus_printed_tolerance + 1e-4
jm_bic_tol <- 5e-3
jm_logit_tol <- mplus_printed_tolerance + 1e-3

jm_path <- function(...) file.path("jstats-mplus", ...)

# --- Reading an Mplus input echo ----------------------------------------------

#' The sections of the input echoed at the top of an Mplus output
#'
#' A section header sits at the start of a line (Mplus echoes it with at most
#' three spaces); everything up to the next header belongs to it. TITLE is
#' dropped: it is free text and may contain BY, ON, CENSORED or even
#' `DATA: FILE IS`.
#' @return A named list of section texts, names upper case (`"MODEL"`,
#'   `"MODEL C1"`, ...).
jm_sections <- function(lines) {
  start <- grep("^INPUT INSTRUCTIONS", lines)[1L]
  # The echo ends where Mplus starts talking: its input warnings, or the
  # line saying the input was read.
  end <- grep("^(INPUT READING TERMINATED NORMALLY|\\*\\*\\* WARNING)", lines)
  end <- end[end > start][1L]
  stopifnot("no input echo in the output" = !is.na(start) && !is.na(end))
  echo <- lines[seq.int(start + 1L, end - 1L)]
  echo <- echo[nzchar(trimws(echo))]
  header <- paste0("^ {0,3}(TITLE|DATA|VARIABLE|DEFINE|ANALYSIS|",
                   "MODEL[A-Za-z0-9 ]*|OUTPUT|SAVEDATA|PLOT|MONTECARLO)\\s*:(.*)$")
  is_header <- grepl(header, echo, ignore.case = TRUE)
  stopifnot("the echo does not start with a section" = isTRUE(is_header[1L]))
  section_id <- cumsum(is_header)
  names_of <- toupper(trimws(gsub("\\s+", " ",
    sub(header, "\\1", echo[is_header], ignore.case = TRUE))))
  bodies <- vapply(split(echo, section_id), \(chunk) {
    chunk[1L] <- sub(header, "\\2", chunk[1L], ignore.case = TRUE)
    paste(chunk, collapse = " ")
  }, character(1))
  sections <- stats::setNames(as.list(bodies), names_of)
  sections[names(sections) != "TITLE"]
}

#' The statements of a section, whitespace squeezed
jm_statements <- function(text) {
  if (is.null(text)) return(character())
  statements <- trimws(gsub("\\s+", " ", strsplit(text, ";", fixed = TRUE)[[1L]]))
  statements[nzchar(statements)]
}

#' The value of an option such as `NAMES ARE ...`, or NULL when absent
jm_option <- function(statements, key) {
  pattern <- sprintf("^%s\\b\\s*(ARE|IS|=)?\\s*", key)
  hit <- grep(pattern, statements, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0L) return(NULL)
  stopifnot("an option is given twice" = length(hit) == 1L)
  sub(pattern, "", hit, ignore.case = TRUE)
}

#' Statements without the class labels (`%OVERALL%`, `%c#1%`) that open them
jm_strip_labels <- function(statements) sub("^(%[^%]*%\\s*)+", "", statements)

#' Expand an Mplus variable list (`y1-y4 x`) to lower-case names
#'
#' A range whose ends are both in `known` is taken in that list's order, as
#' Mplus does for USEVARIABLES; otherwise it is a numbered range.
jm_expand <- function(value, known = NULL) {
  if (length(value) == 0L) return(character())
  stopifnot("one variable list at a time" = length(value) == 1L)
  tokens <- strsplit(tolower(trimws(value)), "\\s+")[[1L]]
  unlist(lapply(tokens, \(token) {
    ends <- strsplit(token, "-", fixed = TRUE)[[1L]]
    if (length(ends) == 1L) return(token)
    if (all(ends %in% known)) {
      at <- match(ends, known)
      return(known[seq.int(at[1L], at[2L])])
    }
    prefix <- sub("[0-9]+$", "", ends)
    stopifnot("a numbered range needs one prefix" = identical(prefix[1L], prefix[2L]))
    numbers <- as.integer(sub("^[^0-9]*", "", ends))
    paste0(prefix[1L], seq.int(numbers[1L], numbers[2L]))
  }), use.names = FALSE)
}

#' The model a run specifies, read from its input echo
#' @return A list: `names`, `analysis` (variables Mplus analysed), `indicators`,
#'   `covariates` (`c ON`), `auxiliary` and `auxiliary_method`, the
#'   measurement type sets, `k`, `class_varying_variances`, and `missing`.
jm_input <- function(lines) {
  sections <- jm_sections(lines)
  variable <- jm_statements(sections$VARIABLE)
  model <- jm_statements(sections$MODEL)
  names <- jm_expand(jm_option(variable, "NAMES"))
  auxiliary_text <- jm_option(variable, "AUXILIARY") %||% ""
  auxiliary_method <- unique(toupper(regmatches(
    auxiliary_text, gregexpr("(?<=\\()[A-Za-z0-9]+(?=\\))", auxiliary_text,
                             perl = TRUE))[[1L]]))
  auxiliary <- jm_expand(trimws(gsub("\\([^)]*\\)", " ", auxiliary_text)), names)
  analysis <- jm_expand(jm_option(variable, "USEVARIABLES"), names)
  if (length(analysis) == 0L) analysis <- setdiff(names, auxiliary)
  classes <- jm_option(variable, "CLASSES")
  class_name <- tolower(trimws(sub("\\(.*$", "", classes)))
  stopifnot("one latent class variable" = !grepl("\\)\\s*\\S", classes))
  # A statement may carry the class labels that open it (`%OVERALL% c ON x`,
  # or `%OVERALL% %c#1% [y1]` when %OVERALL% is empty).
  bare <- jm_strip_labels(model)
  on <- grep(sprintf("^%s ON ", class_name), bare, ignore.case = TRUE, value = TRUE)
  covariates <- jm_expand(sub("^\\S+ ON ", "", on, ignore.case = TRUE), names)
  # Statements after the first class-specific label that are bare variable
  # lists free the variances in that class.
  class_specific <- cumsum(grepl("%[^%]*#[^%]*%", model)) > 0L
  variance_lists <- class_specific & nzchar(bare) &
    !grepl("^\\[|\\b(ON|BY|WITH)\\b|[|@$]", bare, ignore.case = TRUE)
  list(names = names, analysis = analysis,
       indicators = setdiff(analysis, covariates), covariates = covariates,
       auxiliary = auxiliary, auxiliary_method = auxiliary_method,
       categorical = jm_expand(jm_option(variable, "CATEGORICAL"), names),
       nominal = jm_expand(jm_option(variable, "NOMINAL"), names),
       count = jm_expand(jm_option(variable, "COUNT"), names),
       k = as.integer(sub(".*\\(\\s*([0-9]+)\\s*\\).*", "\\1", classes)),
       class_varying_variances = any(variance_lists),
       missing = jm_option(variable, "MISSING"),
       data_file = sub("^\\s*FILE\\s*(IS|=)?\\s*\"?([^\"]+?)\"?\\s*$", "\\2",
                       jm_statements(sections$DATA)[1L], ignore.case = TRUE))
}

#' The run's data, columns named as in NAMES, Mplus missing codes as NA
jm_read_data <- function(input) {
  na_strings <- "NA"
  missing <- input$missing
  if (!is.null(missing)) {
    stopifnot("only MISSING ARE ALL (codes) or a symbol is handled" =
                grepl("^(ALL\\s*\\(.*\\)|[*.])$", missing, ignore.case = TRUE))
    na_strings <- if (grepl("^ALL", missing, ignore.case = TRUE)) {
      strsplit(trimws(sub("^ALL\\s*\\((.*)\\)$", "\\1", missing,
                          ignore.case = TRUE)), "[ ,]+")[[1L]]
    } else missing
  }
  data <- utils::read.table(jm_path(input$data_file), header = FALSE,
                            na.strings = na_strings)
  stopifnot("the data file has one column per NAMES entry" =
              ncol(data) == length(input$names))
  names(data) <- input$names
  data
}

# --- Reading Mplus results ------------------------------------------------------

#' Logit coefficients of a `C#k ON` / `Intercepts` block starting after `from`
#' @return A data frame: `class` (Mplus label), `term` (lower case, the
#'   intercept as `(Intercept)`), `estimate`; each against Mplus's last class.
jm_logits_after <- function(lines, from) {
  body <- lines[seq.int(from + 1L, length(lines))]
  stop_at <- which(grepl("^\\S", body))[1L]
  body <- body[seq_len(stop_at - 1L)]
  heading <- ifelse(grepl("^\\s*C#[0-9]+\\s+ON\\s*$", body), "on",
              ifelse(grepl("^\\s*Intercepts\\s*$", body), "intercept", NA))
  entry <- grepl("^\\s+\\S+\\s+-?[0-9]+\\.[0-9]+(\\s+\\S+){3}\\s*$", body)
  block <- cumsum(!is.na(heading))
  heads <- body[!is.na(heading)]
  rows <- lapply(which(entry), \(i) {
    fields <- strsplit(trimws(body[i]), "\\s+")[[1L]]
    head_line <- heads[block[i]]
    if (grepl("Intercepts", head_line)) {
      data.frame(class = as.integer(sub("C#", "", fields[1L])),
                 term = "(Intercept)", estimate = as.numeric(fields[2L]))
    } else {
      data.frame(class = as.integer(sub("^\\s*C#([0-9]+).*$", "\\1", head_line)),
                 term = tolower(fields[1L]), estimate = as.numeric(fields[2L]))
    }
  })
  result <- do.call(rbind, rows)
  stopifnot("no logit coefficients after the anchor" = NROW(result) > 0L)
  result
}

#' The one-step membership logits of MODEL RESULTS
jm_model_logits <- function(lines) {
  results <- grep("^MODEL RESULTS", lines)[1L]
  anchor <- grep("^Categorical Latent Variables\\s*$", lines)
  jm_logits_after(lines, anchor[anchor > results][1L])
}

#' The R3STEP logits (first parameterization, against the last class)
jm_r3step_logits <- function(lines) {
  anchor <- grep("^THE 3-STEP \\(FIXED LOGITS\\) PROCEDURE", lines)[1L]
  stopifnot("no R3STEP block" = !is.na(anchor))
  jm_logits_after(lines, anchor)
}

#' Mplus's fit summary
jm_mplus_fit <- function(lines) {
  list(log_likelihood = printed_value(lines, "^\\s+H0 Value\\s"),
       n_parameters = printed_value(lines, "^Number of Free Parameters"),
       bic = printed_value(lines, "^\\s+Bayesian \\(BIC\\)"))
}

# --- Matching classes -------------------------------------------------------------

jm_permutations <- function(k) {
  if (k <= 1L) return(list(seq_len(k)))
  unlist(lapply(seq_len(k), \(first)
    lapply(jm_permutations(k - 1L), \(rest) c(first, setdiff(seq_len(k), first)[rest]))),
    recursive = FALSE)
}

#' This package's class for each Mplus class, from the posterior class counts
#' @return A list: `order` (fit class of each Mplus class), `gap` (largest
#'   count difference under it) and `runner_up` (under the next best order).
jm_match_classes <- function(mplus_counts, fit_counts) {
  candidates <- jm_permutations(length(mplus_counts))
  gaps <- vapply(candidates, \(p) max(abs(fit_counts[p] - mplus_counts)), numeric(1))
  ranked <- sort(gaps)
  list(order = candidates[[which.min(gaps)]], gap = ranked[1L],
       runner_up = if (length(ranked) > 1L) ranked[2L] else Inf)
}

#' Fit logits of every class against the class matched to Mplus's last one
#' @return A data frame laid out as [jm_logits_after()]'s.
jm_fit_logits <- function(coefficients, order) {
  full <- cbind(coefficients, 0)
  k <- length(order)
  reference <- full[, order[k]]
  do.call(rbind, lapply(seq_len(k - 1L), \(j)
    data.frame(class = j, term = rownames(full),
               estimate = unname(full[, order[j]] - reference))))
}

jm_compare_logits <- function(theirs, mine, run) {
  joined <- merge(theirs, mine, by = c("class", "term"), suffixes = c("_mplus", "_fit"))
  expect_identical(nrow(joined), nrow(theirs), label = paste(run, "logit terms"))
  expect_identical(nrow(joined), nrow(mine), label = paste(run, "logit terms"))
  expect_lt(max(abs(joined$estimate_mplus - joined$estimate_fit)), jm_logit_tol,
            label = paste(run, "logits"))
}

# --- The runs ---------------------------------------------------------------------

# `model`: "lca" for discrete indicators, otherwise the covariance structure of
# the profile model. Mplus's default (no class-specific variance statements)
# is EEI; variances freed in every class make it VVI.
jm_runs <- data.frame(
  run = c(
    "lpa-constraint-new", "lpa-cov-tech14", "lpa-eei-k3-tech11",
    "lpa-ex710-tech14", "lpa-k3-tech11", "lpa-k3-vvi-tech14", "lpa-tech14",
    "lca-entropy", "lca-k1starts-tech11", "lca-k3-tech10", "lca-k3-tech11",
    "lca-k3-tech14", "lca-k4-tech11", "lca-k5-tech11", "lca-tech10",
    "lca-tech14", "lca-tech14-lrtstarts", "lca-tech14-tech11",
    "miterations-uncapped", "savedata-cprob", "training-partial-off",
    "lca-4cat-tech10", "lca-mixed-tech10", "lca-mixed-tech14", "lca-poly-tech10",
    "lca-nominal-tech11", "lca-nominal-k3-tech11", "lca-count-tech11",
    "lca-cov-entropy", "lca-cov-tech10", "lca-covariate-tech11",
    "lca-cov-k3-tech14",
    "lca-du3step-multi", "lca-de3step-multi", "lca-r3step-joint",
    "ug-ex7.3-de3step", "ug-ex7.3-du3step", "ug-ex7.3-r3step",
    "ug-ex7.3-two-step-du3step", "ug-ex7.3-two-step-du3step-both"),
  model = c("EEI", "EEI", "EEI", "VVI", "VVI", "VVI", "EEI", rep("lca", 33L))
)
jm_runs$file <- ifelse(startsWith(jm_runs$run, "ug-"),
                       paste0("mplus-", jm_runs$run, ".out.txt"),
                       paste0("mplus-demo-", jm_runs$run, ".out.txt"))

jm_fit_run <- function(model, input, data) {
  if (identical(model, "lca")) {
    lca(data, input$indicators, input$k, count = input$count,
        profile_covariates = input$covariates,
        n_starts = jm_settings$n_starts, seed = jm_settings$seed,
        tol = jm_settings$tol, max_iter = jm_settings$max_iter)
  } else {
    lpa(data, input$indicators, input$k, model = model,
        profile_covariates = input$covariates,
        n_starts = jm_settings$n_starts, seed = jm_settings$seed,
        tol = jm_settings$tol, max_iter = jm_settings$max_iter)
  }
}

#' Fit one run and keep what the comparisons need
jm_prepare <- function(run, model, file) {
  lines <- read_mplus_lines(jm_path(file))
  input <- jm_input(lines)
  data <- jm_read_data(input)
  fit <- jm_fit_run(model, input, data)
  counts <- mplus_posterior_counts(lines, input$k)
  list(run = run, model = model, lines = lines, input = input, data = data,
       fit = fit, match = jm_match_classes(counts, colSums(fit$subject_posteriors)))
}

jm_check_specification <- function(case) {
  input <- case$input
  run <- case$run
  discrete <- union(union(input$categorical, input$nominal), input$count)
  if (identical(case$model, "lca")) {
    expect_setequal(discrete, input$indicators)
  } else {
    expect_identical(length(discrete), 0L, label = paste(run, "discrete items"))
    expect_identical(input$class_varying_variances, identical(case$model, "VVI"),
                     label = paste(run, "class-varying variances"))
  }
}

jm_cases <- Map(jm_prepare, jm_runs$run, jm_runs$model, jm_runs$file)

invisible(lapply(jm_cases, \(case) {
  test_that(sprintf("%s: latents reaches Mplus's maximum", case$run), {
    jm_check_specification(case)
    expected <- jm_mplus_fit(case$lines)
    fit <- case$fit
    expect_true(fit$converged, label = paste(case$run, "converged"))
    expect_identical(as.numeric(fit$n_parameters), expected$n_parameters,
                     label = paste(case$run, "free parameters"))
    expect_lt(abs(fit$log_likelihood - expected$log_likelihood), jm_ll_tol,
              label = paste(case$run, "log likelihood"))
    # Never a lower maximum than Mplus's.
    expect_gte(fit$log_likelihood, expected$log_likelihood - jm_ll_tol,
               label = paste(case$run, "log likelihood floor"))
    expect_lt(abs(fit$bic_individual - expected$bic), jm_bic_tol,
              label = paste(case$run, "BIC"))
  })

  test_that(sprintf("%s: classification quality matches Mplus", case$run), {
    lines <- case$lines
    quality <- grep("^CLASSIFICATION QUALITY", lines)[1L]
    skip_if(is.na(quality), "Mplus printed no classification quality")
    k <- case$input$k
    # The class match is unambiguous: under it the posterior class counts
    # differ by less than 1e-4 per observation on average (Mplus stops short
    # of the exact maximum), and far less than under any other labelling.
    expect_lt(case$match$gap / nrow(case$data), 1e-4,
              label = paste(case$run, "class counts per observation"))
    expect_gt(case$match$runner_up, 10 * case$match$gap)
    theirs <- list(entropy = mplus_entropy_after(lines, quality),
                   avepp = mplus_avepp_after(lines, quality, k))
    mine <- latents_classification(case$fit, "individuals")
    gap <- classification_discrepancy(theirs, mine, case$match$order)
    expect_lt(gap[["entropy"]], mplus_printed_tolerance + 2e-3,
              label = paste(case$run, "entropy"))
    expect_lt(gap[["avepp"]], mplus_printed_tolerance + 2e-3,
              label = paste(case$run, "AvePP"))
  })

  if (length(case$input$covariates) > 0L) {
    test_that(sprintf("%s: membership logits match Mplus", case$run), {
      jm_compare_logits(jm_model_logits(case$lines),
                        jm_fit_logits(case$fit$profile_coefficients, case$match$order),
                        case$run)
    })
  }
}))

test_that("savedata-cprob: posteriors match the CPROB file Mplus saved", {
  case <- jm_cases[["savedata-cprob"]]
  saved <- utils::read.table(jm_path("mplus-demo-savedata-cprob-out.dat"),
                             col.names = c(paste0("u", 1:5), "cprob1", "cprob2", "mlcc"))
  expect_identical(nrow(saved), nrow(case$data))
  # The saved file carries the indicators too: it is the same rows in order.
  expect_identical(as.matrix(saved[paste0("u", 1:5)]),
                   as.matrix(case$data[paste0("u", 1:5)]) + 0)
  theirs <- as.matrix(saved[c("cprob1", "cprob2")])
  mine <- case$fit$subject_posteriors[, case$match$order]
  posterior_gap <- max(abs(theirs - mine))
  expect_lt(posterior_gap, mplus_printed_tolerance + 1e-4)
  quality <- grep("^CLASSIFICATION QUALITY", case$lines)[1L]
  gap <- classification_discrepancy(
    list(entropy = mplus_entropy_after(case$lines, quality),
         avepp = mplus_avepp_after(case$lines, quality, 2L)),
    latents_classification(case$fit, "individuals"), case$match$order)
  expect_lt(gap[["entropy"]], mplus_printed_tolerance + posterior_gap)
  expect_lt(gap[["avepp"]], mplus_printed_tolerance + posterior_gap)
})

# R3STEP: the measurement model above, then r3step() on the auxiliary
# covariates. The class match of the measurement fit carries over.
jm_check_r3step <- function(run) {
  case <- jm_cases[[run]]
  expect_identical(case$input$auxiliary_method, "R3STEP")
  table <- r3step(case$fit, case$data, covariates = case$input$auxiliary)
  k <- case$input$k
  terms <- unique(table$term)
  coefficients <- vapply(sprintf("profile_%d", seq_len(k - 1L)), \(outcome)
    vapply(terms, \(term) table$estimate[table$outcome == outcome & table$term == term],
           numeric(1)), numeric(length(terms)))
  coefficients <- matrix(coefficients, nrow = length(terms),
                         dimnames = list(terms, NULL))
  jm_compare_logits(jm_r3step_logits(case$lines),
                    jm_fit_logits(coefficients, case$match$order), run)
}

test_that("lca-r3step-joint: R3STEP logits match Mplus", {
  jm_check_r3step("lca-r3step-joint")
})

test_that("ug-ex7.3-r3step: R3STEP logits match Mplus", {
  # Regression: the fixed-error likelihood has a near-separation ridge here
  # (step-one classification errors 0.0815 and 0.0048). BFGS from zero settled
  # on it (-log L 81.669) and r3step() stopped with latents_no_converge; Mplus's
  # maximum is C#1 ON X1 = 3.756, intercept -0.899 (-log L 62.606).
  jm_check_r3step("ug-ex7.3-r3step")
})

# Beyond the profile and class models: a mixture regression and a latent
# transition model, specifications stated here and checked against the echo.

test_that("mixreg-classspecific: mixture_regression() reaches Mplus's maximum", {
  lines <- read_mplus_lines(jm_path("mplus-demo-mixreg-classspecific.out.txt"))
  input <- jm_input(lines)
  model <- jm_statements(jm_sections(lines)$MODEL)
  # y ON x with class-specific slopes and intercepts; the residual variance is
  # not mentioned in a class, so it is class-invariant: variance = "equal".
  expect_true(all(c("%OVERALL% y ON x", "%c#1% y ON x*0.8", "%c#2% y ON x*-0.8")
                  %in% model))
  expect_false(input$class_varying_variances)
  data <- jm_read_data(input)
  fit <- mixture_regression(y ~ x, data, n_classes = input$k, variance = "equal",
                            n_starts = jm_settings$n_starts, seed = jm_settings$seed,
                            tol = jm_settings$tol, max_iter = jm_settings$max_iter)
  expected <- jm_mplus_fit(lines)
  log_likelihood <- as.numeric(logLik(fit))
  expect_identical(as.numeric(fit$n_parameters), expected$n_parameters)
  expect_lt(abs(log_likelihood - expected$log_likelihood), jm_ll_tol)
  expect_gte(log_likelihood, expected$log_likelihood - jm_ll_tol)
  expect_lt(abs(stats::BIC(fit) - expected$bic), jm_bic_tol)
})

test_that("lta-classvarying-measurement: lta() reaches Mplus's maximum", {
  lines <- read_mplus_lines(jm_path("mplus-demo-lta-classvarying-measurement.out.txt"))
  sections <- jm_sections(lines)
  variable <- jm_statements(sections$VARIABLE)
  # Two binary waves, each with its own measurement model (MODEL c1 and MODEL
  # c2 free the thresholds separately), and c2 ON c1: lta() with
  # measurement = "occasion".
  expect_identical(jm_option(variable, "CLASSES"), "c1(2) c2(2)")
  expect_true("c2 ON c1" %in% jm_strip_labels(jm_statements(sections$MODEL)))
  expect_true(all(c("MODEL C1", "MODEL C2") %in% names(sections)))
  names <- jm_expand(jm_option(variable, "NAMES"))
  expect_setequal(jm_expand(jm_option(variable, "CATEGORICAL"), names), names)
  wide <- utils::read.table(jm_path("mplus-demo-lta-classvarying-measurement.dat"),
                            col.names = names)
  occasion <- function(time, columns) {
    data.frame(id = seq_len(nrow(wide)), time = time,
               stats::setNames(wide[columns], c("u1", "u2", "u3")))
  }
  long <- rbind(occasion(1L, c("u11", "u12", "u13")),
                occasion(2L, c("u21", "u22", "u23")))
  items <- c("u1", "u2", "u3")
  fit <- lta(long, items, "id", n_profiles = 2L, time = "time",
             categorical = items, measurement = "occasion",
             n_starts = jm_settings$n_starts, seed = jm_settings$seed,
             tol = jm_settings$tol, max_iter = jm_settings$max_iter)
  expected <- jm_mplus_fit(lines)
  expect_true(fit$converged)
  expect_identical(as.numeric(fit$n_parameters), expected$n_parameters)
  expect_lt(abs(fit$log_likelihood - expected$log_likelihood), jm_ll_tol)
  expect_gte(fit$log_likelihood, expected$log_likelihood - jm_ll_tol)
  # Mplus's BIC counts persons, this package's `bic` the independent units.
  expect_lt(abs(fit$bic - expected$bic), jm_bic_tol)
})
