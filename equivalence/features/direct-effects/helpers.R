# Reading the direct-effect kit's outputs and fitting each case with latents.
# Sourced (sys.source) by equivalence/test-features-direct-effects.R and by
# compare.R; needs montecarlo_relabel(), read_mplus_lines() and printed_value()
# from equivalence/helper-montecarlo.R and helper-mplus-classification.R.
#
# Every quantity is named the same way for all three programs:
# "mean:<indicator>:<class>", "variance:<indicator>:<class>" (class "shared"
# under EEI), "effect:<indicator>:<covariate>" and "proportion:<class>", where
# <class> is the reference's class number after alignment.

de_cases <- function(root) readRDS(file.path(root, "kit", "cases.rds"))

#' Fit one case as the programs were asked to, from several seeds; keep the best
de_fit <- function(case, seeds = 1:3, n_starts = 20L) {
  spec <- case$spec
  if (nrow(case$data) >= 10000L) {
    seeds <- seeds[seq_len(min(2L, length(seeds)))]
    n_starts <- 10L
  }
  fits <- lapply(seeds, function(seed) {
    suppressMessages(lpa(case$data, spec$vars, spec$k, model = spec$model,
                         n_starts = n_starts, seed = seed, tol = 1e-10, max_iter = 5000,
                         direct_effects = spec$direct_effects))
  })
  log_likelihoods <- vapply(fits, `[[`, numeric(1), "log_likelihood")
  best <- fits[[which.max(log_likelihoods)]]
  best$seed_log_likelihoods <- log_likelihoods
  best
}

#' Reference parameter names in Latent GOLD's write order (parameters=first,
#' `y <- 1 | Cluster + z`): class logits 2..k, then per indicator its class
#' intercepts and its effects, then the variances (by class under VVI)
de_latentgold_names <- function(spec) {
  k <- spec$k
  location <- unlist(lapply(spec$vars, function(indicator) {
    covariates <- spec$effects$covariate[spec$effects$indicator == indicator]
    c(sprintf("mean:%s:%d", indicator, seq_len(k)),
      sprintf("effect:%s:%s", indicator, covariates))
  }))
  variances <- if (identical(spec$model, "VVI") && k > 1L) {
    as.vector(t(outer(spec$vars, seq_len(k), \(v, j) sprintf("variance:%s:%d", v, j))))
  } else sprintf("variance:%s:shared", spec$vars)
  c(sprintf("logit:%d", seq_len(k)[-1L]), location, variances)
}

#' One row of Latent GOLD's write= file, as numbers
de_write_row <- function(lines, name) {
  fields <- strsplit(grep(sprintf(",%s,", name), lines, value = TRUE, fixed = TRUE)[1L],
                     ",", fixed = TRUE)[[1L]]
  stopifnot("no such row in the write file" = length(fields) > 5L)
  as.numeric(fields[-seq_len(5L)])
}

#' Everything compared from one Latent GOLD run
de_latentgold <- function(name, case, root) {
  spec <- case$spec
  returned <- file.path(root, "returned")
  lines <- gsub("\r", "", iconv(readLines(file.path(returned, paste0(name, ".lst")),
                                          warn = FALSE), "latin1", "UTF-8"))
  written <- gsub("\r", "", readLines(file.path(returned, paste0(name, "_write.txt")),
                                      warn = FALSE))
  saved <- utils::read.delim(file.path(returned, paste0(name, "_posteriors.txt")),
                             check.names = FALSE)
  parameter_names <- de_latentgold_names(spec)
  estimates <- de_write_row(written, "Parameters")
  errors <- de_write_row(written, "Parameters-SE")
  stopifnot("Latent GOLD wrote another number of parameters than the syntax has" =
              length(estimates) == length(parameter_names) &&
              length(errors) == length(parameter_names))
  names(estimates) <- names(errors) <- parameter_names
  profile <- de_write_row(written, "Profile")
  profile_errors <- de_write_row(written, "Profile-SE")
  sizes <- stats::setNames(profile[seq_len(spec$k)], sprintf("proportion:%d", seq_len(spec$k)))
  size_errors <- stats::setNames(profile_errors[seq_len(spec$k)], names(sizes))
  keep <- !startsWith(parameter_names, "logit")
  list(program = "latentgold",
       log_likelihood = de_write_row(written, "Statistics")[[1L]],
       n_parameters = as.integer(printed_value(lines, "^Number of parameters \\(Npar\\)")),
       bic_printed = printed_value(lines, "^BIC \\(based on LL\\)"),
       estimates = c(estimates[keep], sizes), errors = c(errors[keep], size_errors),
       id = as.integer(saved$id),
       posteriors = as.matrix(saved[paste0("Cluster#", seq_len(spec$k))]),
       terminated = TRUE, replicated = NA)
}

#' Parameter numbers of an Mplus TECH1 output: one row per (class, matrix,
#' row, column) holding a free parameter
de_mplus_tech1 <- function(lines) {
  start <- grep("^\\s*TECHNICAL 1 OUTPUT", lines)[1L]
  stop <- grep("^\\s*STARTING VALUES", lines)
  stop <- stop[stop > start][1L]
  section <- lines[start:(stop - 1L)]
  class_at <- grep("PARAMETER SPECIFICATION FOR LATENT CLASS", section)
  matrix_at <- grep("^\\s+(NU|LAMBDA|THETA|ALPHA|BETA|PSI|TAU|GAMMA)(\\(C\\))?\\s*$", section)
  is_number <- function(tokens) grepl("^-?[0-9]+$", tokens)
  rows <- lapply(matrix_at, function(at) {
    block <- utils::tail(class_at[class_at < at], 1L)
    class_label <- sub(".*LATENT CLASS ", "", trimws(section[[block]]))
    ends <- c(matrix_at, class_at, length(section) + 1L)
    end <- min(ends[ends > at]) - 1L
    body <- section[(at + 1L):end]
    body <- body[nzchar(trimws(body)) & !grepl("^\\s*_+", body)]
    columns <- character()
    entries <- list()
    # Wide matrices come in chunks of five columns: a header line names a
    # chunk's columns, the rows that follow fill them.
    invisible(lapply(body, function(line) {
      tokens <- strsplit(trimws(line), "\\s+")[[1L]]
      if (!any(is_number(tokens))) {
        columns <<- tokens
      } else {
        label <- if (is_number(tokens[[1L]])) "" else tokens[[1L]]
        values <- as.integer(tokens[is_number(tokens)])
        entries[[length(entries) + 1L]] <<- data.frame(
          class = class_label, matrix = trimws(section[[at]]), row = label,
          column = columns[seq_along(values)], number = values)
      }
    }))
    do.call(rbind, entries)
  })
  table <- do.call(rbind, rows)
  table[table$number > 0L, , drop = FALSE]
}

#' Mplus parameter numbers in the common naming
de_mplus_names <- function(tech1, spec) {
  in_class <- function(matrix, row, column, class) {
    hit <- tech1$number[tech1$class == as.character(class) & tech1$matrix == matrix &
                          tech1$row == row & tech1$column == column]
    if (length(hit) == 0L) NA_integer_ else hit[[1L]]
  }
  upper <- toupper(spec$vars)
  numbers <- unlist(lapply(seq_len(spec$k), function(j) {
    means <- vapply(upper, \(v) {
      nu <- in_class("NU", "", v, j)
      if (is.na(nu)) in_class("ALPHA", "", v, j) else nu
    }, integer(1))
    variances <- vapply(upper, \(v) {
      theta <- in_class("THETA", v, v, j)
      if (is.na(theta)) in_class("PSI", v, v, j) else theta
    }, integer(1))
    c(stats::setNames(means, sprintf("mean:%s:%d", spec$vars, j)),
      stats::setNames(variances, sprintf("variance:%s:%d", spec$vars, j)))
  }))
  effects <- vapply(seq_len(nrow(spec$effects)), \(e) {
    in_class("BETA", toupper(spec$effects$indicator[[e]]),
             toupper(spec$effects$covariate[[e]]), 1L)
  }, integer(1))
  names(effects) <- sprintf("effect:%s:%s", spec$effects$indicator, spec$effects$covariate)
  logits <- tech1[tech1$matrix == "ALPHA(C)", , drop = FALSE]
  c(numbers, effects,
    stats::setNames(logits$number, sprintf("logit:%s", sub("C#", "", logits$column))))
}

#' Everything compared from one Mplus run
de_mplus <- function(name, case, root) {
  spec <- case$spec
  returned <- file.path(root, "returned")
  path <- file.path(returned, paste0(name, ".out"))
  lines <- read_mplus_lines(path)
  terminated <- any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines))
  if (!terminated) return(list(program = "mplus", terminated = FALSE, lines = lines))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(returned, paste0(name, "-results.dat")), quiet = TRUE)
  stopifnot("results file does not end with the fit block" =
              isTRUE(all.equal(values[2L * npar + 1L], npar)))
  numbers <- de_mplus_names(de_mplus_tech1(lines), spec)
  stopifnot("a parameter is missing from TECH1" = !anyNA(numbers),
            "TECH1 numbers do not cover the free parameters" =
              setequal(unique(numbers), seq_len(npar)))
  estimates <- stats::setNames(values[numbers], names(numbers))
  errors <- stats::setNames(values[npar + numbers], names(numbers))
  # A covariate Mplus read as s z (DEFINE) has the effect b / s.
  rescaled <- startsWith(names(estimates), "effect:")
  estimates[rescaled] <- estimates[rescaled] * spec$mplus_scale
  errors[rescaled] <- errors[rescaled] * spec$mplus_scale
  # EEI: one variance per indicator, repeated in every class block.
  if (identical(spec$model, "EEI") || spec$k == 1L) {
    shared <- grepl("^variance:.*:1$", names(estimates))
    variance_names <- sub(":1$", ":shared", names(estimates)[shared])
    keep <- !grepl("^variance:", names(estimates))
    estimates <- c(estimates[keep], stats::setNames(estimates[shared], variance_names))
    errors <- c(errors[keep], stats::setNames(errors[shared], variance_names))
  }
  # Class proportions from the logits (last class the reference).
  logit_names <- grep("^logit:", names(estimates), value = TRUE)
  logits <- c(estimates[logit_names], 0)
  proportions <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  names(proportions) <- sprintf("proportion:%d", seq_len(spec$k))
  keep <- !startsWith(names(estimates), "logit")
  saved <- utils::read.table(file.path(returned, paste0(name, "-posteriors.dat")),
                             na.strings = "*")
  columns <- montecarlo_save_columns(lines)
  stopifnot("save file does not match its layout" = length(columns) == ncol(saved))
  names(saved) <- columns
  list(program = "mplus", log_likelihood = values[2L * npar + 2L], n_parameters = npar,
       bic_results = values[2L * npar + 4L], estimates = c(estimates[keep], proportions),
       errors = errors[keep], id = as.integer(saved$ID),
       posteriors = if (spec$k == 1L) matrix(1, nrow(saved), 1L) else
         as.matrix(saved[grep("^CPROB", columns)]),
       terminated = TRUE,
       replicated = any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
}

#' A latents fit's quantities in the common naming, classes in the
#' reference's order (`order[j]` is the latents class matching reference
#' class j)
de_latents_values <- function(fit, spec, order) {
  inference <- parameter_inference(fit)
  rank <- match(seq_len(spec$k), order)
  class_of <- function(outcome) rank[as.integer(sub("profile_", "", outcome))]
  means <- subset(inference, parameter == "mean")
  variances <- subset(inference, parameter == "variance")
  effects <- subset(inference, parameter == "direct_effect")
  proportions <- subset(inference, parameter == "probability" & level == "profile")
  variance_class <- ifelse(variances$outcome == "shared", "shared",
                           as.character(class_of(variances$outcome)))
  if (spec$k == 1L || identical(spec$model, "EEI")) variance_class <- "shared"
  keys <- c(sprintf("mean:%s:%d", means$term, class_of(means$outcome)),
            sprintf("variance:%s:%s", variances$term, variance_class),
            sprintf("effect:%s", effects$term),
            sprintf("proportion:%d", class_of(proportions$outcome)))
  estimates <- c(means$estimate, variances$estimate, effects$estimate, proportions$estimate)
  errors <- c(means$standard_error, variances$standard_error, effects$standard_error,
              proportions$standard_error)
  keep <- !duplicated(keys)
  list(estimates = stats::setNames(estimates[keep], keys[keep]),
       errors = stats::setNames(errors[keep], keys[keep]))
}

#' latents posteriors in the rows of a reference's save file (the fit keeps
#' the data's row order; `rows` are the data's ids)
de_posteriors_by_id <- function(fit, id, rows) {
  stopifnot("the reference saved another set of ids" =
              setequal(id, rows) && !anyDuplicated(id) &&
              nrow(fit$subject_posteriors) == length(rows))
  fit$subject_posteriors[match(id, rows), , drop = FALSE]
}

#' Every compared quantity of one case against one reference, one row each
#'
#' Tolerances (each row says which applies):
#' - log likelihood: both references iterate to a relative change below 1e-9
#'   (Mplus RLOGCRITERION, Latent GOLD's NR tolerance 1e-12 is tighter), and
#'   latents to tol 1e-10; 1e-6 relative covers three such stopping errors.
#' - estimates: the optimizers stop at different points of the same flat
#'   maximum. With posterior gap d (largest row difference after alignment)
#'   a class-level quantity moves by at most about d times its scale; we
#'   allow 1e-4 relative plus 10 d (d measured per case, printed in the row).
#' - standard errors: computed at slightly different points by different
#'   numerical Hessians; 1e-3 relative plus 10 d.
#' - posteriors: d itself is reported. Where three weakly separated classes
#'   leave the likelihood flat, rows move by up to 1.04e-4 at likelihoods
#'   equal to 8e-8 (vvi3_03; the Monte Carlo kit measured 1.8e-4), so d must
#'   stay below 5e-4.
#' Mplus's log likelihood differs from Latent GOLD's and latents' at agreeing
#' estimates by up to 3.5e-8 relative (larger at larger n, either sign), well
#' inside the 1e-6.
de_compare <- function(case, fit, reference) {
  spec <- case$spec
  rows <- list()
  add <- function(quantity, item, theirs, mine, tolerance) {
    rows[[length(rows) + 1L]] <<- data.frame(
      case = case$name, config = spec$config, program = reference$program,
      quantity = quantity, item = item, reference = unname(theirs), latents = unname(mine),
      difference = abs(unname(mine) - unname(theirs)), tolerance = unname(tolerance))
  }
  mine_posteriors <- de_posteriors_by_id(fit, reference$id, case$data$id)
  alignment <- montecarlo_relabel(mine_posteriors, reference$posteriors)
  gap <- alignment$gap
  mine <- de_latents_values(fit, spec, alignment$order)
  add("posterior", "max row", 0, gap, 5e-4)
  add("log_likelihood", "", reference$log_likelihood, fit$log_likelihood,
      1e-6 * abs(reference$log_likelihood))
  add("n_parameters", "", reference$n_parameters, fit$n_parameters, 0)
  n <- nrow(case$data)
  bic <- -2 * reference$log_likelihood + reference$n_parameters * log(n)
  add("bic", "from reference LL", bic, fit$bic, 2e-6 * abs(reference$log_likelihood))
  if (!is.null(reference[["bic_results"]])) {
    add("bic", "results file", reference[["bic_results"]], fit$bic,
        2e-6 * abs(reference$log_likelihood))
  }
  if (!is.null(reference[["bic_printed"]])) {
    add("bic", "printed", reference[["bic_printed"]], fit$bic,
        5e-5 + 2e-6 * abs(reference$log_likelihood))
  }
  items <- names(reference$estimates)
  stopifnot("latents lacks a reference quantity" = all(items %in% names(mine$estimates)))
  scale <- pmax(1, abs(reference$estimates))
  kind <- sub(":.*", "", items)
  invisible(lapply(seq_along(items), \(i) {
    add(paste0(kind[[i]], "_estimate"), items[[i]], reference$estimates[[i]],
        mine$estimates[[items[[i]]]], (1e-4 + 10 * gap) * scale[[i]])
  }))
  # Latent GOLD writes -1 for a quantity it does not estimate (the one class
  # of k = 1).
  error_items <- names(reference$errors)[reference$errors >= 0]
  invisible(lapply(error_items, \(item) {
    add(paste0(sub(":.*", "", item), "_se"), item, reference$errors[[item]],
        mine$errors[[item]], (1e-3 + 10 * gap) * abs(reference$errors[[item]]))
  }))
  # predict() on the training data reproduces the fit's posteriors.
  predicted <- predict(fit, newdata = case$data, type = "posterior")
  predicted <- as.matrix(predicted[grep("^probability_profile_", names(predicted))])
  add("predict", "max row vs reference", 0,
      max(abs(predicted[match(reference$id, case$data$id), alignment$order, drop = FALSE] -
                reference$posteriors)), 5e-4)
  do.call(rbind, rows)
}

#' The whole comparison of every case with both programs
de_compare_all <- function(root, cases = de_cases(root)) {
  per_case <- lapply(cases, function(case) {
    fit <- de_fit(case)
    references <- list(
      latentgold = if (file.exists(file.path(root, "returned", paste0(case$name, ".lst"))))
        de_latentgold(case$name, case, root),
      mplus = if (file.exists(file.path(root, "returned", paste0(case$name, ".out"))))
        de_mplus(case$name, case, root))
    references <- Filter(Negate(is.null), references)
    usable <- Filter(\(reference) isTRUE(reference$terminated), references)
    log_likelihoods <- vapply(usable, `[[`, numeric(1), "log_likelihood")
    best <- max(log_likelihoods)
    # Quantities are compared only where a program reached the maximum latents
    # reports (within the log-likelihood tolerance). A program above it is a
    # latents failure, asserted from `maxima`; one below it stopped at a lower
    # maximum, which `maxima` records.
    same <- abs(log_likelihoods - fit$log_likelihood) <= 1e-6 * abs(fit$log_likelihood)
    comparison <- do.call(rbind, lapply(usable[same], \(reference) {
      de_compare(case, fit, reference)
    }))
    maxima <- data.frame(
      case = case$name, config = case$spec$config,
      program = c(names(usable), "latents"),
      log_likelihood = c(log_likelihoods, fit$log_likelihood),
      best_reference = best,
      replicated = c(vapply(usable, \(r) isTRUE(r$replicated), logical(1)), NA),
      compared = c(same, NA),
      seed_spread = c(rep(NA_real_, length(usable)),
                      diff(range(fit$seed_log_likelihoods))))
    maxima$gap_to_best <- maxima$log_likelihood - best
    list(comparison = comparison, maxima = maxima,
         failed = setdiff(names(references), names(usable)), fit = fit)
  })
  list(comparison = do.call(rbind, lapply(per_case, `[[`, "comparison")),
       maxima = do.call(rbind, lapply(per_case, `[[`, "maxima")),
       failed = unlist(lapply(names(per_case), \(name) {
         if (length(per_case[[name]]$failed)) paste(name, per_case[[name]]$failed)
       })),
       fits = lapply(per_case, `[[`, "fit"))
}

#' Evaluate `expr`, muffling only latents' documented boundary warning (a
#' degenerate small-sample fit in the reference battery; the battery records
#' it and skips the standard errors there)
de_quietly <- function(expr) {
  withCallingHandlers(expr, latents_boundary = \(w) invokeRestart("muffleWarning"))
}
