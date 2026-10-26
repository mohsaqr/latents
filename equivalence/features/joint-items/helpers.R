# Fitting the joint-items kit with latents and reading what Latent GOLD and
# Mplus returned for it. Sourced by equivalence/test-features-joint-items.R
# (working directory: equivalence/), which also sees the shared equivalence
# helpers (read_mplus_lines, printed_value, montecarlo_relabel).

ji_root <- file.path("features", "joint-items")
ji_cases <- function() readRDS(file.path(ji_root, "kit", "cases.rds"))

#' Explicit precision controls for independent numerical verification
#' @param tol,max_iter Original fitting controls.
#' @return Supplied-start controls, never looser than the original controls.
ji_verification_controls <- function(tol = 1e-10, max_iter = 5000L) {
  list(n_starts = 1L, tol = min(tol, 1e-13), max_iter = max(max_iter, 5000L))
}

#' Fit one seed with the declared controls
#' @param case Kit case.
#' @param seed,n_starts,tol,max_iter Fitting controls.
#' @param local_dependence Whether to join the case's pairs.
#' @param start Optional supplied parameters.
#' @return A fitted model, without diagnostic bookkeeping.
ji_fit_once <- function(case, seed, n_starts, tol, max_iter,
                        local_dependence = TRUE, start = NULL) {
  spec <- case$spec
  common <- list(n_starts = n_starts, tol = tol, max_iter = max_iter, seed = seed,
                 missing = if (isTRUE(spec$missing)) "fiml" else "error",
                 local_dependence = if (local_dependence) spec$pairs else character(),
                 weights = spec$weights, start = start)
  common <- Filter(Negate(is.null), common)
  if (is.null(spec$h)) {
    do.call(lca, c(list(case$data, spec$vars, spec$k), common))
  } else {
    do.call(multilca, c(list(case$data, spec$vars, spec$id, n_profiles = spec$k,
                             n_group_classes = spec$h), common))
  }
}

#' Resume only the fit's own selected parameters, without reference inputs
#' @param case,fit The case and its original fitted model.
#' @param seed,tol,max_iter Original fitting controls.
#' @param local_dependence Whether the original fit joined the case's pairs.
#' @return The fit after a strict supplied-start verification run.
ji_refine <- function(case, fit, seed, tol = 1e-10, max_iter = 5000L,
                      local_dependence = TRUE) {
  controls <- ji_verification_controls(tol, max_iter)
  ji_fit_once(case, seed, controls$n_starts, controls$tol, controls$max_iter,
              local_dependence, starting_values(fit))
}

#' Fit every original seed, then strictly verify each selected fit
#'
#' Production convergence is a per-step gain condition, not a bound on the
#' remaining objective gap. The external comparison uses the declared strict
#' resume recipe and still checks its gap against the independent reference.
#' Every original fit and its original controls remain available separately.
ji_fit <- function(case, n_starts = 15L, seeds = 1:2, local_dependence = TRUE, tol = 1e-10,
                   max_iter = 5000) {
  caught <- character()
  capture <- function(expression) {
    withCallingHandlers(expression, warning = function(w) {
      caught <<- c(caught, class(w)[1L])
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  }
  original <- lapply(seeds, function(seed) capture(
    ji_fit_once(case, seed, n_starts, tol, max_iter, local_dependence)))
  fits <- Map(function(seed, fit) capture(
    ji_refine(case, fit, seed, tol, max_iter, local_dependence)), seeds, original)
  maxima <- vapply(fits, `[[`, numeric(1), "log_likelihood")
  best <- fits[[which.max(maxima)]]
  best$ji_seed_maxima <- maxima
  best$ji_original_fits <- original
  best$ji_verified_fits <- fits
  best$ji_seeds <- seeds
  best$ji_original_controls <- list(n_starts = n_starts, tol = tol, max_iter = max_iter)
  best$ji_verification_controls <- ji_verification_controls(tol, max_iter)
  best$ji_warnings <- unique(caught)
  best
}

#' Original and verified seed fits against the selected fit's reference polish
#' @param case,fit The kit case and its configured fit.
#' @param polish The independent polish from `ji_polish()` (tighter production
#'   reference for the two-level cases, which have no independent group EM).
#' @return One diagnostic row per original seed; no default reliability claim.
ji_verification_diagnostics <- function(case, fit, polish) {
  reference <- polish$reference
  gap <- function(candidate) {
    montecarlo_relabel(candidate$subject_posteriors, reference$posteriors)$gap
  }
  do.call(rbind, lapply(seq_along(fit$ji_seeds), function(index) {
    original <- fit$ji_original_fits[[index]]
    verified <- fit$ji_verified_fits[[index]]
    data.frame(case = case$name, seed = fit$ji_seeds[index],
      original_n_starts = fit$ji_original_controls$n_starts,
      original_tol = fit$ji_original_controls$tol,
      original_max_iter = fit$ji_original_controls$max_iter,
      original_log_likelihood = original$log_likelihood,
      original_converged = original$converged,
      original_gap_to_reference = reference$log_likelihood - original$log_likelihood,
      original_posterior_gap = gap(original),
      verification_n_starts = fit$ji_verification_controls$n_starts,
      verification_tol = fit$ji_verification_controls$tol,
      verification_max_iter = fit$ji_verification_controls$max_iter,
      verification_log_likelihood = verified$log_likelihood,
      verification_converged = verified$converged,
      verification_gap_to_reference = reference$log_likelihood - verified$log_likelihood,
      verification_posterior_gap = gap(verified), reference_method = reference$method)
  }))
}

#' Apply a function over a list on two cores (the machine is shared)
#'
#' Forked workers where they are safe; otherwise a socket cluster whose
#' workers load this source tree and the joint-items files. A worker that
#' fails stops the run with its message. Set `latents.joint_items.serial`
#' to `TRUE` with `options()` to run without worker processes.
ji_parallel <- function(x, fun, ...) {
  results <- if (isTRUE(getOption("latents.joint_items.serial", FALSE))) {
    lapply(x, fun, ...)
  } else if (latents:::.latents_fork_safe()) {
    parallel::mclapply(x, fun, ..., mc.cores = 2L, mc.preschedule = FALSE)
  } else {
    cluster <- parallel::makePSOCKcluster(2L)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    setup <- function(root, files) {
      suppressMessages(pkgload::load_all(root, quiet = TRUE))
      invisible(lapply(files, sys.source, envir = globalenv()))
      NULL
    }
    environment(setup) <- globalenv()
    # The workers hold every joint-items function in their global environment;
    # sending `fun` with its own would serialize the caller's whole frame.
    environment(fun) <- globalenv()
    shared <- list.files(".", "^helper-.*[.]R$", full.names = TRUE)
    own <- file.path(ji_root, c("helpers.R", "reference.R", "compare.R", "generators.R"))
    parallel::clusterCall(cluster, setup, normalizePath(".."),
                          normalizePath(c(shared, own)))
    parallel::parLapply(cluster, x, fun, ...)
  }
  failed <- vapply(results, \(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(failed)) {
    stop(sprintf("%d of %d workers failed: %s", sum(failed), length(results),
                 paste(names(results)[failed], vapply(results[failed], \(r)
                   paste(format(r), collapse = " "), ""), collapse = "; ")), call. = FALSE)
  }
  results
}

#' Fit every case
ji_fit_all <- function(cases, ...) ji_parallel(cases, function(case, ...) ji_fit(case, ...), ...)

# -- Latent GOLD -------------------------------------------------------------------

ji_lg_lines <- function(run) {
  path <- file.path(ji_root, "returned", paste0(run, ".lst"))
  sub("\r$", "", iconv(readLines(path, warn = FALSE), "latin1", "UTF-8"))
}

#' Probability tables of the EstimatedValues-Model section
#'
#' One block per dependent variable: a header row starting with "Cluster",
#' preceded by the item's name (one item) or by the first item's name, its
#' category row and the second item's name (a pair). Each class row holds
#' value / s.e. pairs.
#' @return Named list of k x C matrices; a pair is named "a:b" and its
#'   columns "category_a:category_b".
ji_lg_estimated <- function(lines, k) {
  start <- grep("^EstimatedValues-Model", lines)[1L]
  stopifnot("no EstimatedValues-Model section" = !is.na(start))
  section <- lines[start:length(lines)]
  headers <- which(startsWith(section, "Cluster\t"))
  tables <- lapply(headers, function(at) {
    header <- strsplit(section[at], "\t", fixed = TRUE)[[1L]]
    value_columns <- which(seq_along(header) > 1L & header != "s.e." & nzchar(trimws(header)))
    rows <- lapply(section[at + seq_len(k)], \(row) strsplit(row, "\t", fixed = TRUE)[[1L]])
    stopifnot("class rows are not numbered 1..k" =
                identical(vapply(rows, `[`, "", 1L), as.character(seq_len(k))))
    values <- t(vapply(rows, \(row) as.numeric(row[value_columns]), numeric(length(value_columns))))
    second <- header[value_columns]
    name_line <- \(i) trimws(strsplit(section[i], "\t", fixed = TRUE)[[1L]][2L])
    # A pair's block has a row of the first item's categories two lines up;
    # a single item's block has the blank separator there.
    categories <- strsplit(section[at - 2L], "\t", fixed = TRUE)[[1L]]
    if (sum(nzchar(trimws(categories))) >= 2L) {
      colnames(values) <- paste(trimws(categories[value_columns]), second, sep = ":")
      name <- paste(name_line(at - 3L), name_line(at - 1L), sep = ":")
    } else {
      colnames(values) <- second
      name <- name_line(at - 1L)
    }
    list(name = name, values = values)
  })
  stats::setNames(lapply(tables, `[[`, "values"), vapply(tables, `[[`, "", "name"))
}

#' Class sizes from the Profile section's Size row (Cluster columns)
ji_lg_sizes <- function(lines, k) {
  start <- grep("^Profile", lines)[1L]
  header <- strsplit(lines[start + 1L], "\t", fixed = TRUE)[[1L]]
  first <- which(header == "Cluster")[1L]
  size <- strsplit(grep("^Size\t", lines, value = TRUE)[1L], "\t", fixed = TRUE)[[1L]]
  as.numeric(size[first + 2L * (seq_len(k) - 1L)])
}

#' Everything compared from one Latent GOLD run
ji_latentgold <- function(run, k) {
  lines <- ji_lg_lines(run)
  saved <- utils::read.delim(file.path(ji_root, "returned", paste0(run, "_posteriors.txt")),
                             check.names = FALSE, colClasses = "character")
  list(log_likelihood = printed_value(lines, "^Log-likelihood \\(LL\\)"),
       n_parameters = as.integer(printed_value(lines, "^Number of parameters \\(Npar\\)")),
       bic = printed_value(lines, "^BIC \\(based on LL\\)"),
       n_cases = printed_value(lines, "^Number of cases"),
       sizes = ji_lg_sizes(lines, k),
       estimated = ji_lg_estimated(lines, k),
       id = as.integer(saved[["id"]]),
       posteriors = matrix(as.numeric(as.matrix(saved[paste0("Cluster#", seq_len(k))])),
                           ncol = k))
}

# -- Mplus -------------------------------------------------------------------------

#' Everything compared from one Mplus run
ji_mplus <- function(name, k) {
  path <- file.path(ji_root, "returned", paste0(name, ".out"))
  lines <- read_mplus_lines(path)
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(ji_root, "returned", paste0(name, "-results.dat")), quiet = TRUE)
  stopifnot("results file does not end with the fit block" =
              isTRUE(all.equal(values[2L * npar + 1L], npar)))
  saved <- utils::read.table(file.path(ji_root, "returned", paste0(name, "-posteriors.dat")),
                             na.strings = "*")
  columns <- montecarlo_save_columns(lines)
  names(saved) <- columns
  mlr <- as.integer(any(grepl("H0 Scaling Correction Factor", lines)))
  list(log_likelihood = values[2L * npar + 2L], n_parameters = npar,
       # MLR saves the scaling correction factor after the log likelihood.
       bic = values[2L * npar + 4L + mlr], id = as.integer(saved$ID),
       posteriors = as.matrix(saved[grep("^CPROB", columns)]),
       terminated = any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
       replicated = any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
}

# -- latents side --------------------------------------------------------------------

#' The fit's quantities in one place, classes in a given order
ji_latents <- function(fit, order = seq_len(fit$n_profiles)) {
  two <- !isTRUE(fit$single_level)
  sizes <- if (two) c(fit$group_probabilities %*% fit$profile_probabilities) else
    c(fit$profile_probabilities[1L, ])
  # Response matrices retain the full grid's canonical labels. The public
  # dependence table can label unobserved combinations NA:NA; its independent
  # recomputation remains checked separately in ji_study_one().
  joint <- lapply(fit$response_probabilities[fit$local_dependence$name],
                   \(m) m[order, , drop = FALSE])
  items <- lapply(fit$response_probabilities[setdiff(names(fit$response_probabilities),
                                                     names(joint))],
                  \(m) m[order, , drop = FALSE])
  list(log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters, bic = fit$bic,
       sizes = sizes[order], items = items, joint = joint,
       posteriors = fit$subject_posteriors[, order, drop = FALSE])
}

#' Unobserved combinations of each pair's observed categories
#'
#' Latent GOLD's `a <-> b | Cluster` has C1 C2 - 1 free probabilities per
#' class whether or not every combination is observed, as does latents.
#' Mplus and poLCA are given only the observed combinations in this kit.
ji_unobserved_cells <- function(case) {
  sum(vapply(strsplit(case$spec$pairs, ":", fixed = TRUE), function(p) {
    a <- case$data[[p[1]]]
    b <- case$data[[p[2]]]
    both <- !is.na(a) & !is.na(b)
    length(unique(a[!is.na(a)])) * length(unique(b[!is.na(b)])) -
      nrow(unique(data.frame(a[both], b[both])))
  }, numeric(1)))
}

#' Extract a set of joint tables' probabilities for unobserved combinations
#' @param case A kit case, including its original data and pair specification.
#' @param probabilities Named joint probability matrices, with category labels.
#' @return Probabilities of unobserved cells present in the supplied tables.
ji_unobserved_probabilities <- function(case, probabilities) {
  stopifnot("case and probability tables must be lists" =
              is.list(case) && is.list(probabilities))
  observed <- ref_code(case$data, case$spec$vars, case$spec$pairs)$levels
  unlist(lapply(case$spec$pairs, function(pair) {
    table <- probabilities[[pair]]
    table[, setdiff(colnames(table), observed[[pair]]), drop = FALSE]
  }), use.names = FALSE)
}

#' One comparison row
ji_row <- function(case, reference, quantity, latents, value, tolerance) {
  latents <- as.numeric(latents)
  value <- as.numeric(value)
  stopifnot("compared vectors differ in length" = length(latents) == length(value))
  data.frame(case = case, reference = reference, quantity = quantity,
             n = length(value),
             # An empty comparison (no unobserved combinations) differs by nothing.
             difference = if (length(value) == 0L) 0 else max(abs(latents - value)),
             tolerance = tolerance, stringsAsFactors = FALSE)
}
