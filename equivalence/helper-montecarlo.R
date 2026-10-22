# Reading the Monte Carlo equivalence outputs (montecarlo/), and fitting each
# case with this package. Every comparison reads the package's reported
# values, aligns classes once from the posteriors and reuses that alignment
# for every class-level quantity, and checks ids for completeness.

montecarlo_dir <- file.path("montecarlo")
montecarlo_cases <- function() readRDS(file.path(montecarlo_dir, "kit", "cases.rds"))

#' Fit one case as Mplus and Latent GOLD were asked to fit it
montecarlo_fit <- function(case, seed = 1L, n_starts = 20L) {
  spec <- case$spec
  data <- case$data
  missing <- if (isTRUE(spec$missing)) "fiml" else "error"
  common <- list(n_starts = n_starts, seed = seed, tol = 1e-10, max_iter = 5000,
                 missing = missing)
  if (identical(spec$level, "transition")) {
    return(do.call(lta, c(list(data, spec$vars, "g", n_profiles = spec$k, time = "t",
                               variance_model = spec$variance), common)))
  }
  if (identical(spec$level, "single")) {
    if (isTRUE(spec$categorical)) {
      return(do.call(lca, c(list(data, spec$vars, n_classes = spec$k), common)))
    }
    return(do.call(lpa, c(list(data, spec$vars, n_profiles = spec$k, model = spec$model,
                               profile_covariates = spec$covariates %||% character()),
                          common)))
  }
  if (isTRUE(spec$categorical)) {
    return(do.call(multilca, c(list(data, spec$vars, "g", n_profiles = spec$k,
                                     n_group_classes = spec$h), common)))
  }
  do.call(multilpa, c(list(data, spec$vars, "g", n_profiles = spec$k,
                           n_group_classes = spec$h, model = spec$model), common))
}

#' One permutation of this package's classes for the whole posterior matrix
montecarlo_relabel <- function(mine, theirs) {
  permutations <- function(k) if (k == 1L) list(1L) else
    do.call(c, lapply(seq_len(k), \(i) lapply(permutations(k - 1L),
                                             \(p) c(i, setdiff(seq_len(k), i)[p]))))
  candidates <- permutations(ncol(mine))
  gaps <- vapply(candidates, \(p) max(abs(mine[, p, drop = FALSE] - theirs)), numeric(1))
  list(order = candidates[[which.min(gaps)]], gap = min(gaps))
}

#' Column order of a Mplus save file, from the output's SAVEDATA section
montecarlo_save_columns <- function(lines) {
  start <- grep("Order and format of variables", lines)[1L]
  stop <- grep("Save file format", lines)
  stop <- stop[stop > start][1L]
  stopifnot("no save-file layout in the output" = !is.na(start) && !is.na(stop))
  fields <- trimws(lines[(start + 1L):(stop - 1L)])
  fields <- fields[nzchar(fields)]
  toupper(vapply(strsplit(fields, "\\s+"), `[[`, character(1), 1L))
}

#' Everything compared from one Mplus run
montecarlo_mplus <- function(name, case, root = montecarlo_dir) {
  spec <- case$spec
  path <- file.path(root, "returned", paste0(name, ".out"))
  lines <- read_mplus_lines(path)
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(root, "returned", paste0(name, "-results.dat")),
                 quiet = TRUE)
  stopifnot("results file does not end with the fit block" =
              isTRUE(all.equal(values[2L * npar + 1L], npar)))
  saved <- utils::read.table(file.path(root, "returned",
                                       paste0(name, "-posteriors.dat")), na.strings = "*")
  columns <- montecarlo_save_columns(lines)
  stopifnot("save file does not match its layout" = length(columns) == ncol(saved))
  names(saved) <- columns
  two <- identical(spec$level, "two")
  posteriors <- if (two) {
    # Joint CB x CW probabilities, cb outer: column (h - 1) * K + k.
    joint <- as.matrix(saved[grep("^CPROB", columns)])
    k <- spec$k; h <- spec$h
    list(profiles = vapply(seq_len(k), \(j) rowSums(joint[, (seq_len(h) - 1L) * k + j,
                                                          drop = FALSE]), numeric(nrow(joint))),
         groups = vapply(seq_len(h), \(j) rowSums(joint[, (j - 1L) * k + seq_len(k),
                                                        drop = FALSE]), numeric(nrow(joint))))
  } else list(profiles = as.matrix(saved[grep("^CPROB", columns)]))
  classification <- if (two)
    mplus_classification_twolevel(path, spec$k, spec$h) else
      list(individuals = mplus_classification_single(path, spec$k))
  list(log_likelihood = values[2L * npar + 2L], n_parameters = npar,
       bic = values[2L * npar + 4L],
       entropy_precise = if (two) NA_real_ else values[2L * npar + 6L],
       n_observations = as.integer(printed_value(lines, "^Number of observations")),
       n_clusters = if (two) as.integer(printed_value(lines, "Number of clusters")) else NA_integer_,
       classification = classification,
       counts = if (two) mplus_variable_counts(lines, "CW", spec$k) else
         mplus_posterior_counts(lines, spec$k),
       id = as.integer(saved$ID), group = if (two) as.integer(saved[[ncol(saved)]]) else NULL,
       posteriors = posteriors,
       terminated = any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
       replicated = any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
}

#' Everything compared from one Latent GOLD run
montecarlo_latentgold <- function(name, case, root = montecarlo_dir) {
  spec <- case$spec
  lines <- iconv(readLines(file.path(root, "returned", paste0(name, ".lst")),
                           warn = FALSE), "latin1", "UTF-8")
  saved <- utils::read.delim(file.path(root, "returned",
                                       paste0(name, "_posteriors.txt")), check.names = FALSE)
  latent <- if (identical(spec$level, "transition")) "State" else "Cluster"
  list(log_likelihood = printed_value(lines, "^Log-likelihood \\(LL\\)"),
       n_parameters = as.integer(printed_value(lines, "^Number of parameters \\(Npar\\)")),
       id = as.integer(saved$id),
       profiles = as.matrix(saved[paste0(latent, "#", seq_len(spec$k))]),
       groups = if (identical(spec$level, "two"))
         as.matrix(saved[paste0("GClass#", seq_len(spec$h))]) else NULL,
       group = if (identical(spec$level, "two")) as.integer(saved[[2L]]) else NULL)
}

#' Profile means implied by a reference's posteriors (EM's fixed point)
montecarlo_implied_means <- function(posteriors, values) {
  observed <- !is.na(values)
  vapply(seq_len(ncol(values)), function(j) {
    weights <- posteriors * observed[, j]
    colSums(weights * ifelse(observed[, j], values[, j], 0)) / colSums(weights)
  }, numeric(ncol(posteriors)))
}

#' Fit every case on several cores
#'
#' Forked workers where they are safe; otherwise (macOS with Accelerate's
#' BLAS, which forked processes cannot use once the parent has) a socket
#' cluster whose workers load this source tree, so the tests always exercise
#' the code under test. A worker that returns nothing stops the run.
montecarlo_fit_all <- function(cases, n_starts = 20L) {
  # At most two workers: eight froze the shared development machine.
  cores <- if (.Platform$OS.type == "windows") 1L else
    max(1L, min(2L, parallel::detectCores() - 1L))
  worker <- function(case, n_starts) suppressWarnings(montecarlo_fit(case, n_starts = n_starts))
  if (cores == 1L) return(lapply(cases, worker, n_starts = n_starts))
  fits <- if (latents:::.latents_fork_safe()) {
    parallel::mclapply(cases, worker, n_starts = n_starts, mc.cores = cores,
                       mc.preschedule = FALSE)
  } else {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    setup <- function(root, helper) {
      suppressMessages(pkgload::load_all(root, quiet = TRUE))
      sys.source(helper, envir = globalenv())
      NULL
    }
    # Closures sent to a worker carry their environment; the global one is
    # empty there, so nothing of this call is serialized with them.
    environment(setup) <- globalenv()
    environment(worker) <- globalenv()
    parallel::clusterCall(cluster, setup, normalizePath(".."),
                          normalizePath("helper-montecarlo.R"))
    parallel::parLapply(cluster, cases, worker, n_starts = n_starts)
  }
  lost <- vapply(fits, \(fit) is.null(fit) || inherits(fit, "try-error"), logical(1))
  if (any(lost)) stop(sprintf("%d of %d fits were lost: %s", sum(lost), length(fits),
                              paste(names(fits)[lost], collapse = ", ")), call. = FALSE)
  fits
}
